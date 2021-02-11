/**
 * @file LemonTree.hpp
 * @brief Implementation of the class for learning module networks
 *        using the approach of Lemon Tree.
 * @author Ankit Srivastava <asrivast@gatech.edu>
 *
 * Copyright 2020 Georgia Institute of Technology
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
#ifndef DETAIL_LEMONTREE_HPP_
#define DETAIL_LEMONTREE_HPP_

#include "Ganesh.hpp"
#include "Module.hpp"
#include "PFCluster.hpp"

#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filter/gzip.hpp>


template <typename Data, typename Var, typename Set>
LemonTree<Data, Var, Set>::LemonTree(
  const mxx::comm& comm,
  const Data& data
) : ModuleNetworkLearning<Data, Var, Set>(comm, data)
{
  TIMER_RESET(m_tWrite);
  TIMER_RESET(m_tGanesh);
  TIMER_RESET(m_tConsensus);
  TIMER_RESET(m_tModules);
}

template <typename Data, typename Var, typename Set>
/**
 * @brief Default destructor.
 */
LemonTree<Data, Var, Set>::~LemonTree(
)
{
  if (this->m_comm.is_first()) {
    TIMER_ELAPSED_NONZERO("Time taken in writing the files: ", m_tWrite);
    TIMER_ELAPSED_NONZERO("Time taken in the GaneSH run: ", m_tGanesh);
    TIMER_ELAPSED_NONZERO("Time taken in consensus clustering: ", m_tConsensus);
    TIMER_ELAPSED_NONZERO("Time taken in learning the modules: ", m_tModules);
  }
}

template <typename Data, typename Var, typename Set>
std::list<std::list<Set>>
LemonTree<Data, Var, Set>::clusterVarsGanesh(
  const pt::ptree& ganeshConfigs
) const
{
  auto randomSeed = ganeshConfigs.get<uint32_t>("seed");
  std::mt19937 generator(randomSeed);
  auto initClusters = ganeshConfigs.get<Var>("init_num_clust");
  if ((initClusters == 0) || (initClusters > this->m_data.numVars())) {
    initClusters = this->m_data.numVars() / 2;
  }
  auto numRuns = ganeshConfigs.get<uint32_t>("num_runs");
  auto numSteps = ganeshConfigs.get<uint32_t>("num_steps");
  Ganesh<Data, Var, Set> ganesh(this->m_comm, this->m_data, &generator);
  std::list<std::list<Set>> sampledClusters;
  for (auto r = 0u; r < numRuns; ++r) {
    LOG_MESSAGE(info, "Run %u", r);
    ganesh.initializeRandom(initClusters);
    for (auto s = 0u; s <= numSteps; ++s) {
      LOG_MESSAGE(info, "Step %u", s);
      ganesh.clusterTwoWay();
    }
    // XXX: Lemon Tree writes out only the last sampled cluster
    // and uses that for the downstream tasks per run
    // We can replicate this behavior by storing only the last sample in each run
    // Further, we do not use -burn_in and -sample_steps because they do not
    // have any effect on the outcome
    LOG_MESSAGE(info, "Sampling");
    const auto& primaryClusters = ganesh.primaryClusters();
    std::list<Set> varClusters;
    for (const auto& cluster : primaryClusters) {
      varClusters.push_back(cluster.elements());
    }
    sampledClusters.push_back(varClusters);
  }
  return sampledClusters;
}

template <typename Data, typename Var, typename Set>
void
LemonTree<Data, Var, Set>::writeVarClusters(
  const std::string& clusterFile,
  const std::list<std::list<Set>>& varClusters
) const
{
  LOG_MESSAGE(info, "Writing variable clusters to %s", clusterFile);
  std::ofstream gf(clusterFile);
  auto g = 0u;
  for (auto git = varClusters.begin(); git != varClusters.end(); ++git, ++g) {
    std::string clustersFile(clusterFile + "." + std::to_string(g));
    std::ofstream cf(clustersFile);
    auto c = 0u;
    for (auto cit = git->begin(); cit != git->end(); ++cit, ++c) {
      for (const auto e : *cit) {
        cf << this->m_data.varName(e) << "\t" << c << std::endl;
      }
    }
    gf << clustersFile << std::endl;
  }
}

template <typename Data, typename Var, typename Set>
std::vector<double>
LemonTree<Data, Var, Set>::coclusteringMatrix(
  const std::list<std::list<Set>>&& sampledClusters,
  const double minWeight
) const
{
  std::list<std::unordered_map<Var, uint32_t>> varClusterMaps;
  for (const auto& varClusters : sampledClusters) {
    std::unordered_map<Var, uint32_t> thisMap;
    auto c = 0u;
    for (auto cit = varClusters.begin(); cit != varClusters.end(); ++cit, ++c) {
      for (const auto var : *cit) {
        thisMap[var] = c;
      }
    }
    varClusterMaps.push_back(thisMap);
  }
  const auto n = this->m_data.numVars();
  std::vector<double> coMatrix(n * n);
  for (auto u = 0u; u < n; ++u) {
    for (auto v = u + 1; v < n; ++v) {
      auto cooccurrence = 0u;
      for (const auto& varCluster : varClusterMaps) {
        if (varCluster.at(u) == varCluster.at(v)) {
          ++cooccurrence;
        }
      }
      auto weight = static_cast<double>(cooccurrence) / sampledClusters.size();
      if (weight > minWeight) {
        coMatrix[u * n + v] = weight;
        coMatrix[v * n + u] = weight;
      }
    }
  }
  return coMatrix;
}

template <typename Data, typename Var, typename Set>
std::multimap<Var, Var>
LemonTree<Data, Var, Set>::clusterConsensus(
  const std::list<std::list<Set>>&& varClusters,
  const pt::ptree& consensusConfigs
) const
{
  auto coMatrix = this->coclusteringMatrix(std::move(varClusters), consensusConfigs.get<double>("min_weight"));
  auto tolerance = consensusConfigs.get<double>("tolerance");
  auto maxSteps = consensusConfigs.get<uint32_t>("max_steps");
  auto minClustSize = consensusConfigs.get<uint32_t>("min_clust_size");
  auto minClustScore = consensusConfigs.get<double>("min_clust_score");
  auto result = densePerron(std::move(coMatrix), this->m_data.numVars(), tolerance,
                            maxSteps, minClustSize, minClustScore);
  return result;
}

template <typename Data, typename Var, typename Set>
void
LemonTree<Data, Var, Set>::writeConsensusCluster(
  const std::string& consensusFile,
  const std::multimap<Var, Var>& vertexClusters
) const
{
  LOG_MESSAGE(info, "Writing consensus clusters to %s", consensusFile);
  std::ofstream out(consensusFile);
  // Write output in tab seperated format
  for (const auto& cx : vertexClusters) {
    out << this->m_data.varName(cx.second) << "\t" << static_cast<uint32_t>(cx.first) << std::endl;
  }
}

template <typename Data, typename Var, typename Set>
std::list<std::list<Set>>
LemonTree<Data, Var, Set>::clusterObsGanesh(
  const uint32_t numRuns,
  const uint32_t numSteps,
  const uint32_t burnSteps,
  const uint32_t sampleSteps,
  std::mt19937* const generator,
  const Set& clusterVars
) const
{
  // Initialize Gibbs sampler algorithm for this cluster
  Ganesh<Data, Var, Set> ganesh(this->m_comm, this->m_data, generator);
  ganesh.initializeGiven(std::list<Set>(1, clusterVars));
  std::list<std::list<Set>> sampledClusters;
  for (auto r = 0u; r < numRuns; ++r) {
    auto s = 0u;
    for ( ; s < burnSteps; ++s) {
      LOG_MESSAGE(info, "Step %u (burn in)", s);
      ganesh.clusterSecondary();
    }
    for ( ; s < numSteps; ++s) {
      LOG_MESSAGE(info, "Step %u (sampling)", s);
      ganesh.clusterSecondary();
      if ((numSteps - (s + 1)) % sampleSteps == 0) {
        LOG_MESSAGE(info, "Sampling");
        // There should be only one primary cluster; get a reference to it
        const auto& primaryCluster = ganesh.primaryClusters().front();
        const auto& secondaryClusters = primaryCluster.secondaryClusters();
        std::list<Set> obsClusters;
        for (const auto& cluster : secondaryClusters) {
          obsClusters.push_back(cluster.elements());
        }
        sampledClusters.push_back(obsClusters);
      }
    }
  }
  return sampledClusters;
}

template <typename Data, typename Var, typename Set>
void
LemonTree<Data, Var, Set>::readCandidateParents(
  const std::string& fileName,
  Set& candidateParents
) const
{
  LOG_MESSAGE(info, "Reading candidate parents from %s", fileName);
  std::ifstream regFile(boost::filesystem::canonical(fileName).string());
  std::string name;
  while (std::getline(regFile, name)) {
    auto v = this->m_data.varIndex(name);
    if (v < this->m_data.numVars()) {
      candidateParents.insert(v);
    }
  }
  LOG_MESSAGE(info, "Read %u candidate parents", candidateParents.size());
}

template <typename Data, typename Var, typename Set>
std::list<Module<Data, Var, Set>>
LemonTree<Data, Var, Set>::learnModules(
  const std::multimap<Var, Var>&& coClusters,
  const pt::ptree& modulesConfigs
) const
{
  auto randomSeed = modulesConfigs.get<uint32_t>("seed");
  auto numRuns = modulesConfigs.get<uint32_t>("num_runs");
  auto numSteps = modulesConfigs.get<uint32_t>("num_steps");
  auto burnSteps = modulesConfigs.get<uint32_t>("burn_in");
  auto sampleSteps = modulesConfigs.get<uint32_t>("sample_steps");
  auto scoreBHC = modulesConfigs.get<bool>("use_bayesian_score");
  auto scoreGain = modulesConfigs.get<double>("score_gain");
  auto regFile = modulesConfigs.get<std::string>("reg_file");
  auto betaMax = modulesConfigs.get<double>("beta_reg");
  auto numSplits = modulesConfigs.get<uint32_t>("num_reg");
  std::mt19937 generator(randomSeed);
  std::list<Module<Data, Var, Set>> modules;
  for (auto cit = coClusters.begin(); cit != coClusters.end(); ) {
    // Get the range of variables in this cluster
    auto clusterIts = coClusters.equal_range(cit->first);
    // Add all the variables in the cluster to a set
    auto clusterVars = Set(this->m_data.numVars());
    for (auto vit = clusterIts.first; vit != clusterIts.second; ++vit) {
      clusterVars.insert(vit->second);
    }
    // Create a module for this variable cluster
    modules.emplace_back(std::move(clusterVars), this->m_comm, this->m_data);
    auto& module = modules.back();
    // Sample observation clusters for this module
    auto sampledClusters = this->clusterObsGanesh(numRuns, numSteps, burnSteps, sampleSteps,
                                                  &generator, module.variables());
    // Learn tree structures from the observation clusters
    module.learnTreeStructures(std::move(sampledClusters), scoreBHC, scoreGain);
    cit = clusterIts.second;
  }
  Set candidateParents(this->m_data.numVars());
  if (!regFile.empty()) {
    // Read candidate parents from the given file
    this->readCandidateParents(regFile, candidateParents);
  }
  else {
    // Add all the variables as candidate parents
    for (Var v = 0u; v < candidateParents.max(); ++v) {
      candidateParents.insert(v);
    }
  }
  OptimalBeta ob(0.0, betaMax, 1e-5);
  for (auto& module : modules) {
    module.learnParents(generator, candidateParents, ob, numSplits);
  }
  return modules;
}

template <typename Data, typename Var, typename Set>
void
LemonTree<Data, Var, Set>::writeParents(
  std::ofstream& stream,
  const std::unordered_map<Var, double>& splits,
  const uint32_t moduleIndex,
  const double cutoff
) const
{
  std::set<std::pair<Var, double>, std::greater<std::pair<Var, double>>> splitsSorted(splits.begin(), splits.end());
  for (const auto& split : splitsSorted) {
    if (split.second > cutoff) {
      stream << this->m_data.varName(split.first) << "\t" << moduleIndex << "\t" << split.second << std::endl;
    }
  }
}

template <typename Data, typename Var, typename Set>
void
LemonTree<Data, Var, Set>::writeModules(
  const std::string& modulesFile,
  const std::list<Module<Data, Var, Set>>& modules
) const
{
  std::string allParentsFile = modulesFile + ".allreg.txt";
  LOG_MESSAGE(info, "Writing all parents to %s", allParentsFile);
  std::ofstream apf(allParentsFile);
  std::string topParentsFile = modulesFile + ".topreg.txt";
  LOG_MESSAGE(info, "Writing top 1%% parents to %s", topParentsFile);
  std::ofstream tpf(topParentsFile);
  std::string randParentsFile = modulesFile + ".randomreg.txt";
  LOG_MESSAGE(info, "Writing random parents to %s", randParentsFile);
  std::ofstream rpf(randParentsFile);
  std::string xmlFile = modulesFile + ".xml.gz";
  LOG_MESSAGE(info, "Writing modules to XML file %s", xmlFile);
  std::ofstream file(xmlFile, std::ios_base::out | std::ios_base::binary);
  boost::iostreams::filtering_streambuf<boost::iostreams::output> out;
  out.push(boost::iostreams::gzip_compressor());
  out.push(file);
  std::ostream xmlf(&out);

  // XXX: Try to match JAVA precision
  xmlf.precision(16);
  xmlf << "<?xml version='1.0' encoding='iso-8859-1'?>" << std::endl;
  xmlf << "<ModuleNetwork>" << std::endl;
  // First compute the cutoff for top parents
  std::multiset<double, std::greater<double>> allScores;
  for (const auto& module : modules) {
    for (const auto& parent : module.allParents()) {
      allScores.insert(parent.second);
    }
  }
  auto index = static_cast<uint32_t>(round(0.01 * allScores.size()));
  auto cutoff = *std::next(allScores.begin(), index);
  // Now, write the modules
  auto m = 0u;
  for (const auto& module : modules) {
    module.toXML(xmlf, m);
    this->writeParents(apf, module.allParents(), m);
    this->writeParents(tpf, module.allParents(), m, cutoff);
    this->writeParents(rpf, module.randParents(), m);
    ++m;
  }
  xmlf << "</ModuleNetwork>" << std::endl;
}

template <typename Data, typename Var, typename Set>
ModuleNetwork<Var>
LemonTree<Data, Var, Set>::getNetwork_sequential(
  const pt::ptree& algoConfigs
) const
{
  TIMER_START(m_tGanesh);
  const auto& ganeshConfigs = algoConfigs.get_child("ganesh");
  auto varClusters = this->clusterVarsGanesh(ganeshConfigs);
  TIMER_PAUSE(m_tGanesh);
  auto clusterFile = ganeshConfigs.get<std::string>("output_file");
  if (!clusterFile.empty()) {
    TIMER_START(m_tWrite);
    this->writeVarClusters(clusterFile, varClusters);
    TIMER_PAUSE(m_tWrite);
  }
  TIMER_START(m_tConsensus);
  const auto& consensusConfigs = algoConfigs.get_child("tight_clusters");
  auto coClusters = this->clusterConsensus(std::move(varClusters), consensusConfigs);
  TIMER_PAUSE(m_tConsensus);
  auto consensusFile = consensusConfigs.get<std::string>("output_file");
  if (!consensusFile.empty()) {
    TIMER_START(m_tWrite);
    this->writeConsensusCluster(consensusFile, coClusters);
    TIMER_PAUSE(m_tWrite);
  }
  TIMER_START(m_tModules);
  const auto& modulesConfigs = algoConfigs.get_child("regulators");
  auto modules = this->learnModules(std::move(coClusters), modulesConfigs);
  TIMER_PAUSE(m_tModules);
  auto modulesFile = modulesConfigs.get<std::string>("output_file");
  if (!modulesFile.empty()) {
    TIMER_START(m_tWrite);
    this->writeModules(modulesFile, modules);
    TIMER_PAUSE(m_tWrite);
  }
  return ModuleNetwork<Var>(this->m_data.varNames(this->m_allVars));
}

template <typename Data, typename Var, typename Set>
ModuleNetwork<Var>
LemonTree<Data, Var, Set>::getNetwork_parallel(
  const pt::ptree&
) const
{
  throw NotImplementedError("Lemon Tree: Parallel algorithm is not implemented yet");
  return ModuleNetwork<Var>(this->m_data.varNames(this->m_allVars));
}

#endif // DETAIL_LEMONTREE_HPP_
