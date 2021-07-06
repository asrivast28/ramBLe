/**
 * @file GlobalLearning.hpp
 * @brief Implementation of the classes for global learning algorithms.
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
#ifndef DETAIL_GLOBALLEARNING_HPP_
#define DETAIL_GLOBALLEARNING_HPP_

#include "common/SetUtils.hpp"
#include "utils/Logging.hpp"

#include <boost/math/special_functions/binomial.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/iterator/zip_iterator.hpp>
#include <boost/tuple/tuple.hpp>


template <typename Data, typename Var, typename Set>
/**
 * @brief Constructs the object with the given data.
 */
GlobalLearning<Data, Var, Set>::GlobalLearning(
  const mxx::comm& comm,
  const Data& data,
  const double alpha,
  const Var maxConditioning
) : ConstraintBasedLearning<Data, Var, Set>(comm, data, alpha, maxConditioning),
    m_cachedNeighbors(),
    m_removedEdges()
{
  TIMER_RESET(m_tSync);
  TIMER_RESET(m_tDist);
}

template <typename Data, typename Var, typename Set>
/**
 * @brief Default destructor that prints out timing information.
 */
GlobalLearning<Data, Var, Set>::~GlobalLearning(
)
{
  if (this->m_comm.is_first()) {
    TIMER_ELAPSED_NONZERO("Time taken in redistributing: ", m_tDist);
    TIMER_ELAPSED_NONZERO("Time taken in synchronizing neighbors: ", m_tSync);
  }
}

template <typename Data, typename Var, typename Set>
const Set&
GlobalLearning<Data, Var, Set>::getPC(
  const Var target
) const
{
  return m_cachedNeighbors.at(target);
}

template <typename Data, typename Var, typename Set>
const Set&
GlobalLearning<Data, Var, Set>::getMB(
  const Var target
) const
{
  return m_cachedNeighbors.at(target);
}

template <typename Data, typename Var, typename Set>
/**
 * @brief Function for initializing the data structures used for
 *        learning the skeleton.
 *
 * @param myEdges The array to be initialized with the edges
 *                to be handled by this processor.
 * @param allNeighbors The map containing the neighborhood set
 *                     for all the variables in the network.
 * @param allNeighbors The map used for tracking the removals for
 *                     for all the variables in the network.
 */
void
GlobalLearning<Data, Var, Set>::initializeLearning(
  std::vector<std::tuple<Var, Var, double>>& myEdges,
  std::unordered_map<Var, Set>& allNeighbors,
  std::unordered_map<Var, Set>& removedNeighbors
) const
{
  // First, block decompose all the edges on all the processors
  auto n = this->m_allVars.size();
  mxx::blk_dist dist((n * (n - 1)) / 2, this->m_comm);
  auto myOffset = dist.eprefix_size();
  auto mySize = dist.local_size();

  auto vars = std::vector<Var>(this->m_allVars.begin(), this->m_allVars.end());
  auto k = 0u;
  auto currOffset = 0u;
  while (currOffset + (n - (k + 1)) <= static_cast<uint32_t>(myOffset)) {
    currOffset += (n - (k + 1));
    ++k;
  }
  auto primary = vars.begin() + k;
  auto secondary = primary + 1 + (myOffset - currOffset);

  myEdges.resize(mySize);
  for (auto i = 0u; i < mySize; ++i, ++secondary) {
    if (secondary == vars.end()) {
      // Start a new primary variable if the secondary variable is exhausted
      ++primary;
      secondary = primary + 1;
    }
    // Initialize the p-values
    myEdges[i] = std::make_tuple(*primary, *secondary, 0.0);
  }
  // Then, initialize all the neighbor sets
  for (const auto v : vars) {
    allNeighbors.insert(std::make_pair(v, set_init(this->getCandidates(v), this->m_data.numVars())));
    removedNeighbors.insert(std::make_pair(v, set_init(Set(), this->m_data.numVars())));
  }
}

template <typename Data, typename Var, typename Set>
bool
GlobalLearning<Data, Var, Set>::fixWeightedImbalance(
  std::vector<std::tuple<Var, Var, double>>& myEdges,
  const std::vector<double>& myWeights,
  const double imbalanceThreshold
) const
{
  bool fixed = false;
  double myTotalWeight = std::accumulate(myWeights.cbegin(), myWeights.cend(), 0.0);
  TIMER_START(this->m_tMxx);
  double totalWeight = mxx::allreduce(myTotalWeight, this->m_comm);
  TIMER_PAUSE(this->m_tMxx);
  auto imbalance = 0.0;
  if (totalWeight > 0) {
    TIMER_START(this->m_tMxx);
    auto maxWeight = mxx::allreduce(myTotalWeight, mxx::max<double>(), this->m_comm);
    TIMER_PAUSE(this->m_tMxx);
    auto avgWeight = totalWeight / this->m_comm.size();
    imbalance = (maxWeight / avgWeight) - 1.0;
    //if (this->m_comm.is_first()) {
      //std::cout << "Imbalance: " << imbalance << std::endl;
    //}
  }
  if (std::isgreater(imbalance, imbalanceThreshold)) {
    TIMER_START(this->m_tMxx);
    // Get the weight on previous processors
    auto globalPrefix = mxx::exscan(myTotalWeight, this->m_comm);
    TIMER_PAUSE(this->m_tMxx);
    std::vector<uint64_t> sendCounts(this->m_comm.size(), 0);
    if (std::isgreater(myTotalWeight, 0)) {
      double div = totalWeight / this->m_comm.size();
      auto proc = static_cast<uint32_t>(std::floor(globalPrefix / div));
      auto procFirst = 0u;
      double localPrefix = 0;
      double leftWeight = myTotalWeight;
      std::vector<double> myWeightsPrefix(myWeights.size());
      std::partial_sum(myWeights.cbegin(), myWeights.cend(), myWeightsPrefix.begin());
      for (; std::isgreater(leftWeight, 0) && (proc < this->m_comm.size()); ++proc) {
        double sendWeight = std::min<double>((div * (proc + 1)) - globalPrefix, leftWeight);
        auto procLast = std::distance(myWeightsPrefix.cbegin(),
                                      std::lower_bound(myWeightsPrefix.cbegin(),
                                                       myWeightsPrefix.cend(),
                                                       localPrefix + sendWeight));

        if (std::isless(myWeightsPrefix[procLast], localPrefix + sendWeight)) {
          ++procLast;
        }
        sendCounts[proc] = (procLast - procFirst) + 1;
        leftWeight -= (myWeightsPrefix[procLast] - localPrefix);
        globalPrefix += (myWeightsPrefix[procLast] - localPrefix);
        localPrefix = myWeightsPrefix[procLast];
        procFirst = procLast + 1;
      }
    }
    TIMER_START(this->m_tMxx);
    auto recvCounts = mxx::all2all(sendCounts, this->m_comm);
    TIMER_PAUSE(this->m_tMxx);
    auto myRecv = std::accumulate(recvCounts.begin(), recvCounts.end(), 0u);
    std::vector<std::tuple<Var, Var, double>> newEdges(myRecv);
    const std::tuple<Var, Var, double>* sendBuf = (myEdges.size() > 0) ? &myEdges[0] : nullptr;
    std::tuple<Var, Var, double>* recvBuf = (newEdges.size() > 0) ? &newEdges[0] : nullptr;
    TIMER_START(this->m_tMxx);
    mxx::all2allv(sendBuf, sendCounts, recvBuf, recvCounts, this->m_comm);
    TIMER_PAUSE(this->m_tMxx);
    myEdges = std::move(newEdges);
    fixed = true;
  }
  return fixed;
}

template <typename Data, typename Var, typename Set>
/**
 * @brief Function for storing removed edges if they might be relevant
 *        for directing edges later.
 *
 * @param myRemoved A vector of all the removed edges.
 * @param myDSepSets Vector containing d-separating sets for all the removed edges.
 * @param allNeighbors Neighbor sets for all the variables.
 */
void
GlobalLearning<Data, Var, Set>::storeRemovedEdges(
  std::vector<std::tuple<Var, Var, double>>&& myRemoved,
  std::vector<Set>&& myDSepSets,
  const std::unordered_map<Var, Set>& allNeighbors,
  const bool duplicateEdges
) const
{
  LOG_MESSAGE_IF(myRemoved.size() != myDSepSets.size(),
                 error, "Mismatch between number of edges and d-separating sets.");
  // Only retain information which can be used for directing edges later
  auto newEnd = std::remove_if(boost::make_zip_iterator(boost::make_tuple(myRemoved.begin(), myDSepSets.begin())),
                               boost::make_zip_iterator(boost::make_tuple(myRemoved.end(), myDSepSets.end())),
                               [&allNeighbors] (const boost::tuple<std::tuple<Var, Var, double>, Set>& e)
                                               { return set_intersection(allNeighbors.at(std::get<0>(boost::get<0>(e))),
                                                                         allNeighbors.at(std::get<1>(boost::get<0>(e)))).empty(); });
  myRemoved.erase(boost::get<0>(newEnd.get_iterator_tuple()), myRemoved.end());
  myDSepSets.erase(boost::get<1>(newEnd.get_iterator_tuple()), myDSepSets.end());
  TIMER_START(this->m_tMxx);
  // Get the prefix size of all the d-separating sets
  auto removedSizes = mxx::allgather(myRemoved.size(), this->m_comm);
  TIMER_PAUSE(this->m_tMxx);
  auto removedPrefix = std::accumulate(removedSizes.begin(), removedSizes.begin() + this->m_comm.rank(), 0u);
  std::vector<std::tuple<Var, Var, double, uint32_t>> myRemovedEdges(myRemoved.size());
  std::transform(std::make_move_iterator(myRemoved.begin()), std::make_move_iterator(myRemoved.end()),
                 boost::counting_iterator<uint32_t>(removedPrefix), myRemovedEdges.begin(),
                 [] (const std::tuple<Var, Var, double>&& e, const uint32_t idx)
                    { return std::make_tuple(std::get<0>(e), std::get<1>(e), std::get<2>(e), idx); });
  TIMER_START(this->m_tMxx);
  // Block redistribute all the removed edges
  mxx::stable_distribute_inplace(myRemovedEdges, this->m_comm);
  // Sort the removed edges in parallel
  mxx::comm nonzero_comm(static_cast<MPI_Comm>(this->m_comm));
  if (mxx::any_of(myRemovedEdges.size() == 0, this->m_comm)) {
    nonzero_comm = this->m_comm.split(myRemovedEdges.size() > 0);
  }
  if (myRemovedEdges.size() > 0) {
    auto sortEdges = [] (const std::tuple<Var, Var, double, uint32_t>& a, const std::tuple<Var, Var, double, uint32_t>& b)
                        { return std::make_pair(std::get<0>(a), std::get<1>(a)) < std::make_pair(std::get<0>(b), std::get<1>(b)); };
    if (!duplicateEdges) {
      mxx::sort(myRemovedEdges.begin(), myRemovedEdges.end(), sortEdges, nonzero_comm);
    }
    else {
      mxx::stable_sort(myRemovedEdges.begin(), myRemovedEdges.end(), sortEdges, nonzero_comm);
      auto newEnd = mxx::unique(myRemovedEdges.begin(), myRemovedEdges.end(),
                                [] (const std::tuple<Var, Var, double, uint32_t>& a, const std::tuple<Var, Var, double, uint32_t>& b)
                                   { return std::make_pair(std::get<0>(a), std::get<1>(a)) == std::make_pair(std::get<0>(b), std::get<1>(b)); },
                                nonzero_comm);
      myRemovedEdges.erase(newEnd, myRemovedEdges.end());
    }
  }
  // Gather all the removed edges on all the processors
  auto allRemovedEdges = mxx::allgatherv(std::move(myRemovedEdges), this->m_comm);
  // Also gather all the d-separating sets on all the processors
  // XXX: We can not use Set as part of a std::tuple during communication
  //      Therefore, we gather all the sets separately
  auto allDSepSets = set_allgatherv(std::move(myDSepSets), removedSizes, this->m_data.numVars(), this->m_comm);
  TIMER_PAUSE(this->m_tMxx);
  // Finally, store the removed edges with the corresponding d-separating sets
  m_removedEdges.resize(allRemovedEdges.size());
  std::transform(std::make_move_iterator(allRemovedEdges.begin()), std::make_move_iterator(allRemovedEdges.end()),
                 m_removedEdges.begin(),
                 [&allDSepSets] (const std::tuple<Var, Var, double, uint32_t>&& e)
                                { return std::make_tuple(std::get<0>(e), std::get<1>(e), std::get<2>(e), allDSepSets.at(std::get<3>(e))); });
}

template <typename Data, typename Var, typename Set>
/**
 * @brief Function for constructing skeleton from the neighbor sets.
 *
 * @param allNeighbors Neighborhood sets for all the variables.
 */
BayesianNetwork<Var>
GlobalLearning<Data, Var, Set>::constructSkeleton(
  const std::unordered_map<Var, Set>&& allNeighbors
) const
{
  BayesianNetwork<Var> bn(this->m_data.varNames(this->m_allVars));
  for (const auto& neighbors : allNeighbors) {
    const auto x = neighbors.first;
    for (const auto y : neighbors.second) {
      if (x < y) {
        LOG_MESSAGE(info, "+ Adding the edge %s <-> %s",
                          this->m_data.varName(x), this->m_data.varName(y));
        bn.addEdge(x, y, true);
      }
    }
  }
  m_cachedNeighbors = std::move(allNeighbors);
  return bn;
}

template <typename Data, typename Var, typename Set>
/**
 * @brief Checks if the given v-structure (y-x-z) is a collider.
 *        Also returns the p-value for the collider.
 *
 * This function checks if x is part of the set that rendered
 * y and z independent. If not, y-x-z is a collider.
 */
std::pair<bool, double>
GlobalLearning<Data, Var, Set>::checkCollider(
  const Var y,
  const Var x,
  const Var z
) const
{
  static auto findPair = [] (const std::tuple<Var, Var, double, Set>& t, const std::pair<Var, Var>& e)
                            { return std::make_pair(std::get<0>(t), std::get<1>(t)) < e; };
  // Always try to find the forward edge
  auto edge = (y < z) ? std::make_pair(y, z) : std::make_pair(z, y);
  auto eIt = std::lower_bound(m_removedEdges.begin(), m_removedEdges.end(), edge, findPair);
  if ((eIt == m_removedEdges.end()) ||
      (std::make_pair(std::get<0>(*eIt), std::get<1>(*eIt)) != edge)) {
    auto pv = this->m_data.pValue(y, z);
    LOG_MESSAGE(debug, "Computed p-value for edge %s - %s is %g",
                       this->m_data.varName(y), this->m_data.varName(z), pv);
    return std::make_pair(true, pv);
  }
  else {
    auto collider = !(std::get<3>(*eIt).contains(x));
    auto pv = std::get<2>(*eIt);
    LOG_MESSAGE(debug, "Stored p-value for edge %s - %s is %g",
                       this->m_data.varName(y), this->m_data.varName(z), std::get<2>(*eIt));
    return std::make_pair(collider, pv);
  }
}

template <typename Data, typename Var, typename Set>
PCStableCommon<Data, Var, Set>::PCStableCommon(
  const mxx::comm& comm,
  const Data& data,
  const double alpha,
  const Var maxConditioning
) : GlobalLearning<Data, Var, Set>(comm, data, alpha, maxConditioning)
{
}

template <typename Data, typename Var, typename Set>
std::pair<double, Set>
PCStableCommon<Data, Var, Set>::checkEdge(
  const std::tuple<Var, Var, double>& edge,
  const std::unordered_map<Var, Set>& allNeighbors,
  std::unordered_map<Var, Set>& removedNeighbors,
  const uint32_t setSize,
  const bool checkBackward
) const
{
  auto pv = 0.0;
  auto dsep = set_init(Set(), this->m_data.numVars());
  const auto x = std::get<0>(edge);
  const auto y = std::get<1>(edge);
  auto xNeighbors = allNeighbors.at(x);
  auto yNeighbors = allNeighbors.at(y);
  LOG_MESSAGE(debug, "Checking the edge %s <-> %s, d-separating set of size %u",
                     this->m_data.varName(x), this->m_data.varName(y), setSize);
  auto remove = false;
  // First, check the edge using neighbors of x
  if (xNeighbors.size() > setSize) {
    xNeighbors.erase(y);
    std::tie(pv, dsep) = this->m_data.maxPValueSubset(this->m_alpha, x, y, xNeighbors, setSize, setSize);
    remove = this->m_data.isIndependent(this->m_alpha, pv);
  }
  // Then, check the edge using neighbors of y if all of the following hold:
  // 1. The conditioning set used for checking is not empty
  // 2. The neighbors of x did not remove it already
  // 3. The size of the neighborhood of y is greater than the conditioning set size
  if (checkBackward && !remove && (yNeighbors.size() > setSize)) {
    yNeighbors.erase(x);
    // Further, only check if neighborhood of y has some elements which are
    // not present in the neighborhood of x
    if (!set_difference(yNeighbors, xNeighbors).empty()) {
      std::tie(pv, dsep) = this->m_data.maxPValueSubset(this->m_alpha, x, y, yNeighbors, setSize, setSize);
      remove = this->m_data.isIndependent(this->m_alpha, pv);
    }
  }
  LOG_MESSAGE(debug, "%s and %s are " + std::string(this->m_data.isIndependent(this->m_alpha, pv) ? "independent" : "dependent") +
                     " (p-value = %g)",
                     this->m_data.varName(x), this->m_data.varName(y), pv);
  LOG_MESSAGE_IF(remove, debug, "- Removing the edge %s <-> %s",
                                this->m_data.varName(x), this->m_data.varName(y));
  if (remove) {
    // Mark y for removal from neighborhood of x
    removedNeighbors.at(x).insert(y);
    // Mark x for removal from neighborhood of y
    removedNeighbors.at(y).insert(x);
  }
  return std::make_pair(pv, dsep);
}

template <typename Data, typename Var, typename Set>
BayesianNetwork<Var>
PCStableCommon<Data, Var, Set>::getSkeleton_sequential(
  const bool directEdges
) const
{
  std::vector<std::tuple<Var, Var, double>> allEdges;
  std::unordered_map<Var, Set> allNeighbors;
  std::unordered_map<Var, Set> removedNeighbors;
  this->initializeLearning(allEdges, allNeighbors, removedNeighbors);
  auto maxSize = std::min(this->m_maxConditioning, static_cast<Var>(this->m_allVars.size() - 2));
  for (auto s = 0u; (s <= maxSize) && !allEdges.empty(); ++s) {
    LOG_MESSAGE(debug, "Testing %u edges using sets of size %u", allEdges.size(), s);
    TIMER_DECLARE(tIter);
    for (auto& e : allEdges) {
      auto result = this->checkEdge(e, allNeighbors, removedNeighbors, s, (s > 0));
      std::get<2>(e) = result.first;
      // We need to store the removed edges since the d-separating
      // set may be required for directing edges later
      if (directEdges && (s > 0) &&
          this->m_data.isIndependent(this->m_alpha, result.first)) {
        TIMER_START(this->m_tDirect);
        Var x, y;
        std::tie(x, y, std::ignore) = e;
        // We only use this information for directing colliders of type x-z-y
        // However, there is no need to store it if there are no candidates for z
        // i.e., no common neighbors between x and y
        if (!set_intersection(allNeighbors.at(x), allNeighbors.at(y)).empty()) {
          this->m_removedEdges.push_back(std::make_tuple(x, y, result.first, result.second));
        }
        TIMER_PAUSE(this->m_tDirect);
      }
    }
    for (auto& rn : removedNeighbors) {
      if (!rn.second.empty()) {
        allNeighbors[rn.first] = set_difference(allNeighbors.at(rn.first), rn.second);
        rn.second.clear();
      }
    }
    auto newEnd = std::remove_if(allEdges.begin(), allEdges.end(),
                                 [this, &allNeighbors, &s] (const std::tuple<Var, Var, double>& e)
                                                           { return this->m_data.isIndependent(this->m_alpha, std::get<2>(e)) ||
                                                                    (allNeighbors.at(std::get<0>(e)).size() <= (s + 1) &&
                                                                     allNeighbors.at(std::get<1>(e)).size() <= (s + 1)); });
    allEdges.erase(newEnd, allEdges.end());
    if (this->m_comm.is_first()) {
      TIMER_ELAPSED("Time taken in testing all sets of size " + std::to_string(s) + ": ", tIter);
    }
    if (directEdges) {
      TIMER_START(this->m_tDirect);
      auto newEnd = std::remove_if(this->m_removedEdges.begin(), this->m_removedEdges.end(),
                                   [&allNeighbors] (const std::tuple<Var, Var, double, Set>& t)
                                                   { return set_intersection(allNeighbors.at(std::get<0>(t)),
                                                                             allNeighbors.at(std::get<1>(t))).empty(); });
      this->m_removedEdges.erase(newEnd, this->m_removedEdges.end());
      TIMER_PAUSE(this->m_tDirect);
    }
  }
  if (directEdges) {
    TIMER_START(this->m_tDirect);
    std::sort(this->m_removedEdges.begin(), this->m_removedEdges.end(),
              [] (const std::tuple<Var, Var, double, Set>& a, const std::tuple<Var, Var, double, Set>& b)
                 { return std::make_pair(std::get<0>(a), std::get<1>(a)) < std::make_pair(std::get<0>(b), std::get<1>(b)); });
    TIMER_PAUSE(this->m_tDirect);
  }
  return this->constructSkeleton(std::move(allNeighbors));
}

template <typename Data, typename Var, typename Set>
PCStable<Data, Var, Set>::PCStable(
  const mxx::comm& comm,
  const Data& data,
  const double alpha,
  const Var maxConditioning
) : PCStableCommon<Data, Var, Set>(comm, data, alpha, maxConditioning)
{
}

template <typename Data, typename Var, typename Set>
/**
 * @brief Function that synchronizes the candidate sets by taking an
 *        intersection of the sets across all the processors.
 *
 * @param mySets A map with the neighbors of the variables
 *               remaining on this processor.
 */
void
PCStable<Data, Var, Set>::syncSets(
  std::unordered_map<Var, Set>& mySets
) const
{
  TIMER_START(this->m_tMxx);
  set_allintersect_indexed(mySets, this->m_allVars, this->m_data.numVars(), this->m_comm);
  TIMER_PAUSE(this->m_tMxx);
}

template <typename Data, typename Var, typename Set>
BayesianNetwork<Var>
PCStable<Data, Var, Set>::getSkeleton_parallel(
  const bool directEdges,
  const double imbalanceThreshold
) const
{
  // Initialize the learning in a similar fashion as done sequentially
  // Create nC2 edges - one for each unordered variable pair
  std::vector<std::tuple<Var, Var, double>> myEdges;
  std::unordered_map<Var, Set> allNeighbors;
  std::unordered_map<Var, Set> removedNeighbors;
  this->initializeLearning(myEdges, allNeighbors, removedNeighbors);
  // The following two data structures are used for getting
  // d-separating set and p-value for all the removed edges
  std::vector<std::tuple<Var, Var, double>> myRemoved;
  std::vector<Set> myDSepSets;
  auto maxSize = std::min(this->m_maxConditioning, static_cast<Var>(this->m_allVars.size() - 2));
  for (auto s = 0u; (s <= maxSize) && mxx::any_of(myEdges.size() > 0, this->m_comm); ++s) {
    LOG_MESSAGE(debug, "Testing %u edges using sets of size %u", myEdges.size(), s);
    TIMER_DECLARE(tIter);
    for (auto& e : myEdges) {
      auto result = this->checkEdge(e, allNeighbors, removedNeighbors, s, (s > 0));
      std::get<2>(e) = result.first;
      // We need to store the removed edges since the d-separating
      // set may be required for directing edges later
      if (directEdges && (s > 0) &&
          this->m_data.isIndependent(this->m_alpha, result.first)) {
        TIMER_START(this->m_tDirect);
        Var x, y;
        std::tie(x, y, std::ignore) = e;
        // We only use this information for directing colliders of type x-z-y
        // However, there is no need to store it if there are no candidates for z
        // i.e., no common neighbors between x and y
        if (!set_intersection(allNeighbors.at(x), allNeighbors.at(y)).empty()) {
          myRemoved.push_back(e);
          myDSepSets.push_back(result.second);
        }
        TIMER_PAUSE(this->m_tDirect);
      }
    }
    auto newEnd = std::remove_if(myEdges.begin(), myEdges.end(),
                                 [this] (const std::tuple<Var, Var, double>& e)
                                        { return this->m_data.isIndependent(this->m_alpha, std::get<2>(e)); });
    myEdges.erase(newEnd, myEdges.end());
    for (auto& rn : removedNeighbors) {
      if (!rn.second.empty()) {
        allNeighbors[rn.first] = set_difference(allNeighbors.at(rn.first), rn.second);
        rn.second.clear();
      }
    }
    this->m_comm.barrier();
    TIMER_START(this->m_tSync);
    this->syncSets(allNeighbors);
    TIMER_PAUSE(this->m_tSync);
    newEnd = std::remove_if(myEdges.begin(), myEdges.end(),
                            [&allNeighbors, &s] (const std::tuple<Var, Var, double>& e)
                                                { return allNeighbors.at(std::get<0>(e)).size() <= (s + 1) &&
                                                         allNeighbors.at(std::get<1>(e)).size() <= (s + 1); });
    myEdges.erase(newEnd, myEdges.end());
    if (this->m_comm.is_first()) {
      TIMER_ELAPSED("Time taken in testing all sets of size " + std::to_string(s) + ": ", tIter);
    }
    if (std::isgreaterequal(imbalanceThreshold, 0.0)) {
      std::vector<double> myWeights(myEdges.size(), 0.0);
      for (auto e = 0u; e < myEdges.size(); ++e) {
        auto n1 = allNeighbors.at(std::get<0>(myEdges[e])).size();
        auto n2 = allNeighbors.at(std::get<1>(myEdges[e])).size();
        if (n1 > s + 1) {
          myWeights[e] = boost::math::binomial_coefficient<double>(n1 - 1, s + 1);
        }
        if (n2 > s + 1) {
          myWeights[e] += boost::math::binomial_coefficient<double>(n2 - 1, s + 1);
        }
      }
      TIMER_START(this->m_tDist);
      this->fixWeightedImbalance(myEdges, myWeights, imbalanceThreshold);
      TIMER_PAUSE(this->m_tDist);
    }
  }
  if (directEdges) {
    TIMER_START(this->m_tDirect);
    this->storeRemovedEdges(std::move(myRemoved), std::move(myDSepSets), allNeighbors);
    TIMER_PAUSE(this->m_tDirect);
  }
  return this->constructSkeleton(std::move(allNeighbors));
}

template <typename Data, typename Var, typename Set>
PCStable2<Data, Var, Set>::PCStable2(
  const mxx::comm& comm,
  const Data& data,
  const double alpha,
  const Var maxConditioning
) : PCStableCommon<Data, Var, Set>(comm, data, alpha, maxConditioning)
{
}

template <typename Data, typename Var, typename Set>
BayesianNetwork<Var>
PCStable2<Data, Var, Set>::getSkeleton_parallel(
  const bool directEdges,
  const double imbalanceThreshold
) const
{
  // Initialize the learning in a similar fashion as done sequentially
  // Create nC2 edges - one for each unordered variable pair
  std::vector<std::tuple<Var, Var, double>> myEdges;
  std::unordered_map<Var, Set> allNeighbors;
  std::unordered_map<Var, Set> removedNeighbors;
  this->initializeLearning(myEdges, allNeighbors, removedNeighbors);
  // The following two data structures are used for getting
  // d-separating set and p-value for all the removed edges
  std::vector<std::tuple<Var, Var, double>> myRemoved;
  std::vector<Set> myDSepSets;
  std::set<std::tuple<Var, Var, double>> myBackwardRemoved;
  auto maxSize = std::min(this->m_maxConditioning, static_cast<Var>(this->m_allVars.size() - 2));
  for (auto s = 0u; (s <= maxSize) && mxx::any_of(myEdges.size() > 0, this->m_comm); ++s) {
    LOG_MESSAGE(debug, "Testing %u edges using sets of size %u", myEdges.size(), s);
    TIMER_DECLARE(tIter);
    for (auto& e : myEdges) {
      auto result = this->checkEdge(e, allNeighbors, removedNeighbors, s, false);
      std::get<2>(e) = result.first;
      // We need to store the removed edges since the d-separating
      // set may be required for directing edges later
      if (directEdges && (s > 0) &&
          this->m_data.isIndependent(this->m_alpha, result.first)) {
        TIMER_START(this->m_tDirect);
        Var x, y;
        double pv;
        std::tie(x, y, pv) = e;
        // We only use this information for directing colliders of type x-z-y
        // However, there is no need to store it if there are no candidates for z
        // i.e., no common neighbors between x and y
        if (!set_intersection(allNeighbors.at(x), allNeighbors.at(y)).empty()) {
          myRemoved.push_back((x < y) ? e : std::make_tuple(y, x, pv));
          myDSepSets.push_back(result.second);
        }
        TIMER_PAUSE(this->m_tDirect);
      }
    }
    this->m_comm.barrier();
    if (this->m_comm.is_first()) {
      TIMER_ELAPSED("Time taken in testing all sets of size " + std::to_string(s) + ": ", tIter);
    }
    TIMER_START(this->m_tSync);
    this->syncSets(removedNeighbors);
    TIMER_PAUSE(this->m_tSync);
    // Remove all the edges found to be independent on this processor
    // Also remove all the edges found to be independent on any other processor
    auto newEnd = std::remove_if(myEdges.begin(), myEdges.end(),
                                 [this, &removedNeighbors] (const std::tuple<Var, Var, double>& e)
                                                           { return this->m_data.isIndependent(this->m_alpha, std::get<2>(e)) ||
                                                                    removedNeighbors.at(std::get<0>(e)).contains(std::get<1>(e)); });
    myEdges.erase(newEnd, myEdges.end());
    for (auto& rn : removedNeighbors) {
      if (!rn.second.empty()) {
        allNeighbors[rn.first] = set_difference(allNeighbors.at(rn.first), rn.second);
        rn.second.clear();
      }
    }
    if (s == 0) {
      // Both direction of edges can be tested simultaneously when s > 0
      // Create a reverse edge for all the remaining edges
      auto origSize = myEdges.size();
      myEdges.resize(origSize * 2);
      Var x, y;
      for (auto f = 0u, b = origSize; f < origSize; ++f, ++b) {
        std::tie(x, y, std::ignore) = myEdges[f];
        myEdges[b] = std::make_tuple(y, x, 0.0);
      }
    }
    // Remove the edges that can no longer be tested using the first variable's neighborhood
    newEnd = std::remove_if(myEdges.begin(), myEdges.end(),
                            [&allNeighbors, &s] (const std::tuple<Var, Var, double>& e)
                                                { return allNeighbors.at(std::get<0>(e)).size() <= (s + 1); });
    myEdges.erase(newEnd, myEdges.end());
    std::set<std::tuple<Var, Var, double>> remove;
    // Now, copy all the backward edges, i.e., (x, y) such that x > y
    // and the neighborhood of x is the same as that of y
    std::copy_if(myEdges.begin(), myEdges.end(), std::inserter(remove, remove.begin()),
                 [&allNeighbors] (const std::tuple<Var, Var, double>& e)
                                 { return (std::get<0>(e) > std::get<1>(e)) &&
                                          (set_difference(allNeighbors.at(std::get<0>(e)),
                                                          allNeighbors.at(std::get<1>(e))) == Set({std::get<1>(e)})); });
    newEnd = std::remove_if(myEdges.begin(), myEdges.end(), [&remove] (const std::tuple<Var, Var, double>& e)
                                                                      { return remove.find(e) != remove.end(); });
    myEdges.erase(newEnd, myEdges.end());
    // Check if any of the backward edges removed in previous iterations
    // can be added back to the list of edges
    std::set<std::tuple<Var, Var, double>> addBackwardRemoved;
    Var x, y;
    for (const auto& e : myBackwardRemoved) {
      std::tie(x, y, std::ignore) = e;
      if ((allNeighbors.at(x).size() > (s + 1)) &&
          (set_difference(allNeighbors.at(x), allNeighbors.at(y)) != Set({y}))) {
        addBackwardRemoved.insert(e);
      }
    }
    // Remove the edges to be added back from the set of removed edges
    std::set<std::tuple<Var, Var, double>> remaining;
    std::set_difference(std::make_move_iterator(myBackwardRemoved.begin()),
                        std::make_move_iterator(myBackwardRemoved.end()),
                        addBackwardRemoved.begin(), addBackwardRemoved.end(),
                        std::inserter(remaining, remaining.begin()));
    myBackwardRemoved.clear();
    // Also add the edges removed in this iteration
    std::set_union(std::make_move_iterator(remaining.begin()),
                   std::make_move_iterator(remaining.end()),
                   remove.begin(), remove.end(),
                   std::inserter(myBackwardRemoved, myBackwardRemoved.begin()));
    // Re-add the backward edges which now need to be tested
    myEdges.insert(myEdges.end(), addBackwardRemoved.begin(), addBackwardRemoved.end());
    if (std::isgreaterequal(imbalanceThreshold, 0.0)) {
      std::vector<double> myWeights(myEdges.size());
      for (auto e = 0u; e < myEdges.size(); ++e) {
        auto n = allNeighbors.at(std::get<0>(myEdges[e])).size();
        myWeights[e] = boost::math::binomial_coefficient<double>(n - 1, s + 1);
      }
      TIMER_START(this->m_tDist);
      this->fixWeightedImbalance(myEdges, myWeights, imbalanceThreshold);
      TIMER_PAUSE(this->m_tDist);
    }
  }
  if (directEdges) {
    TIMER_START(this->m_tDirect);
    this->storeRemovedEdges(std::move(myRemoved), std::move(myDSepSets), allNeighbors, true);
    TIMER_PAUSE(this->m_tDirect);
  }
  return this->constructSkeleton(std::move(allNeighbors));
}

#endif // DETAIL_GLOBALLEARNING_HPP_
