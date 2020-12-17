/**
 * @file LemonTree.hpp
 * @brief Declaration of the class for learning module networks
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
#ifndef LEMONTREE_HPP_
#define LEMONTREE_HPP_

#include "ModuleNetworkLearning.hpp"

#include "utils/Timer.hpp"

#include <random>


template <typename Data, typename Var, typename Set>
class Module;

/**
 * @brief Class that implements learning of module networks
 *        using the approach of Lemon Tree.
 *
 * @tparam Data Type of the object which is used for querying the data.
 * @tparam Var Type of variable indices (expected to be an integer type).
 * @tparam Set Type of set container.
 */
template <typename Data, typename Var, typename Set>
class LemonTree : public ModuleNetworkLearning<Data, Var, Set> {
public:
  LemonTree(const mxx::comm&, const Data&);

  ~LemonTree();

protected:
  ModuleNetwork<Var>
  getNetwork_sequential(const pt::ptree&) const;

  ModuleNetwork<Var>
  getNetwork_parallel(const pt::ptree&) const;

private:
  std::list<std::list<Set>>
  clusterVarsGanesh(const pt::ptree&) const;

  void
  writeVarClusters(const std::string&, const std::list<std::list<Set>>&) const;

  std::vector<double>
  coclusteringMatrix(const std::list<std::list<Set>>&&, const double) const;

  std::multimap<Var, Var>
  clusterConsensus(const std::list<std::list<Set>>&&, const pt::ptree&) const;

  void
  writeConsensusCluster(const std::string&, const std::multimap<Var, Var>&) const;

  std::list<std::list<Set>>
  clusterObsGanesh(const uint32_t, const uint32_t, const uint32_t, const uint32_t, std::mt19937* const, const Set&) const;

  std::list<Module<Data, Var, Set>>
  learnModules(const std::multimap<Var, Var>&&, const pt::ptree&) const;

  void
  writeParents(std::ofstream&, const std::unordered_map<Var, double>&, const uint32_t, const double = 1.0) const;

  void
  writeModules(const std::string&, const std::list<Module<Data, Var, Set>>&) const;

private:
  TIMER_DECLARE(m_tWrite, mutable);
  TIMER_DECLARE(m_tGanesh, mutable);
  TIMER_DECLARE(m_tConsensus, mutable);
  TIMER_DECLARE(m_tModules, mutable);
}; // class LemonTree

#include "detail/LemonTree.hpp"

#endif // LEMONTREE_HPP_
