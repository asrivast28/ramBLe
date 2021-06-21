/**
 * @file GlobalLearning.hpp
 * @brief Declaration of the classes for global learning algorithms.
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
#ifndef GLOBALLEARNING_HPP_
#define GLOBALLEARNING_HPP_

#include "ConstraintBasedLearning.hpp"


/**
 * @brief Abstract base class for global causal discovery using constraint-based algorithms.
 *
 * @tparam Data Type of the object which is used for querying the data.
 * @tparam Var Type of variable indices (expected to be an integer type).
 * @tparam Set Type of set container.
 */
template <typename Data, typename Var, typename Set>
class GlobalLearning : public ConstraintBasedLearning<Data, Var, Set> {
public:
  GlobalLearning(const mxx::comm&, const Data&, const double, const Var);

  const Set&
  getPC(const Var) const override;

  const Set&
  getMB(const Var) const override;

  virtual
  ~GlobalLearning();

protected:
  void
  initializeLearning(std::vector<std::tuple<Var, Var, double>>&, std::unordered_map<Var, Set>&, std::unordered_map<Var, Set>&) const;

  bool
  fixWeightedImbalance(std::vector<std::tuple<Var, Var, double>>&, const std::vector<double>&, const double) const;

  void
  storeRemovedEdges(std::vector<std::tuple<Var, Var, double>>&&, std::vector<Set>&&, const std::unordered_map<Var, Set>&) const;

  BayesianNetwork<Var>
  constructSkeleton(const std::unordered_map<Var, Set>&&) const;

private:
  std::pair<bool, double>
  checkCollider(const Var, const Var, const Var) const override;

protected:
  mutable std::unordered_map<Var, Set> m_cachedNeighbors;
  mutable std::vector<std::tuple<Var, Var, double, Set>> m_removedEdges;
  TIMER_DECLARE(m_tSync, mutable);
  TIMER_DECLARE(m_tDist, mutable);
  TIMER_DECLARE(m_tRemoved, mutable);
}; // class GlobalLearning


/**
 * @brief Class that implements PC-stable algorithm for learning BN skeleton,
 *        as described by Colombo et al.
 *
 * @tparam Data Type of the object which is used for querying the data.
 * @tparam Var Type of variable indices (expected to be an integer type).
 * @tparam Set Type of set container.
 */
template <typename Data, typename Var, typename Set>
class PCStableCommon : public GlobalLearning<Data, Var, Set> {
public:
  PCStableCommon(const mxx::comm&, const Data&, const double, const Var);

protected:
  std::pair<double, Set>
  checkEdge(const std::tuple<Var, Var, double>&, const std::unordered_map<Var, Set>&, std::unordered_map<Var, Set>&, const uint32_t, const bool) const;

  BayesianNetwork<Var>
  getSkeleton_sequential(const bool) const override;
}; // class PCStableCommon

template <typename Data, typename Var, typename Set>
class PCStable : public PCStableCommon<Data, Var, Set> {
public:
  PCStable(const mxx::comm&, const Data&, const double = 0.05, const Var = std::numeric_limits<Var>::max());

private:
  void
  syncSets(std::unordered_map<Var, Set>&) const;

  BayesianNetwork<Var>
  getSkeleton_parallel(const bool, const double) const override;
}; // class PCStable

template <typename Data, typename Var, typename Set>
class PCStable2 : public PCStableCommon<Data, Var, Set> {
public:
  PCStable2(const mxx::comm&, const Data&, const double = 0.05, const Var = std::numeric_limits<Var>::max());

private:
  BayesianNetwork<Var>
  getSkeleton_parallel(const bool, const double) const override;
}; // class PCStable2

#include "detail/GlobalLearning.hpp"

#endif // GLOBALLEARNING_HPP_
