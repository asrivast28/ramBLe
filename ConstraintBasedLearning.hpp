/**
 * @file ConstraintBasedLearning.hpp
 * @brief Declaration of the classes for constraint-based learning.
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
#ifndef CONSTRAINTBASEDLEARNING_HPP_
#define CONSTRAINTBASEDLEARNING_HPP_

#include "BayesianNetwork.hpp"

#include "mxx/comm.hpp"

#include <map>
#include <unordered_map>


/**
 * @brief Abstract base class for causal discovery using constraint-based learning.
 *
 * @tparam Data Type of the object which is used for querying the data.
 * @tparam Var Type of variable indices (expected to be an integer type).
 * @tparam Set Type of set container.
 */
template <typename Data, typename Var, typename Set>
class ConstraintBasedLearning {
public:
  ConstraintBasedLearning(const mxx::comm&, const Data&, const Var);

  const Set&
  getPC(const Var) const;

  const Set&
  getMB(const Var) const;

  std::vector<std::tuple<double, Var, Var, Var>>
  findVStructures(const Var) const;

  BayesianNetwork<Var>
  getNetwork(const bool, const bool, const double = 0.0) const;

  virtual
  ~ConstraintBasedLearning() { }

protected:
  Set
  getCandidates(const Var) const;

  virtual
  Set
  getCandidatePC(const Var, Set&&) const = 0;

  virtual
  Set
  getCandidateMB(const Var, Set&&) const = 0;

  virtual
  BayesianNetwork<Var>
  getSkeleton_sequential() const;

  void
  parallelInitialize(std::vector<std::tuple<Var, Var, double>>&, std::unordered_map<Var, Set>&) const;

  std::vector<std::pair<Var, Var>>
  symmetryCorrect(const std::unordered_map<Var, Set>&&, const std::set<std::pair<Var, Var>>&&) const;

  void
  syncSets(std::unordered_map<Var, Set>&) const;

  void
  syncMissingSets(const std::vector<std::tuple<Var, Var, double>>&, std::unordered_map<Var, Set>&) const;

  bool
  fixImbalance(std::vector<std::tuple<Var, Var, double>>&, const double) const;

  virtual
  BayesianNetwork<Var>
  getSkeleton_parallel(const double, std::unordered_map<Var, Set>&, std::unordered_map<Var, Set>&) const;

private:
  Set&
  getCandidatePC_cache(const Var, Set&&) const;

  void
  symmetryCorrectPC(const Var, Set&) const;

  Set&
  getCandidateMB_cache(const Var, Set&&) const;

  void
  symmetryCorrectMB(const Var, Set&) const;

  double
  colliderPValue(const Var, const Var, const Var) const;

  std::vector<std::tuple<double, Var, Var, Var>>
  findVStructures() const;

protected:
  const mxx::comm& m_comm;
  const Data m_data;
  Set m_allVars;
  const Var m_maxConditioning;

private:
  mutable std::unordered_map<Var, Set> m_cachedPC;
  mutable std::unordered_map<Var, Set> m_cachedMB;
  mutable std::unordered_map<Var, bool> m_cachedPCSymmetric;
  mutable std::unordered_map<Var, bool> m_cachedMBSymmetric;
}; // class ConstraintBasedLearning

#include "detail/ConstraintBasedLearning.hpp"

#endif // CONSTRAINTBASEDLEARNING_HPP_
