/**
 * @file LocalLearning.hpp
 * @brief Declaration of the classes for local-to-global causal discovery.
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
#ifndef LOCALLEARNING_HPP_
#define LOCALLEARNING_HPP_

#include "ConstraintBasedLearning.hpp"


/**
 * @brief Abstract base class for local-to-global causal discovery.
 *
 * @tparam Data Type of the object which is used for querying the data.
 * @tparam Var Type of variable indices (expected to be an integer type).
 * @tparam Set Type of set container.
 */
template <typename Data, typename Var, typename Set>
class LocalLearning : public ConstraintBasedLearning<Data, Var, Set> {
public:
  LocalLearning(const mxx::comm&, const Data&, const Var);

  const Set&
  getPC(const Var) const override;

  const Set&
  getMB(const Var) const override;

  virtual
  ~LocalLearning() { }

protected:
  virtual
  Set
  getCandidatePC(const Var, Set&&) const = 0;

  virtual
  Set
  getCandidateMB(const Var, Set&&) const = 0;

  void
  parallelInitialize(std::vector<std::tuple<Var, Var, double>>&, std::unordered_map<Var, Set>&) const;

  void
  syncSets(std::unordered_map<Var, Set>&) const;

  void
  syncMissingSets(const std::vector<std::tuple<Var, Var, double>>&, std::unordered_map<Var, Set>&) const;

  std::vector<std::pair<Var, Var>>
  symmetryCorrect(const std::unordered_map<Var, Set>&&, const std::set<std::pair<Var, Var>>&&) const;

private:
  Set&
  getCandidatePC_cache(const Var, Set&&) const;

  void
  symmetryCorrectPC(const Var, Set&) const;

  Set&
  getCandidateMB_cache(const Var, Set&&) const;

  void
  symmetryCorrectMB(const Var, Set&) const;

  BayesianNetwork<Var>
  getSkeleton_sequential(const bool) const override;

protected:
  mutable std::unordered_map<Var, Set> m_cachedPC;
  mutable std::unordered_map<Var, Set> m_cachedMB;
  mutable std::unordered_map<Var, bool> m_cachedPCSymmetric;
  mutable std::unordered_map<Var, bool> m_cachedMBSymmetric;
}; // class LocalLearning

#include "detail/LocalLearning.hpp"

#endif // LOCALLEARNING_HPP_
