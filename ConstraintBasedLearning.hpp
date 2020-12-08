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
#include <stdexcept>
#include <unordered_map>


/**
 * @brief Exception to be used for not implemented functions.
 */
class NotImplementedError : public std::logic_error {
public:
  NotImplementedError(const std::string&);
}; // class NotImplementedError

/**
 * @brief Abstract base class for causal discovery using constraint-based learning.
 *        All the constraint-based learning algorithm implementations should be a
 *        descendant of this class.
 *
 * @tparam Data Type of the object which is used for querying the data.
 * @tparam Var Type of variable indices (expected to be an integer type).
 * @tparam Set Type of set container.
 */
template <typename Data, typename Var, typename Set>
class ConstraintBasedLearning {
public:
  ConstraintBasedLearning(const mxx::comm&, const Data&, const Var);

  virtual
  const Set&
  getPC(const Var) const = 0;

  virtual
  const Set&
  getMB(const Var) const = 0;

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
  BayesianNetwork<Var>
  getSkeleton_sequential(const bool) const = 0;

  bool
  fixImbalance(std::vector<std::tuple<Var, Var, double>>&, const double) const;

  virtual
  BayesianNetwork<Var>
  getSkeleton_parallel(const bool, const double) const = 0;

  virtual
  std::pair<bool, double>
  checkCollider(const Var, const Var, const Var) const = 0;

private:
  std::vector<std::tuple<double, Var, Var, Var>>
  findVStructures() const;

protected:
  const mxx::comm& m_comm;
  const Data m_data;
  Set m_allVars;
  const Var m_maxConditioning;
}; // class ConstraintBasedLearning

#include "detail/ConstraintBasedLearning.hpp"

#endif // CONSTRAINTBASEDLEARNING_HPP_
