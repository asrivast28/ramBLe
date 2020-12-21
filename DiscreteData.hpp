/**
 * @file DiscreteData.hpp
 * @brief Declaration of the functions used for querying the data.
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
#ifndef DISCRETEDATA_HPP_
#define DISCRETEDATA_HPP_

#include "mxx/comm.hpp"
#include "utils/Timer.hpp"

#include <cstdint>
#include <limits>
#include <set>
#include <string>
#include <vector>


/**
 * @brief Class that provides functionality for querying discrete data.
 *
 * @tparam Counter Type of the object that provides counting queries.
 * @tparam Var Type of the variables (expected to be an integral type).
 */
template <typename Counter, typename Var>
class DiscreteData {
public:
  DiscreteData();

  DiscreteData(const Counter&, const std::vector<std::string>&);

  Var
  numVars() const;

  uint32_t
  numObs() const;

  const std::string&
  varName(const Var) const;

  template <typename Set = std::set<Var>>
  std::vector<std::string>
  varNames(const Set&) const;

  const std::vector<std::string>&
  varNames() const;

  Var
  varIndex(const std::string&) const;

  template <typename Set = std::set<Var>>
  Set
  varIndices(const std::vector<std::string>&) const;

  template <typename Set = std::set<Var>>
  double
  pValue(const Var, const Var, const Set& = Set()) const;

  template <typename Set = std::set<Var>>
  bool
  isIndependent(const double, const Var, const Var, const Set& = Set()) const;

  bool
  isIndependent(const double, const double) const;

  template <template <typename...> class SetType, typename... Args>
  double
  maxPValue(const double, const Var, const Var, const SetType<Var, Args...>&, const Var, const Var = 0u) const;

  template <template <typename...> class SetType, typename... Args>
  double
  maxPValue(const double, const Var, const Var, const SetType<Var, Args...>&, const SetType<Var, Args...>&, const Var) const;

  template <template <typename...> class SetType, typename... Args>
  std::pair<double, SetType<Var, Args...>>
  maxPValueSubset(const double, const Var, const Var, const SetType<Var, Args...>&, const Var, const Var = 0u) const;

  template <typename Set>
  bool
  isIndependentAnySubset(const double, const Var, const Var, const Set&, const Var, const Var = 0u) const;

  template <template <typename...> class SetType, typename... Args>
  bool
  isIndependentAnySubset(const double, const Var, const Var, const SetType<Var, Args...>&, const Var, const mxx::comm&) const;

  template <typename Set>
  bool
  isIndependentAnySubset(const double, const Var, const Var, const Set&, const Set&, const Var) const;

  ~DiscreteData();

private:
  uint32_t
  testThreshold(const uint32_t = 1u) const;

private:
  const Counter m_counter;
  const std::vector<std::string> m_varNames;
  TIMER_DECLARE(m_timer, mutable);
};

#include "detail/GSquare.hpp"
#include "detail/DiscreteData.hpp"

#endif // DISCRETEDATA_HPP_
