/**
 * @file DiscreteData.hpp
 * @brief Declaration of the functions used for querying the data.
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

  DiscreteData(const Counter&, const std::vector<std::string>&, const double = 0.05);

  Var
  numVars() const;

  uint32_t
  numRows() const;

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
  double
  assocScore(const Var, const Var, const Set& = Set()) const;

  template <typename Set = std::set<Var>>
  bool
  isIndependent(const Var, const Var, const Set& = Set()) const;

  bool
  isIndependent(const double) const;

  template <template <typename...> class SetType, typename... Args>
  double
  minAssocScore(const Var, const Var, const SetType<Var, Args...>&, const Var) const;

  template <template <typename...> class SetType, typename... Args>
  double
  minAssocScore(const Var, const Var, const SetType<Var, Args...>&, const SetType<Var, Args...>&, const Var) const;

  template <template <typename...> class SetType, typename... Args>
  std::pair<double, SetType<Var, Args...>>
  minAssocScoreSubset(const Var, const Var, const SetType<Var, Args...>&, const Var) const;

  template <typename Set>
  bool
  isIndependentAnySubset(const Var, const Var, const Set&, const Var) const;

  template <template <typename...> class SetType, typename... Args>
  bool
  isIndependentAnySubset(const Var, const Var, const SetType<Var, Args...>&, const Var, const mxx::comm&) const;

  template <typename Set>
  bool
  isIndependentAnySubset(const Var, const Var, const Set&, const Set&, const Var) const;

  ~DiscreteData();

private:
  uint32_t
  testThreshold(const uint32_t = 1u) const;

private:
  Counter m_counter;
  std::vector<std::string> m_varNames;
  double m_threshold;
  TIMER_DECLARE(m_timer, mutable);
};

#include "detail/GSquare.hpp"
#include "detail/DiscreteData.hpp"

#endif // DISCRETEDATA_HPP_
