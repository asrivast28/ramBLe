/**
 * @file DataQuery.hpp
 * @brief Declaration of the functions used for querying the data.
 */
#ifndef DATAQUERY_HPP_
#define DATAQUERY_HPP_

#include <cstdint>
#include <limits>
#include <set>
#include <string>
#include <vector>


/**
 * @brief Class that provides functionality for querying the data.
 *
 * @tparam Counter Type of the object that provides counting queries.
 * @tparam Var Type of the variables (expected to be an integral type).
 */
template <typename Counter, typename Var>
class DataQuery {
public:
  DataQuery();

  DataQuery(const Counter&, const std::vector<std::string>&, const double = 0.05);

  Var
  numVars() const;

  uint32_t
  numRows() const;

  const std::string&
  varName(const Var) const;

  template <typename Set = std::set<Var>>
  std::vector<std::string>
  varNames(const Set&) const;

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

  template <typename Set>
  double
  minAssocScore(const Var, const Var, const Set&, const uint32_t = std::numeric_limits<uint32_t>::max()) const;

  template <typename Set>
  double
  minAssocScore(const Var, const Var, const Set&, const Set&, const uint32_t = std::numeric_limits<uint32_t>::max()) const;

  template <typename Set>
  std::pair<double, Set>
  minAssocScoreSubset(const Var, const Var, const Set&, const uint32_t = std::numeric_limits<uint32_t>::max()) const;

  template <typename Set>
  bool
  isIndependentAnySubset(const Var, const Var, const Set&, const uint32_t = std::numeric_limits<uint32_t>::max()) const;

  template <typename Set>
  bool
  isIndependentAnySubset(const Var, const Var, const Set&, const Set&, const uint32_t = std::numeric_limits<uint32_t>::max()) const;

private:
  Counter m_counter;
  std::vector<std::string> m_varNames;
  double m_threshold;
};

#include "detail/GSquare.hpp"
#include "detail/DataQuery.hpp"

#endif // DATAQUERY_HPP_
