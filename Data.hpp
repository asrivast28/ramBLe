/**
 * @file Data.hpp
 * @brief Declaration of the functions used for querying the data.
 */
#ifndef DATA_HPP_
#define DATA_HPP_

#include <cstdint>
#include <limits>
#include <set>
#include <string>
#include <vector>


/**
 * @brief Class that provides functionality for querying the data.
 *
 * @tparam CounterType Type of the object that provides counting queries.
 * @tparam VarType Type of the variables (expected to be an integral type).
 */
template <typename CounterType, typename VarType>
class Data {
public:
  Data();

  Data(const CounterType&, const std::vector<std::string>&, const double = 0.05);

  VarType
  numVars() const;

  uint32_t
  numRows() const;

  const std::string&
  varName(const VarType) const;

  template <typename SetType = std::set<VarType>>
  std::vector<std::string>
  varNames(const SetType&) const;

  VarType
  varIndex(const std::string&) const;

  template <typename SetType = std::set<VarType>>
  SetType
  varIndices(const std::vector<std::string>&) const;

  template <typename SetType = std::set<VarType>>
  double
  pValue(const VarType, const VarType, const SetType& = SetType()) const;

  template <typename SetType = std::set<VarType>>
  double
  assocScore(const VarType, const VarType, const SetType& = SetType()) const;

  template <typename SetType = std::set<VarType>>
  bool
  isIndependent(const VarType, const VarType, const SetType& = SetType()) const;

  bool
  isIndependent(const double) const;

  template <typename SetType>
  double
  minAssocScore(const VarType, const VarType, const SetType&, const uint32_t = std::numeric_limits<uint32_t>::max()) const;

  template <typename SetType>
  double
  minAssocScore(const VarType, const VarType, const SetType&, const SetType&, const uint32_t = std::numeric_limits<uint32_t>::max()) const;

  template <typename SetType>
  std::pair<double, SetType>
  minAssocScoreSubset(const VarType, const VarType, const SetType&, const uint32_t = std::numeric_limits<uint32_t>::max()) const;

  template <typename SetType>
  bool
  isIndependentAnySubset(const VarType, const VarType, const SetType&, const uint32_t = std::numeric_limits<uint32_t>::max()) const;

  template <typename SetType>
  bool
  isIndependentAnySubset(const VarType, const VarType, const SetType&, const SetType&, const uint32_t = std::numeric_limits<uint32_t>::max()) const;

private:
  template <typename SetType>
  std::pair<uint32_t, double>
  gSquare(const VarType, const VarType, const SetType&) const;

private:
  CounterType m_counter;
  std::vector<std::string> m_varNames;
  double m_threshold;
};

#include "detail/Data.hpp"

#endif // DATA_HPP_
