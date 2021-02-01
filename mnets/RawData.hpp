/**
 * @file RawData.hpp
 * @brief Declaration of the functions used for querying raw data.
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
#ifndef RAWDATA_HPP_
#define RAWDATA_HPP_

#include "utils/Logging.hpp"

#include <cstdint>
#include <set>
#include <string>
#include <vector>


/**
 * @brief Class that provides functionality for querying raw data.
 *
 * @tparam DataType Type of the underlying raw data.
 * @tparam Var Type of the variables (expected to be an integral type).
 */
template <typename DataType, typename Var>
class RawData {
public:
  RawData(const std::vector<DataType>&, const std::vector<std::string>&, const Var, const Var);

  const std::vector<DataType>&
  raw() const;

  const std::string&
  varName(const Var) const;

  const std::vector<std::string>&
  varNames() const;

  template <typename Set = std::set<Var>>
  std::vector<std::string>
  varNames(const Set&) const;

  Var
  varIndex(const std::string&) const;

  Var
  numVars() const;

  Var
  numObs() const;

  double
  operator()(const Var, const Var) const;

  ~RawData();

private:
  const std::vector<DataType>& m_raw;
  const std::vector<std::string> m_varNames;
  const Var m_nvars;
  const Var m_nobs;
};

template <typename DataType, typename Var>
/**
 * @brief Constructs the data provider object.
 *
 * @param raw A reference to the raw data set.
 * @param varNames Names of the variables in the data set.
 * @param n The number of variables in the data set.
 * @param m The number of observations in the data set.
 */
RawData<DataType, Var>::RawData(
  const std::vector<DataType>& raw,
  const std::vector<std::string>& varNames,
  const Var n,
  const Var m
) : m_raw(raw),
    m_varNames(varNames),
    m_nvars(n),
    m_nobs(m)
{
}

template <typename DataType, typename Var>
/**
 * @brief Default destructor.
 */
RawData<DataType, Var>::~RawData(
)
{
}

template <typename DataType, typename Var>
const std::vector<DataType>&
RawData<DataType, Var>::raw(
) const
{
  return m_raw;
}

template <typename DataType, typename Var>
/**
 * @brief Returns the name of a variable.
 *
 * @param x The index of the query variable.
 *
 * @return The name of the query variable.
 */
const std::string&
RawData<DataType, Var>::varName(
  const Var x
) const
{
  LOG_MESSAGE_IF(x >= m_varNames.size(), error, "Variable index %d out of range.", static_cast<uint32_t>(x));
  return m_varNames[x];
}

template <typename DataType, typename Var>
/**
 * @brief Returns the names of all the variables in the data set.
 */
const std::vector<std::string>&
RawData<DataType, Var>::varNames(
) const
{
  return m_varNames;
}

template <typename Counter, typename Var>
/**
 * @brief Returns the names of all the variables in the given set.
 *
 * @tparam Set The type of container for the variable indices.
 * @param vars The indices of all the query variable.
 *
 * @return The name of all the query variables.
 */
template <typename Set>
std::vector<std::string>
RawData<Counter, Var>::varNames(
  const Set& vars
) const
{
  std::vector<std::string> names(vars.size());
  auto i = 0u;
  for (const auto var : vars) {
    LOG_MESSAGE_IF(var >= m_varNames.size(), error, "Variable index %d out of range.", static_cast<uint32_t>(var));
    names[i++] = m_varNames[var];
  }
  return names;
}

template <typename Counter, typename Var>
/**
 * @brief Returns the index of a variable.
 *
 * @param name The name of the query variable.
 *
 * @return The index of the query variable.
 */
Var
RawData<Counter, Var>::varIndex(
  const std::string& name
) const
{
  Var x = 0u;
  for (const auto& var : m_varNames) {
    if (var.compare(name) == 0) {
      break;
    }
    ++x;
  }
  LOG_MESSAGE_IF(x == numVars(), error, "Variable with name %s not found.", name);
  return x;
}

template <typename DataType, typename Var>
/**
 * @brief Returns the number of variables in the data set.
 */
Var
RawData<DataType, Var>::numVars(
) const
{
  return m_nvars;
}

template <typename DataType, typename Var>
/**
 * @brief Returns the number of observations in the data set.
 */
Var
RawData<DataType, Var>::numObs(
) const
{
  return m_nobs;
}

template <typename DataType, typename Var>
/**
 * @brief Returns the data point at the given index.
 *
 * @param i The variable index.
 * @param j The observation index.
 */
double
RawData<DataType, Var>::operator()(
  const Var i,
  const Var j
) const
{
  return m_raw[i * m_nobs + j];
}

#endif // RAWDATA_HPP_
