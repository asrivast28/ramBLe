/**
 * @file MBDiscovery.hpp
 * @brief Implementation of the MBDiscovery functions.
 */
#ifndef DETAIL_MBDISCOVERY_HPP_
#define DETAIL_MBDISCOVERY_HPP_

#include "../SetUtils.hpp"

#include "utils/Logging.hpp"


template <typename DataType, typename VarType, typename SetType>
/**
 * @brief Constructs the object with the given data.
 *
 * @param data Reference to an object of the DataType.
 */
MBDiscovery<DataType, VarType, SetType>::MBDiscovery(
  const DataType& data
) : m_data(data),
    m_vars()
{
  for (auto i = 0u; i < data.numVars(); ++i) {
    m_vars.insert(m_vars.end(), i);
  }
}

template <typename DataType, typename VarType, typename SetType>
/**
 * @brief Function for getting all the candidates for the given target variable.
 *
 * @param target The index of the target variable.
 *
 * @return The indices of all the variables except the target.
 */
SetType
MBDiscovery<DataType, VarType, SetType>::getCandidates(
  const VarType target
) const
{
  auto candidates = m_vars;
  candidates.erase(target);
  return candidates;
}

template <typename DataType, typename VarType, typename SetType>
/**
 * @brief Finds the candidate PC for a variable, and caches the result.
 *
 * @param target The index of the target variable.
 * @param candidates The indices of all the candidate variables.
 *
 * @return A set containing the indices of all the variables
 *         in the candidate PC of the given target variable.
 */
SetType
MBDiscovery<DataType, VarType, SetType>::getCandidatePC_cache(
  const VarType target,
  SetType candidates
) const
{
  auto cacheIt = m_cachedPC.find(target);
  if (cacheIt == m_cachedPC.end()) {
    auto cpc = this->getCandidatePC(target, std::move(candidates));
    m_cachedPC.insert(cacheIt, std::make_pair(target, cpc));
    return cpc;
  }
  else {
    LOG_MESSAGE(debug, "* Found candidate PC for %s in the cache", this->m_data.varName(target))
    return cacheIt->second;
  }
}

template <typename DataType, typename VarType, typename SetType>
/**
 * @brief Performs symmetry correction on the given PC set for the given variable.
 *
 * @param target The index of the target variable.
 * @param cpc The set containing the indices of all the variables in
 *            the candidate PC set.
 */
void
MBDiscovery<DataType, VarType, SetType>::symmetryCorrectPC(
  const VarType target,
  SetType& cpc
) const
{
  auto initial = cpc;
  for (const VarType x: initial) {
    auto candidatesX = this->getCandidates(x);
    auto cpcX = this->getCandidatePC_cache(x, std::move(candidatesX));
    if (!set_contains(cpcX, target)) {
      LOG_MESSAGE(info, "- Removing %s from the PC of %s (asymmetry)", this->m_data.varName(x), this->m_data.varName(target));
      cpc.erase(x);
    }
  }
}

template <typename DataType, typename VarType, typename SetType>
/**
 * @brief The top level function for getting the correct PC set
 *        for the given target variable.
 *
 * @param target The index of the target variable.
 *
 * @return A set containing the indices of all the variables
 *         in the correct PC set of the target variable.
 */
SetType
MBDiscovery<DataType, VarType, SetType>::getPC(
  const VarType target
) const
{
  auto candidates = this->getCandidates(target);
  auto cpc = this->getCandidatePC_cache(target, std::move(candidates));
  this->symmetryCorrectPC(target, cpc);
  return cpc;
}

template <typename DataType, typename VarType, typename SetType>
/**
 * @brief Finds the candidate MB for a variable, and caches the result.
 *
 * @param target The index of the target variable.
 * @param candidates The indices of all the candidate variables.
 *
 * @return A set containing the indices of all the variables
 *         in the candidate MB of the given target variable.
 */
SetType
MBDiscovery<DataType, VarType, SetType>::getCandidateMB_cache(
  const VarType target,
  SetType candidates
) const
{
  auto cacheIt = m_cachedMB.find(target);
  if (cacheIt == m_cachedMB.end()) {
    auto cpc = this->getCandidateMB(target, std::move(candidates));
    m_cachedMB.insert(cacheIt, std::make_pair(target, cpc));
    return cpc;
  }
  else {
    LOG_MESSAGE(debug, "* Found candidate MB for %s in the cache", this->m_data.varName(target));
    return cacheIt->second;
  }
}

template <typename DataType, typename VarType, typename SetType>
/**
 * @brief Symmetry corrects the candidate MB of the target variable.
 *
 * @param target The index of the target variable.
 * @param cmb The indices of the variables in the candidate MB of the target variable.
 *                The function removes the indices from the candidate set.
 */
void
MBDiscovery<DataType, VarType, SetType>::symmetryCorrectMB(
  const VarType target,
  SetType& cmb
) const
{
  auto initial = cmb;
  for (const VarType x: initial) {
    auto candidatesX = this->getCandidates(x);
    auto cmbX = this->getCandidateMB_cache(x, std::move(candidatesX));
    if (!set_contains(cmbX, target)) {
      LOG_MESSAGE(info, "- Removing %s from the MB of %s (asymmetry)", this->m_data.varName(x), this->m_data.varName(target));
      cmb.erase(x);
    }
  }
}

template <typename DataType, typename VarType, typename SetType>
/**
 * @brief Top level function for getting the MB of the target variable.
 *
 * @param target The index of the target variable.
 */
SetType
MBDiscovery<DataType, VarType, SetType>::getMB(
  const VarType target
) const
{
  auto candidates = this->getCandidates(target);
  auto cmb = this->getCandidateMB_cache(target, std::move(candidates));
  this->symmetryCorrectMB(target, cmb);
  return cmb;
}

#endif // DETAIL_MBDISCOVERY_HPP_
