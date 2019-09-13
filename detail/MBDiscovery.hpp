/**
 * @file MBDiscovery.hpp
 * @brief Implementation of the MBDiscovery functions.
 */
#ifndef DETAIL_MBDISCOVERY_HPP_
#define DETAIL_MBDISCOVERY_HPP_

#include "utils/Logging.hpp"


template <typename DataType, typename VarType>
/**
 * @brief Constructs the object with the given data.
 *
 * @param data Reference to an object of the DataType.
 */
MBDiscovery<DataType, VarType>::MBDiscovery(
  const DataType& data
) : m_data(data),
    m_vars()
{
  for (auto i = 0u; i < data.numVars(); ++i) {
    m_vars.insert(m_vars.end(), i);
  }
}

template <typename DataType, typename VarType>
/**
 * @brief Function for getting all the candidates for the given target variable.
 *
 * @param target The index of the target variable.
 *
 * @return The indices of all the variables except the target.
 */
std::set<VarType>
MBDiscovery<DataType, VarType>::getCandidates(
  const VarType target
) const
{
  auto candidates = m_vars;
  candidates.erase(target);
  return candidates;
}

template <typename DataType, typename VarType>
/**
 * @brief Finds the candidate MB for a variable, and caches the result.
 *
 * @param target The index of the target variable.
 * @param candidates The indices of all the candidate variables.
 *
 * @return A set containing the indices of all the variables
 *         in the candidate MB of the given target variable.
 */
std::set<VarType>
MBDiscovery<DataType, VarType>::getCandidateMB_cache(
  const VarType target,
  std::set<VarType> candidates
) const
{
  auto cacheIt = m_cachedMB.find(target);
  if (cacheIt == m_cachedMB.end()) {
    auto cpc = this->getCandidateMB(target, std::move(candidates));
    m_cachedMB.insert(cacheIt, std::make_pair(target, cpc));
    return cpc;
  }
  else {
    LOG_MESSAGE(debug, "Found candidate MB for %s in the cache", this->m_data.varName(target));
    return cacheIt->second;
  }
}

template <typename DataType, typename VarType>
/**
 * @brief Symmetry corrects the candidate MB of the target variable.
 *
 * @param target The index of the target variable.
 * @param cmb The indices of the variables in the candidate MB of the target variable.
 *                The function removes the indices from the candidate set.
 */
void
MBDiscovery<DataType, VarType>::symmetryCorrectMB(
  const VarType target,
  std::set<VarType>& cmb
) const
{
  auto initial = cmb;
  for (const VarType x: initial) {
    auto candidatesX = this->getCandidates(x);
    auto cmbX = this->getCandidateMB_cache(x, std::move(candidatesX));
    if (cmbX.find(target) == cmbX.end()) {
      LOG_MESSAGE(info, "Symmetry Correction: Removing %s from the candidate MB", this->m_data.varName(x));
      cmb.erase(x);
    }
  }
}

template <typename DataType, typename VarType>
/**
 * @brief Top level function for getting the MB of the target variable.
 *
 * @param target The index of the target variable.
 */
std::set<VarType>
MBDiscovery<DataType, VarType>::getMB(
  const VarType target
) const
{
  auto candidates = this->getCandidates(target);
  auto cmb = this->getCandidateMB_cache(target, std::move(candidates));
  this->symmetryCorrectMB(target, cmb);
  return cmb;
}

#endif // DETAIL_MBDISCOVERY_HPP_
