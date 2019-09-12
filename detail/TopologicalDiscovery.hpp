/**
 * @file TopologicalDiscovery.hpp
 * @brief Implementation of the TopologicalDiscovery class and all the
 *        classes that are derived from it.
 */
#ifndef DETAIL_TOPOLOGICALDISCOVERY_HPP_
#define DETAIL_TOPOLOGICALDISCOVERY_HPP_

#include "utils/Logging.hpp"

#include <algorithm>


template <typename DataType, typename VarType>
/**
 * @brief Constructs the object with the given data.
 *
 * @param data Reference to an object of the DataType.
 */
TopologicalDiscovery<DataType, VarType>::TopologicalDiscovery(
  const DataType& data
) : MBDiscovery<DataType, VarType>(data)
{
}

template <typename DataType, typename VarType>
/**
 * @brief Finds the candidate PC for every variable, and caches the result.
 *
 * @param target The index of the target variable.
 * @param candidates The indices of all the candidate variables.
 *
 * @return A set containing the indices of all the variables
 *         in the candidate PC of the given target variable.
 */
std::set<VarType>
TopologicalDiscovery<DataType, VarType>::getCandidatePC(
  const VarType target,
  std::set<VarType> candidates
) const
{
  auto cacheIt = m_cachedPC.find(target);
  if (cacheIt == m_cachedPC.end()) {
    auto cpc = this->getCandidatePC_impl(target, std::move(candidates));
    m_cachedPC.insert(cacheIt, std::make_pair(target, cpc));
    return cpc;
  }
  else {
    LOG_MESSAGE(debug, "Found candidate PC for %s in the cache", this->m_data.varName(target));
    return cacheIt->second;
  }
}

template <typename DataType, typename VarType>
/**
 * @brief Removes false positives from the given candidate PC set
 *        for the given target variable.
 *
 * @param target The index of the target variable.
 * @param cpc The set containing the indices of all the variables in
 *            the candidate PC set.
 *
 * @return A set containing the indices of all the variables removed
 *         from the candidate PC set.
 */
std::set<VarType>
TopologicalDiscovery<DataType, VarType>::removeFalsePC(
  const VarType target,
  std::set<VarType>& cpc
) const
{
  std::set<VarType> removed;
  auto initial = cpc;
  for (const VarType x: initial) {
    cpc.erase(x);
    LOG_MESSAGE(debug, "False Positive: Testing %s for removal", this->m_data.varName(x));
    if (this->m_data.isIndependentAnySubset(target, x, cpc)) {
      LOG_MESSAGE(info, "False Positive: Removing %s from the candidate PC of %s", this->m_data.varName(x), this->m_data.varName(target));
      removed.insert(x);
    }
    else {
      cpc.insert(x);
    }
  }
  return removed;
}

template <typename DataType, typename VarType>
/**
 * @brief Performs symmetry correction on the given PC set for the given variable.
 *
 * @param target The index of the target variable.
 * @param cpc The set containing the indices of all the variables in
 *            the candidate PC set.
 */
void
TopologicalDiscovery<DataType, VarType>::symmetryCorrectPC(
  const VarType target,
  std::set<VarType>& cpc
) const
{
  auto initial = cpc;
  for (const VarType x: initial) {
    auto candidatesX = this->getCandidates(x);
    auto cpcX = this->getCandidatePC(x, std::move(candidatesX));
    if (cpcX.find(target) == cpcX.end()) {
      LOG_MESSAGE(info, "Symmetry Correction: Removing %s from the candidate PC of %s", this->m_data.varName(x), this->m_data.varName(target));
      cpc.erase(x);
    }
  }
}

template <typename DataType, typename VarType>
/**
 * @brief The top level function for getting the correct PC set
 *        for the given target variable.
 *
 * @param target The index of the target variable.
 *
 * @return A set containing the indices of all the variables
 *         in the correct PC set of the target variable.
 */
std::set<VarType>
TopologicalDiscovery<DataType, VarType>::getCorrectPC(
  const VarType target
) const
{
  auto candidates = this->getCandidates(target);
  auto cpc = this->getCandidatePC(target, std::move(candidates));
  this->symmetryCorrectPC(target, cpc);
  return cpc;
}

template <typename DataType, typename VarType>
/**
 * @brief The top level function for getting candidate MB for the given
 *        target variable, as per the algorithm proposed by Pena et al.
 *
 * @param target The index of the target variable.
 * @param candidates The indices of all the candidate variables.
 *
 * @return A set containing the indices of all the variables
 *         in the MB of the given target variable.
 */
std::set<VarType>
TopologicalDiscovery<DataType, VarType>::getCandidateMB(
  const VarType target,
  std::set<VarType> candidates
) const
{
  LOG_MESSAGE(info, "Getting MB for %s", this->m_data.varName(target));
  std::set<VarType> cmb;
  auto cpc = this->getCorrectPC(target);
  for (const VarType y: cpc) {
    LOG_MESSAGE(info, "Parent/Child: Adding %s to the candidate MB", this->m_data.varName(y));
    cmb.insert(y);
    auto cpcY = this->getCorrectPC(y);
    for (const VarType x: cpcY) {
      if ((x != target) && (cpc.find(x) == cpc.end())) {
        candidates.erase(x);
        LOG_MESSAGE(debug, "Checking %s for addition to MB", this->m_data.varName(x));
        auto ret = this->m_data.minAssocScoreSubset(target, x, candidates);
        if (this->m_data.isIndependent(ret.first)) {
          LOG_MESSAGE(debug, "%s found independent of the target, given a subset of the candidates", this->m_data.varName(x));
          auto& z = ret.second;
          z.insert(y);
          if (!this->m_data.isIndependent(target, x, z)) {
            LOG_MESSAGE(info, "Spouse: Adding %s to the candidate MB", this->m_data.varName(x));
            cmb.insert(x);
          }
        }
        candidates.insert(x);
      }
    }
  }
  return cmb;
}

template <typename DataType, typename VarType>
MMPC<DataType, VarType>::MMPC(
  const DataType& data
) : TopologicalDiscovery<DataType, VarType>(data)
{
}

template <typename DataType, typename VarType>
std::set<VarType>
MMPC<DataType, VarType>::getCandidatePC_impl(
  const VarType target,
  std::set<VarType> candidates
) const
{
  LOG_MESSAGE(info, "MMPC: Getting PC for %s", this->m_data.varName(target));
  std::set<VarType> cpc;
  bool changed = true;
  while ((candidates.size() > 0) && changed) {
    changed = false;
    // Find the variable which maximizes the minimum association score with the target,
    // given any subset of the current candidate PC
    VarType x;
    double scoreX = std::numeric_limits<double>::lowest();
    for (const VarType y: candidates) {
      LOG_MESSAGE(debug, "MMPC: Evaluating %s for the next candidate", this->m_data.varName(y));
      auto scoreY = this->m_data.minAssocScore(target, y, cpc);
      if (std::isless(scoreX, scoreY)) {
        x = y;
        scoreX = scoreY;
      }
    }
    LOG_MESSAGE(debug, "MMPC: %s chosen as the best candidate", this->m_data.varName(x));
    // Add the variable to the candidate PC if it is not
    // independedent of the target
    if (!this->m_data.isIndependent(scoreX)) {
      LOG_MESSAGE(info, "MMPC: Adding %s to the candidate PC", this->m_data.varName(x));
      cpc.insert(x);
      changed = true;
    }
    candidates.erase(x);
  }
  // Remove false positives from the candidate PC
  this->removeFalsePC(target, cpc);
  LOG_MESSAGE(info, "%s", std::string(60, '-'));
  return cpc;
}

template <typename DataType, typename VarType>
HITON<DataType, VarType>::HITON(
  const DataType& data
) : TopologicalDiscovery<DataType, VarType>(data)
{
}

template <typename DataType, typename VarType>
std::set<VarType>
HITON<DataType, VarType>::getCandidatePC_impl(
  const VarType target,
  std::set<VarType> candidates
) const
{
  LOG_MESSAGE(info, "HITON-PC: Getting PC for %s", this->m_data.varName(target));
  std::set<VarType> cpc;
  while (candidates.size() > 0) {
    // Find the variable which maximizes the marginal association score with the target
    VarType x;
    double scoreX = std::numeric_limits<double>::lowest();
    for (const VarType y: candidates) {
      LOG_MESSAGE(debug, "HITON-PC: Evaluating %s for the next candidate", this->m_data.varName(y));
      double scoreY = this->m_data.assocScore(target, y);
      if (std::isless(scoreX, scoreY)) {
        x = y;
        scoreX = scoreY;
      }
    }
    LOG_MESSAGE(debug, "HITON-PC: %s chosen as the best candidate", this->m_data.varName(x));
    // Add the variable to the candidate PC
    LOG_MESSAGE(info, "HITON-PC: Adding %s to the candidate PC of %s", this->m_data.varName(x), this->m_data.varName(target));
    cpc.insert(x);
    candidates.erase(x);
    // Remove false positives from the candidate PC
    this->removeFalsePC(target, cpc);
  }
  LOG_MESSAGE(info, "%s", std::string(60, '-'));
  return cpc;
}

template <typename DataType, typename VarType>
SemiInterleavedHITON<DataType, VarType>::SemiInterleavedHITON(
  const DataType& data
) : TopologicalDiscovery<DataType, VarType>(data)
{
}

template <typename DataType, typename VarType>
std::set<VarType>
SemiInterleavedHITON<DataType, VarType>::getCandidatePC_impl(
  const VarType target,
  std::set<VarType> candidates
) const
{
  LOG_MESSAGE(info, "SI-HITON-PC: Getting PC for %s", this->m_data.varName(target));
  std::set<VarType> cpc;
  while (candidates.size() > 0) {
    // Find the variable which maximizes the marginal association score with the target
    VarType x;
    double scoreX = std::numeric_limits<double>::lowest();
    std::set<VarType> remove;
    for (const VarType y: candidates) {
      LOG_MESSAGE(debug, "SI-HITON-PC: Evaluating %s for the next candidate", this->m_data.varName(y));
      double scoreY = this->m_data.assocScore(target, y);
      if (this->m_data.isIndependent(scoreY)) {
        LOG_MESSAGE(debug, "SI-HITON-PC: Marking %s for removal from the candidates", this->m_data.varName(y));
        // Can not be added to the candidate PC, mark for removal
        remove.insert(y);
        continue;
      }
      if (std::isless(scoreX, scoreY)) {
        x = y;
        scoreX = scoreY;
      }
    }
    // Remove all the candidates which can not be added
    for (const VarType y: remove) {
      candidates.erase(y);
    }
    if (candidates.empty()) {
      continue;
    }
    LOG_MESSAGE(debug, "SI-HITON-PC: %s chosen as the best candidate", this->m_data.varName(x));
    // Add the variable to the candidate PC
    LOG_MESSAGE(info, "SI-HITON-PC: Adding %s to the candidate PC of %s", this->m_data.varName(x), this->m_data.varName(target));
    cpc.insert(x);
    candidates.erase(x);
    // Remove false positives from the candidate PC
    this->removeFalsePC(target, cpc);
  }
  LOG_MESSAGE(info, "%s", std::string(60, '-'));
  return cpc;
}

template <typename DataType, typename VarType>
GetPC<DataType, VarType>::GetPC(
  const DataType& data
) : TopologicalDiscovery<DataType, VarType>(data)
{
}

template <typename DataType, typename VarType>
std::set<VarType>
GetPC<DataType, VarType>::getCandidatePC_impl(
  const VarType target,
  std::set<VarType> candidates
) const
{
  LOG_MESSAGE(info, "GetPC: Getting PC for %s", this->m_data.varName(target));
  std::set<VarType> cpc;
  bool changed = true;
  while ((candidates.size() > 0) && changed) {
    changed = false;
    // Find the variable which maximizes the minimum association score with the target,
    // given any subset of the current candidate PC
    VarType x;
    double scoreX = std::numeric_limits<double>::lowest();
    std::set<VarType> remove;
    for (const VarType y: candidates) {
      LOG_MESSAGE(debug, "GetPC: Evaluating %s for the next candidate", this->m_data.varName(y));
      auto scoreY = this->m_data.minAssocScore(target, y, cpc);
      if (this->m_data.isIndependent(scoreY)) {
        LOG_MESSAGE(debug, "GetPC: Marking %s for removal from the candidates", this->m_data.varName(y));
        // Can not be added to the candidate PC, mark for removal
        remove.insert(y);
        continue;
      }
      if (std::isless(scoreX, scoreY)) {
        x = y;
        scoreX = scoreY;
      }
    }
    // Remove all the candidates which can not be added
    for (const VarType y: remove) {
      candidates.erase(y);
    }
    if (candidates.empty()) {
      continue;
    }
    LOG_MESSAGE(debug, "GetPC: %s chosen as the best candidate", this->m_data.varName(x));
    // Add the variable to the candidate PC if it is not
    // independedent of the target
    if (!this->m_data.isIndependent(scoreX)) {
      LOG_MESSAGE(info, "GetPC: Adding %s to the candidate PC of %s", this->m_data.varName(x), this->m_data.varName(target));
      cpc.insert(x);
      changed = true;
    }
    candidates.erase(x);
    // Remove false positives from the candidate PC
    this->removeFalsePC(target, cpc);
  }
  LOG_MESSAGE(info, "%s", std::string(60, '-'));
  return cpc;
}

#endif // DETAIL_TOPOLOGICALDISCOVERY_HPP_
