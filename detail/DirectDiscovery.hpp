/**
 * @file DirectDiscovery.hpp
 * @brief Implementation of the DirectDiscovery as well as
 *        the derived class functions.
 */
#ifndef DETAIL_DIRECTDISCOVERY_HPP_
#define DETAIL_DIRECTDISCOVERY_HPP_

#include "utils/Logging.hpp"

#include <numeric>


template <typename DataType, typename VarType>
/**
 * @brief Constructs the object with the given data.
 *
 * @param data Reference to an object of the DataType.
 */
DirectDiscovery<DataType, VarType>::DirectDiscovery(
  const DataType& data
) : MBDiscovery<DataType, VarType>(data)
{
}

template <typename DataType, typename VarType>
/**
 * @brief Function that shrinks the given candidate MB.
 *
 * @param target The index of the target variable.
 * @param cmb The indices of the variables in the candidate MB of the target variable.
 *                The function removes the indices from the candidate set.
 *
 * @return The indices of the variables that were removed from the candidate MB.
 */
std::set<VarType>
DirectDiscovery<DataType, VarType>::shrinkMB(
  const VarType target,
  std::set<VarType>& cmb
) const
{
  std::set<VarType> removed;
  if (cmb.empty()) {
    return removed;
  }
  auto initial = cmb;
  for (const VarType x: initial) {
    cmb.erase(x);
    LOG_MESSAGE(debug, "Shrink Phase: Evaluating %s for removal from the MB of %s", this->m_data.varName(x), this->m_data.varName(target));
    if (this->m_data.isIndependent(target, x, cmb)) {
      LOG_MESSAGE(info, "- Removing %s from the MB of %s (shrink)", this->m_data.varName(x), this->m_data.varName(target));
      removed.insert(x);
    }
    else {
      cmb.insert(x);
    }
  }
  return removed;
}

template <typename DataType, typename VarType>
/**
 * @brief The top level function for getting the candidate PC for the given
 *        target variable, using the MB of the variable.
 *
 * @param target The index of the target variable.
 * @param candidates The indices of all the candidate variables.
 *
 * @return A set containing the indices of all the variables
 *         in the PC of the given target variable.
 */
std::set<VarType>
DirectDiscovery<DataType, VarType>::getCandidatePC(
  const VarType target,
  std::set<VarType> candidates
) const
{
  LOG_MESSAGE(info, "Direct Discovery: Getting PC from MB for %s", this->m_data.varName(target));
  std::set<VarType> cpc;
  auto mb = this->getMB(target);
  for (const VarType y: mb) {
    LOG_MESSAGE(debug, "Checking %s for addition to MB", this->m_data.varName(y));
    auto mbTest = mb;
    auto mbY = this->getMB(y);
    // Pick the smaller of the two MBs
    if (mbY.size() > mb.size()) {
      mbTest.erase(y);
    }
    else {
      mbTest = mbY;
      mbTest.erase(target);
    }
    if (!this->m_data.isIndependentAnySubset(target, y, mbTest)) {
      LOG_MESSAGE(info, "+ Adding %s to the PC of %s", this->m_data.varName(y), this->m_data.varName(target));
      cpc.insert(y);
    }
  }
  return cpc;
}

template <typename DataType, typename VarType>
GSMB<DataType, VarType>::GSMB(
  const DataType& data
) : DirectDiscovery<DataType, VarType>(data)
{
}

template <typename DataType, typename VarType>
std::set<VarType>
GSMB<DataType, VarType>::getCandidateMB(
  const VarType target,
  std::set<VarType> candidates
) const
{
  LOG_MESSAGE(info, "%s", std::string(60, '-'));
  LOG_MESSAGE(info, "GSMB: Getting MB for %s", this->m_data.varName(target));
  std::set<VarType> cmb;
  bool changed = true;
  while ((candidates.size() > 0) && changed) {
    changed = false;
    auto thisCandidates = candidates;
    while (thisCandidates.size() > 0) {
      // Find the variable with the maximum marginal
      // association score with target
      VarType x;
      double scoreX = 0.0;
      for (const VarType y: thisCandidates) {
        LOG_MESSAGE(debug, "GSMB: Evaluating %s for the next candidate", this->m_data.varName(y));
        double scoreY = this->m_data.assocScore(target, y);
        if (std::isless(scoreX, scoreY)) {
          x = y;
          scoreX = scoreY;
        }
      }
      LOG_MESSAGE(debug, "GSMB: %s chosen as the best candidate", this->m_data.varName(x));
      // Add the variable to the candidate MB if it is not
      // independedent of the target, given the current MB
      if (!this->m_data.isIndependent(target, x, cmb)) {
        LOG_MESSAGE(info, "+ Adding %s to the MB of %s", this->m_data.varName(x), this->m_data.varName(target));
        cmb.insert(x);
        candidates.erase(x);
        changed = true;
      }
      thisCandidates.erase(x);
    }
  }
  this->shrinkMB(target, cmb);
  LOG_MESSAGE(info, "%s", std::string(60, '-'));
  return cmb;
}

template <typename DataType, typename VarType>
IAMB<DataType, VarType>::IAMB(
  const DataType& data
) : DirectDiscovery<DataType, VarType>(data)
{
}

template <typename DataType, typename VarType>
std::set<VarType>
IAMB<DataType, VarType>::getCandidateMB(
  const VarType target,
  std::set<VarType> candidates
) const
{
  LOG_MESSAGE(info, "%s", std::string(60, '-'));
  LOG_MESSAGE(info, "IAMB: Getting MB for %s", this->m_data.varName(target));
  std::set<VarType> cmb;
  bool changed = true;
  while ((candidates.size() > 0) && changed) {
    changed = false;
    // Find the variable with the maximum association score with target,
    // given the current candidate MB
    VarType x;
    double scoreX = 0.0;
    for (const VarType y: candidates) {
      LOG_MESSAGE(debug, "IAMB: Evaluating %s for the next candidate", this->m_data.varName(y));
      double scoreY = this->m_data.assocScore(target, y, cmb);
      if (std::isless(scoreX, scoreY)) {
        x = y;
        scoreX = scoreY;
      }
    }
    LOG_MESSAGE(debug, "IAMB: %s chosen as the best candidate", this->m_data.varName(x));
    // Add the variable to the candidate MB if it is not
    // independedent of the target
    if (!this->m_data.isIndependent(scoreX)) {
      LOG_MESSAGE(info, "+ Adding %s to the MB of %s", this->m_data.varName(x), this->m_data.varName(target));
      cmb.insert(x);
      candidates.erase(x);
      changed = true;
    }
  }
  this->shrinkMB(target, cmb);
  LOG_MESSAGE(info, "%s", std::string(60, '-'));
  return cmb;
}

template <typename DataType, typename VarType>
InterIAMB<DataType, VarType>::InterIAMB(
  const DataType& data
) : DirectDiscovery<DataType, VarType>(data)
{
}

template <typename DataType, typename VarType>
std::set<VarType>
InterIAMB<DataType, VarType>::getCandidateMB(
  const VarType target,
  std::set<VarType> candidates
) const
{
  LOG_MESSAGE(info, "%s", std::string(60, '-'));
  LOG_MESSAGE(info, "InterIAMB: Getting MB for %s", this->m_data.varName(target));
  std::set<VarType> cmb;
  bool changed = true;
  while ((candidates.size() > 0) && changed) {
    changed = false;
    // Find the variable with the maximum association score with target,
    // given the current candidate MB
    VarType x;
    double scoreX = 0.0;
    for (const VarType y: candidates) {
      LOG_MESSAGE(debug, "InterIAMB: Evaluating %s for the next candidate", this->m_data.varName(y));
      double scoreY = this->m_data.assocScore(target, y, cmb);
      if (std::isless(scoreX, scoreY)) {
        x = y;
        scoreX = scoreY;
      }
    }
    LOG_MESSAGE(debug, "InterIAMB: %s chosen as the best candidate", this->m_data.varName(x));
    // Add the variable to the candidate MB if it is not
    // independedent of the target
    if (!this->m_data.isIndependent(scoreX)) {
      LOG_MESSAGE(info, "+ Adding %s to the MB of %s", this->m_data.varName(x), this->m_data.varName(target));
      cmb.insert(x);
      candidates.erase(x);
      changed = true;
    }
    if (changed) {
      auto removed = this->shrinkMB(target, cmb);
      // If the last added variable was removed then the candidate MB
      // did not really change
      if (removed != std::set<VarType>{x}) {
        std::set<VarType> temp;
        std::set_union(candidates.begin(), candidates.end(), removed.begin(), removed.end(), std::inserter(temp, temp.begin()));
        candidates = temp;
      }
      else {
        changed = false;
      }
    }
  }
  LOG_MESSAGE(info, "%s", std::string(60, '-'));
  return cmb;
}

#endif // DETAIL_DIRECTDISCOVERY_HPP_
