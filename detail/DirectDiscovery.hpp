/**
 * @file DirectDiscovery.hpp
 * @brief Implementation of the DirectDiscovery as well as
 *        the derived class functions.
 */
#ifndef DETAIL_DIRECTDISCOVERY_HPP_
#define DETAIL_DIRECTDISCOVERY_HPP_

#include "../SetUtils.hpp"

#include "utils/Logging.hpp"

#include <numeric>


template <typename Data, typename Var, typename Set>
/**
 * @brief Constructs the object with the given data.
 *
 * @param data Reference to an object of the Data.
 */
DirectDiscovery<Data, Var, Set>::DirectDiscovery(
  const Data& data
) : ConstraintBasedDiscovery<Data, Var, Set>(data)
{
}

template <typename Data, typename Var, typename Set>
/**
 * @brief Function that shrinks the given candidate MB.
 *
 * @param target The index of the target variable.
 * @param cmb The indices of the variables in the candidate MB of the target variable.
 *            The function removes the indices from the candidate set.
 *
 * @return The indices of the variables that were removed from the candidate MB.
 */
Set
DirectDiscovery<Data, Var, Set>::shrinkMB(
  const Var target,
  Set& cmb
) const
{
  auto removed = set_init(Set(), this->m_data.numVars());
  if (cmb.empty()) {
    return removed;
  }
  auto initial = cmb;
  for (const Var x: initial) {
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

template <typename Data, typename Var, typename Set>
/**
 * @brief Evaluate if the second variable is in the PC of the first variable.
 *
 * @param x The first variable.
 * @param y The second variable.
 * @param mbX The Markov blanket of the first variable.
 * @param mbY The Markov blanket of the second variable.
 *
 * @return true if the second variable is in the PC of the first variable,
 *         otherwise return false.
 */
bool
DirectDiscovery<Data, Var, Set>::evaluateCandidatePC(
  const Var x,
  const Var y,
  const Set& mbX,
  const Set& mbY
) const
{
  LOG_MESSAGE(debug, "Direct Discovery: Evaluating %s for addition to the PC of %s", this->m_data.varName(y), this->m_data.varName(x));
  auto mbTest = mbX;
  // Pick the smaller of the two MBs
  if (mbY.size() > mbX.size()) {
    mbTest.erase(y);
  }
  else {
    mbTest = mbY;
    mbTest.erase(x);
  }
  return !this->m_data.isIndependentAnySubset(x, y, mbTest);
}

template <typename Data, typename Var, typename Set>
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
Set
DirectDiscovery<Data, Var, Set>::getCandidatePC(
  const Var target,
  Set candidates
) const
{
  LOG_MESSAGE(info, "Direct Discovery: Getting PC from MB for %s", this->m_data.varName(target));
  auto cpc = set_init(Set(), this->m_data.numVars());
  auto mb = this->getMB(target);
  for (const Var y: mb) {
    if (this->evaluateCandidatePC(target, y, mb, this->getMB(y))) {
      LOG_MESSAGE(info, "+ Adding %s to the PC of %s", this->m_data.varName(y), this->m_data.varName(target));
      cpc.insert(y);
    }
  }
  return cpc;
}

template <typename Data, typename Var, typename Set>
GSMB<Data, Var, Set>::GSMB(
  const Data& data
) : DirectDiscovery<Data, Var, Set>(data)
{
}

template <typename Data, typename Var, typename Set>
Set
GSMB<Data, Var, Set>::getCandidateMB(
  const Var target,
  Set candidates
) const
{
  LOG_MESSAGE(info, "%s", std::string(60, '-'));
  LOG_MESSAGE(info, "GSMB: Getting MB for %s", this->m_data.varName(target));
  auto cmb = set_init(Set(), this->m_data.numVars());
  bool changed = true;
  while ((candidates.size() > 0) && changed) {
    changed = false;
    auto curr = candidates.begin();
    while (!changed && (curr != candidates.end())) {
      auto x = *curr;
      LOG_MESSAGE(debug, "GSMB: Evaluating %s for addition to the MB", this->m_data.varName(x));
      auto score = this->m_data.assocScore(target, x, cmb);
      if (!this->m_data.isIndependent(score)) {
        LOG_MESSAGE(info, "+ Adding %s to the MB of %s (score = %g)", this->m_data.varName(x), this->m_data.varName(target), score);
        cmb.insert(x);
        candidates.erase(x);
        changed = true;
      }
      ++curr;
    }
  }
  this->shrinkMB(target, cmb);
  LOG_MESSAGE(info, "%s", std::string(60, '-'));
  return cmb;
}

template <typename Data, typename Var, typename Set>
IAMB<Data, Var, Set>::IAMB(
  const Data& data
) : DirectDiscovery<Data, Var, Set>(data)
{
}

template <typename Data, typename Var, typename Set>
Set
IAMB<Data, Var, Set>::getCandidateMB(
  const Var target,
  Set candidates
) const
{
  LOG_MESSAGE(info, "%s", std::string(60, '-'));
  LOG_MESSAGE(info, "IAMB: Getting MB for %s", this->m_data.varName(target));
  auto cmb = set_init(Set(), this->m_data.numVars());
  bool changed = true;
  while ((candidates.size() > 0) && changed) {
    changed = false;
    // Find the variable with the maximum association score with target,
    // given the current candidate MB
    Var x = this->m_data.numVars();
    double scoreX = 0.0;
    for (const Var y: candidates) {
      LOG_MESSAGE(debug, "IAMB: Evaluating %s for addition to the MB", this->m_data.varName(y));
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

template <typename Data, typename Var, typename Set>
InterIAMB<Data, Var, Set>::InterIAMB(
  const Data& data
) : DirectDiscovery<Data, Var, Set>(data)
{
}

template <typename Data, typename Var, typename Set>
Set
InterIAMB<Data, Var, Set>::getCandidateMB(
  const Var target,
  Set candidates
) const
{
  LOG_MESSAGE(info, "%s", std::string(60, '-'));
  LOG_MESSAGE(info, "InterIAMB: Getting MB for %s", this->m_data.varName(target));
  auto cmb = set_init(Set(), this->m_data.numVars());
  bool changed = true;
  while ((candidates.size() > 0) && changed) {
    changed = false;
    // Find the variable with the maximum association score with target,
    // given the current candidate MB
    Var x = this->m_data.numVars();
    double scoreX = 0.0;
    for (const Var y: candidates) {
      LOG_MESSAGE(debug, "InterIAMB: Evaluating %s for addition to the MB", this->m_data.varName(y));
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
      if (removed != Set{x}) {
        candidates = set_union(candidates, removed);
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
