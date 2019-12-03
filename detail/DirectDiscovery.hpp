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
  const mxx::comm& comm,
  const Data& data,
  const Var maxConditioning
) : ConstraintBasedDiscovery<Data, Var, Set>(comm, data, maxConditioning)
{
}

template <typename Data, typename Var, typename Set>
std::pair<Var, double>
DirectDiscovery<Data, Var, Set>::pickBestCandidate(
  const Var target,
  const Set& candidates,
  const Set& cmb
) const
{
  Var x = this->m_data.numVars();
  double scoreX = 0.0;
  for (const Var y: candidates) {
    LOG_MESSAGE(debug, "Grow: Evaluating %s for addition to the MB", this->m_data.varName(y));
    double scoreY = this->m_data.assocScore(target, y, cmb);
    if (std::isless(scoreX, scoreY)) {
      x = y;
      scoreX = scoreY;
    }
  }
  LOG_MESSAGE(debug, "Grow: %s chosen as the best candidate", this->m_data.varName(x));
  return std::make_pair(x, scoreX);
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
    LOG_MESSAGE(debug, "Shrink: Evaluating %s for removal from the MB of %s", this->m_data.varName(x), this->m_data.varName(target));
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
  LOG_MESSAGE(debug, "Neighbors: Evaluating %s for addition to the PC of %s", this->m_data.varName(y), this->m_data.varName(x));
  auto mbTest = mbX;
  // Pick the smaller of the two MBs
  if (mbY.size() > mbX.size()) {
    mbTest.erase(y);
  }
  else {
    mbTest = mbY;
    mbTest.erase(x);
  }
  return !this->m_data.isIndependentAnySubset(x, y, mbTest, this->m_maxConditioning);
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
  LOG_MESSAGE(info, "Neighbors: Getting PC from MB for %s", this->m_data.varName(target));
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
  const mxx::comm& comm,
  const Data& data,
  const Var maxConditioning
) : DirectDiscovery<Data, Var, Set>(comm, data, maxConditioning)
{
}

template <typename Data, typename Var, typename Set>
std::pair<Var, double>
GSMB<Data, Var, Set>::pickBestCandidate(
  const Var target,
  const Set& candidates,
  const Set& cmb
) const
{
  for (const Var y: candidates) {
    LOG_MESSAGE(debug, "Grow: Evaluating %s for addition to the MB", this->m_data.varName(y));
    double score = this->m_data.assocScore(target, y, cmb);
    if (!this->m_data.isIndependent(score)) {
      LOG_MESSAGE(debug, "Grow: %s chosen as the best candidate", this->m_data.varName(y));
      return std::make_pair(y, score);
    }
  }
  return std::make_pair(this->m_data.numVars(), 0.0);
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
  Var x = this->m_data.numVars();
  double scoreX = 0.0;
  while ((candidates.size() > 0) && changed) {
    changed = false;
    std::tie(x, scoreX) = this->pickBestCandidate(target, candidates, cmb);
    if (!this->m_data.isIndependent(scoreX)) {
      LOG_MESSAGE(info, "+ Adding %s to the MB of %s (score = %g)", this->m_data.varName(x), this->m_data.varName(target), scoreX);
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
IAMB<Data, Var, Set>::IAMB(
  const mxx::comm& comm,
  const Data& data,
  const Var maxConditioning
) : DirectDiscovery<Data, Var, Set>(comm, data, maxConditioning)
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
  Var x = this->m_data.numVars();
  double scoreX = 0.0;
  while ((candidates.size() > 0) && changed) {
    changed = false;
    std::tie(x, scoreX) = this->pickBestCandidate(target, candidates, cmb);
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
  const mxx::comm& comm,
  const Data& data,
  const Var maxConditioning
) : DirectDiscovery<Data, Var, Set>(comm, data, maxConditioning)
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
  Var x = this->m_data.numVars();
  double scoreX = 0.0;
  while ((candidates.size() > 0) && changed) {
    changed = false;
    std::tie(x, scoreX) = this->pickBestCandidate(target, candidates, cmb);
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
