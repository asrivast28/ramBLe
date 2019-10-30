/**
 * @file TopologicalDiscovery.hpp
 * @brief Implementation of the TopologicalDiscovery class and all the
 *        classes that are derived from it.
 */
#ifndef DETAIL_TOPOLOGICALDISCOVERY_HPP_
#define DETAIL_TOPOLOGICALDISCOVERY_HPP_

#include "../SetUtils.hpp"

#include "utils/Logging.hpp"

#include <algorithm>


template <typename Data, typename Var, typename Set>
/**
 * @brief Constructs the object with the given data.
 *
 * @param data Reference to an object of the Data.
 */
TopologicalDiscovery<Data, Var, Set>::TopologicalDiscovery(
  const Data& data,
  const mxx::comm& comm
) : ConstraintBasedDiscovery<Data, Var, Set>(data, comm)
{
}

template <typename Data, typename Var, typename Set>
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
Set
TopologicalDiscovery<Data, Var, Set>::removeFalsePC(
  const Var target,
  Set& cpc
) const
{
  auto removed = set_init(Set(), this->m_data.numVars());
  auto initial = cpc;
  for (const Var x: initial) {
    cpc.erase(x);
    LOG_MESSAGE(debug, "False Positive: Testing %s for removal", this->m_data.varName(x));
    if (this->m_data.isIndependentAnySubset(target, x, cpc)) {
      LOG_MESSAGE(info, "- Removing %s from the PC of %s (FP)", this->m_data.varName(x), this->m_data.varName(target));
      removed.insert(x);
    }
    else {
      cpc.insert(x);
    }
  }
  return removed;
}

template <typename Data, typename Var, typename Set>
/**
 * @brief The top level function for getting the candidate MB for the given
 *        target variable, using the PC sets, as per the algorithm proposed
 *        by Pena et al.
 *
 * @param target The index of the target variable.
 * @param candidates The indices of all the candidate variables.
 *
 * @return A set containing the indices of all the variables
 *         in the MB of the given target variable.
 */
Set
TopologicalDiscovery<Data, Var, Set>::getCandidateMB(
  const Var target,
  Set candidates
) const
{
  LOG_MESSAGE(info, "Topological Discovery: Getting MB from PC for %s", this->m_data.varName(target));
  auto cmb = set_init(Set(), this->m_data.numVars());
  auto pc = this->getPC(target);
  for (const Var y: pc) {
    LOG_MESSAGE(info, "+ Adding %s to the MB of %s (parent/child)", this->m_data.varName(y), this->m_data.varName(target));
    cmb.insert(y);
    auto pcY = this->getPC(y);
    for (const Var x: pcY) {
      if ((x != target) && !set_contains(pc, x)) {
        candidates.erase(x);
        LOG_MESSAGE(debug, "Evaluating %s for addition to the MB", this->m_data.varName(x));
        auto ret = this->m_data.minAssocScoreSubset(target, x, candidates);
        if (this->m_data.isIndependent(ret.first)) {
          LOG_MESSAGE(debug, "%s found independent of the target, given a subset of the candidates", this->m_data.varName(x));
          auto& z = ret.second;
          z.insert(y);
          if (!this->m_data.isIndependent(target, x, z)) {
            LOG_MESSAGE(info, "+ Adding %s to the MB of %s (spouse)", this->m_data.varName(x), this->m_data.varName(target));
            cmb.insert(x);
          }
        }
        candidates.insert(x);
      }
    }
  }
  return cmb;
}

template <typename Data, typename Var, typename Set>
MMPC<Data, Var, Set>::MMPC(
  const Data& data,
  const mxx::comm& comm
) : TopologicalDiscovery<Data, Var, Set>(data, comm)
{
}

template <typename Data, typename Var, typename Set>
Set
MMPC<Data, Var, Set>::getCandidatePC(
  const Var target,
  Set candidates
) const
{
  LOG_MESSAGE(info, "%s", std::string(60, '-'));
  LOG_MESSAGE(info, "MMPC: Getting PC for %s", this->m_data.varName(target));
  auto cpc = set_init(Set(), this->m_data.numVars());
  bool changed = true;
  while ((candidates.size() > 0) && changed) {
    changed = false;
    // Find the variable which maximizes the minimum association score with the target,
    // given any subset of the current candidate PC
    Var x = this->m_data.numVars();
    double scoreX = std::numeric_limits<double>::lowest();
    for (const Var y: candidates) {
      LOG_MESSAGE(debug, "MMPC: Evaluating %s for addition to the PC", this->m_data.varName(y));
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
      LOG_MESSAGE(info, "+ Adding %s to the PC of %s", this->m_data.varName(x), this->m_data.varName(target));
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

template <typename Data, typename Var, typename Set>
HITON<Data, Var, Set>::HITON(
  const Data& data,
  const mxx::comm& comm
) : TopologicalDiscovery<Data, Var, Set>(data, comm)
{
}

template <typename Data, typename Var, typename Set>
Set
HITON<Data, Var, Set>::getCandidatePC(
  const Var target,
  Set candidates
) const
{
  LOG_MESSAGE(info, "%s", std::string(60, '-'));
  LOG_MESSAGE(info, "HITON-PC: Getting PC for %s", this->m_data.varName(target));
  auto cpc = set_init(Set(), this->m_data.numVars());
  while (candidates.size() > 0) {
    // Find the variable which maximizes the marginal association score with the target
    Var x = this->m_data.numVars();
    double scoreX = std::numeric_limits<double>::lowest();
    for (const Var y: candidates) {
      LOG_MESSAGE(debug, "HITON-PC: Evaluating %s for addition to the PC", this->m_data.varName(y));
      double scoreY = this->m_data.assocScore(target, y);
      if (std::isless(scoreX, scoreY)) {
        x = y;
        scoreX = scoreY;
      }
    }
    LOG_MESSAGE(debug, "HITON-PC: %s chosen as the best candidate", this->m_data.varName(x));
    // Add the variable to the candidate PC
    LOG_MESSAGE(info, "+ Adding %s to the PC of %s", this->m_data.varName(x), this->m_data.varName(target));
    cpc.insert(x);
    candidates.erase(x);
    // Remove false positives from the candidate PC
    this->removeFalsePC(target, cpc);
  }
  LOG_MESSAGE(info, "%s", std::string(60, '-'));
  return cpc;
}

template <typename Data, typename Var, typename Set>
SemiInterleavedHITON<Data, Var, Set>::SemiInterleavedHITON(
  const Data& data,
  const mxx::comm& comm
) : TopologicalDiscovery<Data, Var, Set>(data, comm)
{
}

template <typename Data, typename Var, typename Set>
Set
SemiInterleavedHITON<Data, Var, Set>::getCandidatePC(
  const Var target,
  Set candidates
) const
{
  LOG_MESSAGE(info, "%s", std::string(60, '-'));
  LOG_MESSAGE(info, "SI-HITON-PC: Getting PC for %s", this->m_data.varName(target));
  auto cpc = set_init(Set(), this->m_data.numVars());
  while (candidates.size() > 0) {
    // Find the variable which maximizes the marginal association score with the target
    Var x = this->m_data.numVars();
    double scoreX = std::numeric_limits<double>::lowest();
    auto remove = set_init(Set(), this->m_data.numVars());
    for (const Var y: candidates) {
      LOG_MESSAGE(debug, "SI-HITON-PC: Evaluating %s for addition to the PC", this->m_data.varName(y));
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
    for (const Var y: remove) {
      candidates.erase(y);
    }
    if (candidates.empty()) {
      continue;
    }
    LOG_MESSAGE(debug, "SI-HITON-PC: %s chosen as the best candidate", this->m_data.varName(x));
    // Add the variable to the candidate PC
    LOG_MESSAGE(info, "+ Adding %s to the PC of %s", this->m_data.varName(x), this->m_data.varName(target));
    cpc.insert(x);
    candidates.erase(x);
    // Remove false positives from the candidate PC
    this->removeFalsePC(target, cpc);
  }
  LOG_MESSAGE(info, "%s", std::string(60, '-'));
  return cpc;
}

template <typename Data, typename Var, typename Set>
GetPC<Data, Var, Set>::GetPC(
  const Data& data,
  const mxx::comm& comm
) : TopologicalDiscovery<Data, Var, Set>(data, comm)
{
}

template <typename Data, typename Var, typename Set>
Set
GetPC<Data, Var, Set>::getCandidatePC(
  const Var target,
  Set candidates
) const
{
  LOG_MESSAGE(info, "%s", std::string(60, '-'));
  LOG_MESSAGE(info, "GetPC: Getting PC for %s", this->m_data.varName(target));
  auto cpc = set_init(Set(), this->m_data.numVars());
  bool changed = true;
  while ((candidates.size() > 0) && changed) {
    changed = false;
    // Find the variable which maximizes the minimum association score with the target,
    // given any subset of the current candidate PC
    Var x = this->m_data.numVars();
    double scoreX = std::numeric_limits<double>::lowest();
    auto remove = set_init(Set(), this->m_data.numVars());
    for (const Var y: candidates) {
      LOG_MESSAGE(debug, "GetPC: Evaluating %s for addition to the PC", this->m_data.varName(y));
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
    for (const Var y: remove) {
      candidates.erase(y);
    }
    if (candidates.empty()) {
      continue;
    }
    LOG_MESSAGE(debug, "GetPC: %s chosen as the best candidate", this->m_data.varName(x));
    // Add the variable to the candidate PC if it is not
    // independedent of the target
    if (!this->m_data.isIndependent(scoreX)) {
      LOG_MESSAGE(info, "+ Adding %s to the PC of %s", this->m_data.varName(x), this->m_data.varName(target));
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
