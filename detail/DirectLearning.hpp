/**
 * @file DirectLearning.hpp
 * @brief Implementation of the classes for direct learning algorithms.
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
#ifndef DETAIL_DIRECTLEARNING_HPP_
#define DETAIL_DIRECTLEARNING_HPP_

#include "../SetUtils.hpp"

#include "utils/Logging.hpp"

#include <algorithm>


template <typename Data, typename Var, typename Set>
/**
 * @brief Constructs the object with the given data.
 *
 * @param data Reference to an object of the Data.
 */
DirectLearning<Data, Var, Set>::DirectLearning(
  const mxx::comm& comm,
  const Data& data,
  const Var maxConditioning
) : ConstraintBasedLearning<Data, Var, Set>(comm, data, maxConditioning)
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
DirectLearning<Data, Var, Set>::removeFalsePC(
  const Var target,
  Set& cpc
) const
{
  auto removed = set_init(Set(), this->m_data.numVars());
  auto initial = cpc;
  for (const Var x : initial) {
    cpc.erase(x);
    LOG_MESSAGE(debug, "False Positive: Testing %s for removal", this->m_data.varName(x));
    if (this->m_data.isIndependentAnySubset(target, x, cpc, this->m_maxConditioning)) {
      LOG_MESSAGE(info, "- Removing %s from the PC of %s (FP)", this->m_data.varName(x), this->m_data.varName(target));
      removed.insert(x);
    }
    cpc.insert(x);
  }
  cpc = set_difference(cpc, removed);
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
DirectLearning<Data, Var, Set>::getCandidateMB(
  const Var target,
  Set&& candidates
) const
{
  LOG_MESSAGE(info, "Blankets: Getting MB from PC for %s", this->m_data.varName(target));
  auto cmb = set_init(Set(), this->m_data.numVars());
  const auto& pc = this->getPC(target);
  for (const Var y : pc) {
    LOG_MESSAGE(info, "+ Adding %s to the MB of %s (parent/child)",
                      this->m_data.varName(y), this->m_data.varName(target));
    cmb.insert(y);
    const auto& pcY = this->getPC(y);
    for (const Var x : pcY) {
      if ((x != target) && !set_contains(pc, x)) {
        candidates.erase(x);
        LOG_MESSAGE(debug, "Evaluating %s for addition to the MB", this->m_data.varName(x));
        auto ret = this->m_data.maxPValueSubset(target, x, candidates, this->m_maxConditioning);
        if (this->m_data.isIndependent(ret.first)) {
          LOG_MESSAGE(debug, "%s found independent of the target, given a subset of the candidates", this->m_data.varName(x));
          auto& z = ret.second;
          z.insert(y);
          if (!this->m_data.isIndependent(target, x, z)) {
            LOG_MESSAGE(info, "+ Adding %s to the MB of %s (spouse)",
                              this->m_data.varName(x), this->m_data.varName(target));
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
/**
 * @brief Utility function for updating the maximum p-values of
 *        candidates for addition to the PC sets of the target.
 *
 * @param target The index of the target variable.
 * @param maxPValues The array containing the <max p-value, candidate> pairs.
 *                   Candidates found to be independent of the target after
 *                   the update are removed from the array.
 * @param cpc The current candidate PC set.
 * @param setNext Set containing the index of the next variable to be added
 *                to the candidate PC set (expected to contain exactly one element).
 */
void
DirectLearning<Data, Var, Set>::updateMaxPValues(
  const Var target,
  std::vector<std::pair<double, Var>>& maxPValues,
  const Set& cpc,
  const Set& setNext
) const
{
  for (auto& mpv : maxPValues) {
    auto pvY = mpv.first;
    auto y = mpv.second;
    LOG_MESSAGE(debug, "Updating max p-value for %s (previous p-value = %g)", this->m_data.varName(y), pvY);
    pvY = std::max(pvY, this->m_data.maxPValue(target, y, cpc, setNext, this->m_maxConditioning));
    LOG_MESSAGE(debug, "%s is ", this->m_data.isIndependent(pvY) ? "independent of" : "dependent on",
                       " the target %s (updated p-value = %g)",
                       this->m_data.varName(y), this->m_data.varName(target), pvY);
    mpv.first = pvY;
  }
  // Remove the independent elements to retain only the plausible candidates
  auto last = std::remove_if(maxPValues.begin(), maxPValues.end(), [this] (const std::pair<double, Var>& mpv)
                                                                          { return this->m_data.isIndependent(mpv.first); });
  maxPValues.erase(last, maxPValues.end());
}

template <typename Data, typename Var, typename Set>
MMPC<Data, Var, Set>::MMPC(
  const mxx::comm& comm,
  const Data& data,
  const Var maxConditioning
) : DirectLearning<Data, Var, Set>(comm, data, maxConditioning)
{
}

template <typename Data, typename Var, typename Set>
Set
MMPC<Data, Var, Set>::getCandidatePC(
  const Var target,
  Set&& candidates
) const
{
  LOG_MESSAGE(info, "%s", std::string(60, '-'));
  LOG_MESSAGE(info, "MMPC: Getting PC for %s", this->m_data.varName(target));
  std::vector<std::pair<double, Var>> maxPValues(candidates.size());
  std::transform(candidates.begin(), candidates.end(), maxPValues.begin(),
                 [] (const Var y) { return std::make_pair(0.0, y); });
  auto setNext = set_init(Set(), this->m_data.numVars());
  auto cpc = set_init(Set(), this->m_data.numVars());
  this->updateMaxPValues(target, maxPValues, cpc, setNext);
  while (!maxPValues.empty()) {
    // Choose the variable which minimizes the maximum p-value with the target,
    // given any subset of the current candidate PC
    auto mpv = std::min_element(maxPValues.begin(), maxPValues.end());
    auto x = mpv->second;
    LOG_MESSAGE(info, "+ Adding %s to the PC of %s (p-value = %g)",
                      this->m_data.varName(x), this->m_data.varName(target), mpv->first);
    setNext.insert(x);
    // Remove the min element by replacing it with the last element
    *mpv = maxPValues.back();
    maxPValues.pop_back();
    this->updateMaxPValues(target, maxPValues, cpc, setNext);
    cpc.insert(x);
    setNext.erase(x);
  }
  // Remove false positives from the candidate PC
  this->removeFalsePC(target, cpc);
  LOG_MESSAGE(info, "%s", std::string(60, '-'));
  return cpc;
}

template <typename Data, typename Var, typename Set>
HITON<Data, Var, Set>::HITON(
  const mxx::comm& comm,
  const Data& data,
  const Var maxConditioning
) : DirectLearning<Data, Var, Set>(comm, data, maxConditioning)
{
}

template <typename Data, typename Var, typename Set>
Set
HITON<Data, Var, Set>::getCandidatePC(
  const Var target,
  Set&& candidates
) const
{
  LOG_MESSAGE(info, "%s", std::string(60, '-'));
  LOG_MESSAGE(info, "HITON-PC: Getting PC for %s", this->m_data.varName(target));
  auto cpc = set_init(Set(), this->m_data.numVars());
  while (candidates.size() > 0) {
    // Find the variable which minimizes the marginal p-value with the target
    Var x = this->m_data.numVars();
    double pvX = std::numeric_limits<double>::max();
    for (const Var y : candidates) {
      LOG_MESSAGE(debug, "HITON-PC: Evaluating %s for addition to the PC", this->m_data.varName(y));
      double pvY = this->m_data.pValue(target, y);
      if (std::isgreater(pvX, pvY)) {
        x = y;
        pvX = pvY;
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
  const mxx::comm& comm,
  const Data& data,
  const Var maxConditioning
) : DirectLearning<Data, Var, Set>(comm, data, maxConditioning)
{
}

template <typename Data, typename Var, typename Set>
Set
SemiInterleavedHITON<Data, Var, Set>::getCandidatePC(
  const Var target,
  Set&& candidates
) const
{
  LOG_MESSAGE(info, "%s", std::string(60, '-'));
  LOG_MESSAGE(info, "SI-HITON-PC: Getting PC for %s", this->m_data.varName(target));
  std::vector<std::pair<double, Var>> maxPValues(candidates.size());
  std::transform(candidates.begin(), candidates.end(), maxPValues.begin(),
                 [] (const Var y) { return std::make_pair(0.0, y); });
  auto setEmpty = set_init(Set(), this->m_data.numVars());
  auto cpc = set_init(Set(), this->m_data.numVars());
  this->updateMaxPValues(target, maxPValues, cpc, setEmpty);
  // Sort the maxPValues array before iterating over it
  std::sort(maxPValues.begin(), maxPValues.end());
  for (const auto& mpv : maxPValues) {
    // Choose the variables in the ascending order of marginal p-value with the target
    // as a proxy for decreasing associativity
    auto x = mpv.second;
    LOG_MESSAGE(debug, "SI-HITON-PC: %s chosen as the best candidate", this->m_data.varName(x));
    if (!this->m_data.isIndependentAnySubset(target, x, cpc, this->m_maxConditioning)) {
      // Add the variable to the candidate PC
      LOG_MESSAGE(info, "+ Adding %s to the candidate PC of %s (p-value = %g)",
                        this->m_data.varName(x), this->m_data.varName(target), mpv.first);
      cpc.insert(x);
    }
  }
  // Remove false positives from the candidate PC
  this->removeFalsePC(target, cpc);
  LOG_MESSAGE(info, "%s", std::string(60, '-'));
  return cpc;
}

template <typename Data, typename Var, typename Set>
GetPC<Data, Var, Set>::GetPC(
  const mxx::comm& comm,
  const Data& data,
  const Var maxConditioning
) : DirectLearning<Data, Var, Set>(comm, data, maxConditioning)
{
}

template <typename Data, typename Var, typename Set>
Set
GetPC<Data, Var, Set>::getCandidatePC(
  const Var target,
  Set&& candidates
) const
{
  LOG_MESSAGE(info, "%s", std::string(60, '-'));
  LOG_MESSAGE(info, "GetPC: Getting PC for %s", this->m_data.varName(target));
  auto cpc = set_init(Set(), this->m_data.numVars());
  bool changed = true;
  while ((candidates.size() > 0) && changed) {
    changed = false;
    // Find the variable which minimizes the maximum p-value with the target,
    // given any subset of the current candidate PC
    Var x = this->m_data.numVars();
    double pvX = std::numeric_limits<double>::max();
    auto remove = set_init(Set(), this->m_data.numVars());
    for (const Var y : candidates) {
      LOG_MESSAGE(debug, "GetPC: Evaluating %s for addition to the PC", this->m_data.varName(y));
      auto pvY = this->m_data.maxPValue(target, y, cpc, this->m_maxConditioning);
      if (this->m_data.isIndependent(pvY)) {
        LOG_MESSAGE(debug, "GetPC: Marking %s for removal from the candidates", this->m_data.varName(y));
        // Can not be added to the candidate PC, mark for removal
        remove.insert(y);
        continue;
      }
      if (std::isgreater(pvX, pvY)) {
        x = y;
        pvX = pvY;
      }
    }
    // Remove all the candidates which can not be added
    for (const Var y : remove) {
      candidates.erase(y);
    }
    if (candidates.empty()) {
      continue;
    }
    LOG_MESSAGE(debug, "GetPC: %s chosen as the best candidate", this->m_data.varName(x));
    // Add the variable to the candidate PC if it is not
    // independedent of the target
    if (!this->m_data.isIndependent(pvX)) {
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

#endif // DETAIL_DIRECTLEARNING_HPP_
