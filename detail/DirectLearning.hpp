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

#include "common/SetUtils.hpp"
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
  const double alpha,
  const Var maxConditioning
) : LocalLearning<Data, Var, Set>(comm, data, alpha, maxConditioning),
    m_cachedCandidatePC()
{
  TIMER_RESET(m_tForward);
  TIMER_RESET(m_tBackward);
  TIMER_RESET(m_tDist);
  TIMER_RESET(m_tSymmetry);
  TIMER_RESET(m_tSync);
  TIMER_RESET(m_tNeighbors);
}

template <typename Data, typename Var, typename Set>
/**
 * @brief Default destructor that prints out timing information.
 */
DirectLearning<Data, Var, Set>::~DirectLearning(
)
{
  if (this->m_comm.is_first()) {
    TIMER_ELAPSED_NONZERO("Time taken in adding candidate neighbors: ", m_tForward);
    TIMER_ELAPSED_NONZERO("Time taken in removing false positive neighbors: ", m_tBackward);
    TIMER_ELAPSED_NONZERO("Time taken in redistributing: ", m_tDist);
    TIMER_ELAPSED_NONZERO("Time taken in symmetry correcting the neighbors: ", m_tSymmetry);
    TIMER_ELAPSED_NONZERO("Time taken in synchronizing the neighbors: ", m_tSync);
    TIMER_ELAPSED_NONZERO("Time taken in getting the neighbors: ", m_tNeighbors);
  }
}

template <typename Data, typename Var, typename Set>
/**
 * @brief The top level function for getting and caching
 *        the candidate PC sets.
 *
 * @param target The index of the target variable.
 * @param candidates The indices of all the candidate variables.
 *
 * @return A set containing the indices of all the variables
 *         in the candidate PC of the given target variable.
 */
Set
DirectLearning<Data, Var, Set>::getCandidatePC(
  const Var target,
  Set&& candidates
) const
{
  const auto cpc = this->getCandidatePC_impl(target, std::move(candidates));
  // Cache the candidate PCs, i.e., the PC sets before symmetry correction
  m_cachedCandidatePC[target] = cpc;
  return cpc;
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
    if (this->m_data.isIndependentAnySubset(this->m_alpha, target, x, cpc, this->m_maxConditioning)) {
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
std::set<std::pair<Var, Var>>
DirectLearning<Data, Var, Set>::backwardPhase(
  std::unordered_map<Var, Set>& myNeighbors
) const
{
  std::set<std::pair<Var, Var>> removed;
  // Remove false positive candidates from the PC sets
  // of all the local primary variables
  for (auto& neighbors : myNeighbors) {
    auto falsePC = this->removeFalsePC(neighbors.first, neighbors.second);
    for (const auto x : falsePC) {
      removed.insert(std::make_pair(neighbors.first, x));
    }
  }
  return removed;
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
    pvY = std::max(pvY, this->m_data.maxPValue(this->m_alpha, target, y, cpc, setNext, this->m_maxConditioning));
    LOG_MESSAGE(debug, "%s is " + std::string(this->m_data.isIndependent(this->m_alpha, pvY) ? "independent of" : "dependent on") +
                       " the target %s (updated p-value = %g)",
                       this->m_data.varName(y), this->m_data.varName(target), pvY);
    mpv.first = pvY;
  }
  // Remove the independent elements to retain only the plausible candidates
  auto last = std::remove_if(maxPValues.begin(), maxPValues.end(), [this] (const std::pair<double, Var>& mpv)
                                                                          { return this->m_data.isIndependent(this->m_alpha, mpv.first); });
  maxPValues.erase(last, maxPValues.end());
}

template <typename Data, typename Var, typename Set>
/**
 * @brief Utility function for updating the p-values in the local pairs
 *        candidates for addition to the PC sets of the primary variables.
 *
 * @param myPV A list of the p-values corresponding to all the local pairs.
 * @param myNeighbors A map with the candidate PCs of the primary variables on this processor.
 * @param nextAdditions A map containing the index of the next variable to be added
 *                      to the candidate PC set of every primary variable on this processor
 *                      (expected to contain exactly one element).
 */
void
DirectLearning<Data, Var, Set>::updateMyPValues(
  std::vector<std::tuple<Var, Var, double>>& myPV,
  const std::unordered_map<Var, Set>& myNeighbors,
  const std::unordered_map<Var, Set>& nextAdditions
) const
{
  for (auto& pv : myPV) {
    auto target = std::get<0>(pv);
    auto y = std::get<1>(pv);
    auto pvY = std::get<2>(pv);
    pvY = std::max(pvY, this->m_data.maxPValue(this->m_alpha, target, y, myNeighbors.at(target), nextAdditions.at(target), this->m_maxConditioning));
    LOG_MESSAGE(debug, "%s is " + std::string(this->m_data.isIndependent(this->m_alpha, pvY) ? "independent of" : "dependent on") +
                       " the target %s (updated p-value = %g)",
                       this->m_data.varName(y), this->m_data.varName(target), pvY);
    std::get<2>(pv) = pvY;
  }
  // Remove the independent elements to retain only the plausible candidates
  auto last = std::remove_if(myPV.begin(), myPV.end(), [this] (const std::tuple<Var, Var, double>& mpv)
                                                              { return this->m_data.isIndependent(this->m_alpha, std::get<2>(mpv)); });
  myPV.erase(last, myPV.end());
}

template <typename Data, typename Var, typename Set>
/**
 * @brief Function that grows candidate PCs of the primary variables on this processor.
 *
 * @param myPV A list of the p-values corresponding to all the local pairs.
 * @param comparePV A function for choosing the p-value using prefix scan.
 * @param reverse If the prefix scan should be performed in reverse as well.
 *                This is set to false if the first value in the segment should be picked.
 * @param myNeighbors A map with the candidate PCs of the primary variables on this processor.
 *
 * @return Set of all the tuples which were added to the candidate PCs on this processor.
 */
template <typename Compare>
std::set<std::tuple<Var, Var, double>>
DirectLearning<Data, Var, Set>::forwardPhase(
  const std::vector<std::tuple<Var, Var, double>>& myPV,
  const Compare& comparePV,
  const bool reverse,
  std::unordered_map<Var, Set>& myNeighbors
) const
{
  std::set<std::tuple<Var, Var, double>> added;
  std::vector<std::tuple<Var, Var, double>> minPV(myPV.size());
  // First, do a forward segmented parallel prefix with primary variable defining the segment boundaries
  // This will get the secondary variable with the minimum p-value to the corresponding primary variable boundary
  mxx::global_scan(myPV.begin(), myPV.end(), minPV.begin(), comparePV, false, this->m_comm);
  if (reverse) {
    // Then, do a reverse segmented parallel prefix with the same segments as before
    // This will effectively broadcast the secondary variable with the minimum p-value within the segments
    mxx::global_scan_inplace(minPV.rbegin(), minPV.rend(), comparePV, false, this->m_comm.reverse());
  }
  // There might be multiple local copies of the minimum p-value corresponding to every segment
  // Retain only one per segment
  auto comparePrimary = [] (const std::tuple<Var, Var, double>& a, const std::tuple<Var, Var, double>& b)
                           { return std::get<0>(a) == std::get<0>(b); };
  auto uniqueEnd = std::unique(minPV.begin(), minPV.end(), comparePrimary);
  for (auto it = minPV.begin(); it != uniqueEnd; ++it) {
    if (!this->m_data.isIndependent(this->m_alpha, std::get<2>(*it))) {
      // Add y to the blanket of x
      LOG_MESSAGE(info, "+ Adding %s to the PC of %s (p-value = %g)",
                  this->m_data.varName(std::get<1>(*it)), this->m_data.varName(std::get<0>(*it)), std::get<2>(*it));
      myNeighbors[std::get<0>(*it)].insert(std::get<1>(*it));
      added.insert(*it);
    }
  }
  return added;
}

template <typename Data, typename Var, typename Set>
/**
 * @brief Function for getting the undirected skeleton network in parallel.
 *
 * @param imbalanceThreshold Specifies the amount of imbalance the algorithm should tolerate.
 */
BayesianNetwork<Var>
DirectLearning<Data, Var, Set>::getSkeleton_parallel(
  const bool,
  const double imbalanceThreshold
) const
{
  TIMER_START(this->m_tNeighbors);
  std::vector<std::tuple<Var, Var, double>> myPV;
  std::unordered_map<Var, Set> myNeighbors;
  this->parallelInitialize(myPV, myNeighbors);
  // Remember all the local MB ordered pairs
  std::set<std::pair<Var, Var>> myAdded;
  this->forwardBackward(std::move(myPV), myNeighbors, myAdded, imbalanceThreshold);
  m_cachedCandidatePC = myNeighbors;

  /* Symmetry correction */
  this->m_comm.barrier();
  TIMER_START(this->m_tSymmetry);
  auto myPairs = this->symmetryCorrect(std::move(myNeighbors), std::move(myAdded));
  this->m_comm.barrier();
  TIMER_PAUSE(this->m_tSymmetry);
  /* End of Symmetry Correction */

  TIMER_START(this->m_tNeighbors);
  decltype(this->m_cachedPC) allNeighbors;
  for (const auto x : this->m_allVars) {
    allNeighbors[x] = set_init(Set(), this->m_data.numVars());
    // Initialize candidate sets for all the missing variables on this process
    if (m_cachedCandidatePC.find(x) == m_cachedCandidatePC.end()) {
      m_cachedCandidatePC[x] = set_init(Set(), this->m_data.numVars());
    }
  }
  for (const auto& p : myPairs) {
    allNeighbors[p.first].insert(p.second);
    allNeighbors[p.second].insert(p.first);
  }
  // Sync all the neighbors across all the processors
  TIMER_START(this->m_tSync);
  this->syncSets(m_cachedCandidatePC);
  this->syncSets(allNeighbors);
  TIMER_PAUSE(this->m_tSync);

  // We can now create the Bayesian network independently on every processor
  auto varNames = this->m_data.varNames(this->m_allVars);
  BayesianNetwork<Var> bn(varNames);
  for (const auto x : this->m_allVars) {
    for (const auto y : allNeighbors.at(x)) {
      if (x < y) {
        LOG_MESSAGE_IF(this->m_comm.is_first(), info, "+ Adding the edge %s <-> %s",
                       this->m_data.varName(x), this->m_data.varName(y));
        bn.addEdge(x, y, true);
      }
    }
  }
  // Correctly set the cached version of all the PC sets
  // All the sets have been symmetry corrected
  for (const auto x : allNeighbors) {
    this->m_cachedPCSymmetric[x.first] = true;
  }
  this->m_cachedPC = std::move(allNeighbors);
  return bn;
}

template <typename Data, typename Var, typename Set>
/**
 * @brief The top level function for getting the superset of MB
 *        for the given variable.
 *
 * @param target The index of the target variable.
 *
 * @return A set containing the indices of all the variables
 *         in the superset of MB of the given target variable.
 */
Set
DirectLearning<Data, Var, Set>::getMBSuperset(
  const Var target
) const
{
  // Use the non symmetry corrected PC sets for getting
  // the superset of MB sets for each variable
  // XXX: We are assuming that all the PCs would have
  //      already been learned by the time this gets called
  const auto& cpc = m_cachedCandidatePC.at(target);
  auto cmb = cpc;
  for (const Var y : cpc) {
    cmb = set_union(cmb, m_cachedCandidatePC.at(y));
  }
  cmb.erase(target);
  return cmb;
}

template <typename Data, typename Var, typename Set>
/**
 * @brief The top level function for getting the MB for the given
 *        target variable, using the PC sets, as per the algorithm proposed
 *        by Pena et al.
 *
 * @param target The index of the target variable.
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
  // Start with the superset of MB of the target,
  // if it is available
  if (m_cachedCandidatePC.find(target) != m_cachedCandidatePC.end()) {
    candidates = this->getMBSuperset(target);
  }
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
        auto ret = this->m_data.maxPValueSubset(this->m_alpha, target, x, candidates, this->m_maxConditioning);
        if (this->m_data.isIndependent(this->m_alpha, ret.first)) {
          LOG_MESSAGE(debug, "%s found independent of the target, given a subset of the candidates", this->m_data.varName(x));
          auto& z = ret.second;
          z.insert(y);
          if (!this->m_data.isIndependent(this->m_alpha, target, x, z)) {
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
 * @brief Checks if the given v-structure (y-x-z) is a collider.
 *        Also returns the p-value for the collider.
 */
std::pair<bool, double>
DirectLearning<Data, Var, Set>::checkCollider(
  const Var y,
  const Var x,
  const Var z
) const
{
  static auto smallerSet = [] (const Set& first, const Set& second) -> const auto&
                              { return (first.size() < second.size()) ? first : second; };
  auto setX = set_init(Set(), this->m_data.numVars());
  setX.insert(x);
  auto mbY = this->getMBSuperset(y);
  if (mbY.contains(z)) {
    mbY.erase(z);
  }
  mbY.erase(x);
  auto mbZ = this->getMBSuperset(z);
  if (mbZ.contains(y)) {
    mbZ.erase(y);
  }
  mbZ.erase(x);
  auto pv = this->m_data.maxPValue(this->m_alpha, y, z, smallerSet(mbY, mbZ), setX, this->m_maxConditioning);
  auto collider = !this->m_data.isIndependent(this->m_alpha, pv);
  return std::make_pair(collider, pv);
}

template <typename Data, typename Var, typename Set>
MMPC<Data, Var, Set>::MMPC(
  const mxx::comm& comm,
  const Data& data,
  const double alpha,
  const Var maxConditioning
) : DirectLearning<Data, Var, Set>(comm, data, alpha, maxConditioning)
{
}

template <typename Data, typename Var, typename Set>
Set
MMPC<Data, Var, Set>::getCandidatePC_impl(
  const Var target,
  Set&& candidates
) const
{
  LOG_MESSAGE(info, "%s", std::string(60, '-'));
  LOG_MESSAGE(info, "MMPC: Getting PC for %s", this->m_data.varName(target));
  TIMER_START(this->m_tForward);
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
  TIMER_PAUSE(this->m_tForward);
  // Remove false positives from the candidate PC
  TIMER_START(this->m_tBackward);
  this->removeFalsePC(target, cpc);
  TIMER_PAUSE(this->m_tBackward);
  LOG_MESSAGE(info, "%s", std::string(60, '-'));
  return cpc;
}

template <typename Data, typename Var, typename Set>
void
MMPC<Data, Var, Set>::forwardBackward(
  std::vector<std::tuple<Var, Var, double>>&& myPV,
  std::unordered_map<Var, Set>& myNeighbors,
  std::set<std::pair<Var, Var>>& myAdded,
  const double imbalanceThreshold
) const
{
  /* Forward Phase */
  TIMER_START(this->m_tForward);
  auto comparePV = [] (const std::tuple<Var, Var, double>& a, const std::tuple<Var, Var, double>& b)
                      { return (std::get<0>(a) == std::get<0>(b)) ?
                               (std::islessequal(std::get<2>(a), std::get<2>(b)) ? a : b) : b; };
  bool changed = true;
  auto nextAdditions = myNeighbors;
  for (auto& na : nextAdditions) {
    na.second.clear();
  }
  while (changed) {
    this->updateMyPValues(myPV, myNeighbors, nextAdditions);
    for (auto& na : nextAdditions) {
      if (!na.second.empty()) {
        myNeighbors[na.first] = set_union(myNeighbors.at(na.first), na.second);
        na.second.clear();
      }
    }
    auto added = this->forwardPhase(myPV, comparePV, true, nextAdditions);
    // Track candidate PC changes for all the primary variables across processors
    auto changes = set_init(Set(), this->m_data.numVars());
    for (const auto& as : added) {
      changes.insert(std::get<0>(as));
    }
    set_allunion(changes, this->m_comm);
    if (!changes.empty()) {
      if (!added.empty()) {
        // Record the added tuples belonging to this processor
        for (const auto& mpv : myPV) {
          if (added.find(mpv) != added.end()) {
            myAdded.insert(std::make_pair(std::get<0>(mpv), std::get<1>(mpv)));
          }
        }
      }
      // Remove added tuples from future consideration, if they were on this processor
      // Also remove the tuples corresponding to primary variables whose MB did not change
      auto addedOrUnchanged = [&added, &changes] (const std::tuple<Var, Var, double>& mpv)
                                                 { return (added.find(mpv) != added.end()) ||
                                                          !changes.contains(std::get<0>(mpv)); };
      auto newEnd = std::remove_if(myPV.begin(), myPV.end(), addedOrUnchanged);
      myPV.erase(newEnd, myPV.end());
      if (imbalanceThreshold > 1.0) {
        TIMER_START(this->m_tDist);
        if (this->fixImbalance(myPV, imbalanceThreshold)) {
          TIMER_START(this->m_tSync);
          this->syncMissingSets(myPV, myNeighbors);
          for (const auto& mn : myNeighbors) {
            if (nextAdditions.find(mn.first) == nextAdditions.end()) {
              nextAdditions.insert(std::make_pair(mn.first, set_init(Set(), this->m_data.numVars())));
            }
          }
          this->syncSets(nextAdditions);
          TIMER_PAUSE(this->m_tSync);
        }
        TIMER_PAUSE(this->m_tDist);
      }
    }
    else {
      changed = false;
    }
  }
  TIMER_START(this->m_tSync);
  this->syncSets(myNeighbors);
  TIMER_PAUSE(this->m_tSync);
  TIMER_PAUSE(this->m_tForward);
  /* End of Forward Phase */

  /* Backward Phase */
  TIMER_START(this->m_tBackward);
  this->backwardPhase(myNeighbors);
  TIMER_PAUSE(this->m_tBackward);
  /* End of Backward Phase */
}

template <typename Data, typename Var, typename Set>
HITON<Data, Var, Set>::HITON(
  const mxx::comm& comm,
  const Data& data,
  const double alpha,
  const Var maxConditioning
) : DirectLearning<Data, Var, Set>(comm, data, alpha, maxConditioning)
{
}

template <typename Data, typename Var, typename Set>
Set
HITON<Data, Var, Set>::getCandidatePC_impl(
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
BayesianNetwork<Var>
HITON<Data, Var, Set>::getSkeleton_parallel(
  const bool,
  const double
) const
{
  throw NotImplementedError("HITON: Parallel algorithm is not implemented yet");
  return BayesianNetwork<Var>(this->m_data.varNames(this->m_allVars));
}

template <typename Data, typename Var, typename Set>
SemiInterleavedHITON<Data, Var, Set>::SemiInterleavedHITON(
  const mxx::comm& comm,
  const Data& data,
  const double alpha,
  const Var maxConditioning
) : DirectLearning<Data, Var, Set>(comm, data, alpha, maxConditioning)
{
}

template <typename Data, typename Var, typename Set>
Set
SemiInterleavedHITON<Data, Var, Set>::getCandidatePC_impl(
  const Var target,
  Set&& candidates
) const
{
  LOG_MESSAGE(info, "%s", std::string(60, '-'));
  LOG_MESSAGE(info, "SI-HITON-PC: Getting PC for %s", this->m_data.varName(target));
  TIMER_START(this->m_tForward);
  std::vector<std::pair<double, Var>> maxPValues(candidates.size());
  std::transform(candidates.begin(), candidates.end(), maxPValues.begin(),
                 [] (const Var y) { return std::make_pair(0.0, y); });
  auto setNext = set_init(Set(), this->m_data.numVars());
  auto cpc = set_init(Set(), this->m_data.numVars());
  this->updateMaxPValues(target, maxPValues, cpc, setNext);
  // Sort the maxPValues array before iterating over it since we want to
  // consider the variables in the ascending order of marginal p-value with
  // the target as a proxy for decreasing associativity
  std::sort(maxPValues.begin(), maxPValues.end());
  while (!maxPValues.empty()) {
    // Choose the first variable which is dependent on the target,
    // given any subset of the current candidate PC
    auto mpv = maxPValues.begin();
    auto x = mpv->second;
    // Add the variable to the candidate PC
    LOG_MESSAGE(info, "+ Adding %s to the candidate PC of %s (p-value = %g)",
                      this->m_data.varName(x), this->m_data.varName(target), mpv->first);
    setNext.insert(x);
    // Erase the selected element
    maxPValues.erase(mpv);
    this->updateMaxPValues(target, maxPValues, cpc, setNext);
    cpc.insert(x);
    setNext.erase(x);
  }
  TIMER_PAUSE(this->m_tForward);
  // Remove false positives from the candidate PC
  TIMER_START(this->m_tBackward);
  this->removeFalsePC(target, cpc);
  TIMER_PAUSE(this->m_tBackward);
  LOG_MESSAGE(info, "%s", std::string(60, '-'));
  return cpc;
}

template <typename Data, typename Var, typename Set>
void
SemiInterleavedHITON<Data, Var, Set>::forwardBackward(
  std::vector<std::tuple<Var, Var, double>>&& myPV,
  std::unordered_map<Var, Set>& myNeighbors,
  std::set<std::pair<Var, Var>>& myAdded,
  const double imbalanceThreshold
) const
{
  /* Forward Phase */
  TIMER_START(this->m_tForward);
  bool changed = true;
  auto nextAdditions = myNeighbors;
  for (auto& na : nextAdditions) {
    na.second.clear();
  }
  // First, arrange the candidates for every primary variables
  // in the order of decreasing marginal associativity (or increasing p-values)
  // This is the order in which the sequential algorithm considers the candidates
  this->updateMyPValues(myPV, myNeighbors, nextAdditions);
  auto sortPV = [] (const std::tuple<Var, Var, double>& a, const std::tuple<Var, Var, double>& b)
                   { return (std::get<0>(a) == std::get<0>(b)) ?
                            (std::islessgreater(std::get<2>(a), std::get<2>(b)) ? std::isless(std::get<2>(a), std::get<2>(b)) : (std::get<1>(a) < std::get<1>(b))) :
                            (std::get<0>(a) < std::get<0>(b)); };
  mxx::sort(myPV.begin(), myPV.end(), sortPV, this->m_comm);
  auto comparePV = [] (const std::tuple<Var, Var, double>& a, const std::tuple<Var, Var, double>& b)
                      { return (std::get<0>(a) == std::get<0>(b)) ? a : b; };
  while (changed) {
    auto added = this->forwardPhase(myPV, comparePV, false, nextAdditions);
    // Track candidate PC changes for all the primary variables across processors
    auto changes = set_init(Set(), this->m_data.numVars());
    for (const auto& as : added) {
      changes.insert(std::get<0>(as));
    }
    set_allunion(changes, this->m_comm);
    if (!changes.empty()) {
      if (!added.empty()) {
        // Record the added tuples belonging to this processor
        for (const auto& mpv : myPV) {
          if (added.find(mpv) != added.end()) {
            myAdded.insert(std::make_pair(std::get<0>(mpv), std::get<1>(mpv)));
          }
        }
      }
      // Remove added tuples from future consideration, if they were on this processor
      // Also remove the tuples corresponding to primary variables whose MB did not change
      auto addedOrUnchanged = [&added, &changes] (const std::tuple<Var, Var, double>& mpv)
                                                 { return (added.find(mpv) != added.end()) ||
                                                          !changes.contains(std::get<0>(mpv)); };
      auto newEnd = std::remove_if(myPV.begin(), myPV.end(), addedOrUnchanged);
      myPV.erase(newEnd, myPV.end());
      if (imbalanceThreshold > 1.0) {
        TIMER_START(this->m_tDist);
        if (this->fixImbalance(myPV, imbalanceThreshold)) {
          TIMER_START(this->m_tSync);
          this->syncMissingSets(myPV, myNeighbors);
          for (const auto& mn : myNeighbors) {
            if (nextAdditions.find(mn.first) == nextAdditions.end()) {
              nextAdditions.insert(std::make_pair(mn.first, set_init(Set(), this->m_data.numVars())));
            }
          }
          this->syncSets(nextAdditions);
          TIMER_PAUSE(this->m_tSync);
        }
        TIMER_PAUSE(this->m_tDist);
      }
      this->updateMyPValues(myPV, myNeighbors, nextAdditions);
      for (auto& na : nextAdditions) {
        if (!na.second.empty()) {
          myNeighbors[na.first] = set_union(myNeighbors.at(na.first), na.second);
          na.second.clear();
        }
      }
    }
    else {
      changed = false;
    }
  }
  TIMER_START(this->m_tSync);
  this->syncSets(myNeighbors);
  TIMER_PAUSE(this->m_tSync);
  TIMER_PAUSE(this->m_tForward);
  /* End of Forward Phase */

  /* Backward Phase */
  TIMER_START(this->m_tBackward);
  this->backwardPhase(myNeighbors);
  TIMER_PAUSE(this->m_tBackward);
  /* End of Backward Phase */
}

template <typename Data, typename Var, typename Set>
GetPC<Data, Var, Set>::GetPC(
  const mxx::comm& comm,
  const Data& data,
  const double alpha,
  const Var maxConditioning
) : DirectLearning<Data, Var, Set>(comm, data, alpha, maxConditioning)
{
}

template <typename Data, typename Var, typename Set>
Set
GetPC<Data, Var, Set>::getCandidatePC_impl(
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
      auto pvY = this->m_data.maxPValue(this->m_alpha, target, y, cpc, this->m_maxConditioning);
      if (this->m_data.isIndependent(this->m_alpha, pvY)) {
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
    if (!this->m_data.isIndependent(this->m_alpha, pvX)) {
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

template <typename Data, typename Var, typename Set>
BayesianNetwork<Var>
GetPC<Data, Var, Set>::getSkeleton_parallel(
  const bool,
  const double
) const
{
  throw NotImplementedError("GetPC: Parallel algorithm is not implemented yet");
  return BayesianNetwork<Var>(this->m_data.varNames(this->m_allVars));
}

#endif // DETAIL_DIRECTLEARNING_HPP_
