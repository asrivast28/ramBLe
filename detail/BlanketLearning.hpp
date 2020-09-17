/**
 * @file BlanketLearning.hpp
 * @brief Implementation of the classes for blanket learning algorithms.
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
#ifndef DETAIL_BLANKETLEARNING_HPP_
#define DETAIL_BLANKETLEARNING_HPP_

#include "../SetUtils.hpp"

#include "utils/Logging.hpp"

#include <numeric>


template <typename Data, typename Var, typename Set>
/**
 * @brief Constructs the object with the given data.
 *
 * @param data Reference to an object of the Data.
 */
BlanketLearning<Data, Var, Set>::BlanketLearning(
  const mxx::comm& comm,
  const Data& data,
  const Var maxConditioning
) : ConstraintBasedLearning<Data, Var, Set>(comm, data, maxConditioning)
{
  TIMER_RESET(m_tGrow);
  TIMER_RESET(m_tShrink);
  TIMER_RESET(m_tDist);
  TIMER_RESET(m_tSymmetry);
  TIMER_RESET(m_tSync);
  TIMER_RESET(m_tBlankets);
  TIMER_RESET(m_tNeighbors);
}

template <typename Data, typename Var, typename Set>
BlanketLearning<Data, Var, Set>::~BlanketLearning(
)
{
  if (this->m_comm.is_first()) {
    TIMER_ELAPSED_NONZERO("Time taken in growing the candidate blankets: ", m_tGrow);
    TIMER_ELAPSED_NONZERO("Time taken in shrinking the candidate blankets: ", m_tShrink);
    TIMER_ELAPSED_NONZERO("Time taken in redistributing: ", m_tDist);
    TIMER_ELAPSED_NONZERO("Time taken in symmetry correcting the blankets: ", m_tSymmetry);
    TIMER_ELAPSED_NONZERO("Time taken in synchronizing the blankets: ", m_tSync);
    TIMER_ELAPSED_NONZERO("Time taken in getting the blankets: ", m_tBlankets);
    TIMER_ELAPSED_NONZERO("Time taken in getting the neighbors: ", m_tNeighbors);
  }
}

template <typename Data, typename Var, typename Set>
std::pair<Var, double>
BlanketLearning<Data, Var, Set>::pickBestCandidate(
  const Var target,
  const Set& candidates,
  const Set& cmb
) const
{
  Var x = this->m_data.numVars();
  double pvX = std::numeric_limits<double>::max();
  for (const Var y : candidates) {
    LOG_MESSAGE(debug, "Grow: Evaluating %s for addition to the MB", this->m_data.varName(y));
    double pvY = this->m_data.pValue(target, y, cmb);
    if (std::isgreater(pvX, pvY)) {
      x = y;
      pvX = pvY;
    }
  }
  LOG_MESSAGE_IF(x < this->m_data.numVars(), debug, "Grow: %s chosen as the best candidate", this->m_data.varName(x));
  return std::make_pair(x, pvX);
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
BlanketLearning<Data, Var, Set>::shrinkMB(
  const Var target,
  Set& cmb
) const
{
  auto removed = set_init(Set(), this->m_data.numVars());
  if (cmb.empty()) {
    return removed;
  }
  auto initial = cmb;
  for (const Var x : initial) {
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
BlanketLearning<Data, Var, Set>::evaluateCandidatePC(
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

  return !this->m_data.isIndependentAnySubset(x, y, mbTest, this->m_maxConditioning, this->m_comm);
}

template <typename Data, typename Var, typename Set>
/**
 * @brief The top level function for getting the candidate PC for the given
 *        target variable, using the MB of the variable.
 *
 * @param target The index of the target variable.
 *
 * @return A set containing the indices of all the variables
 *         in the PC of the given target variable.
 */
Set
BlanketLearning<Data, Var, Set>::getCandidatePC(
  const Var target,
  Set&&
) const
{
  LOG_MESSAGE(info, "Neighbors: Getting PC from MB for %s",
              this->m_data.varName(target));
  auto cpc = set_init(Set(), this->m_data.numVars());
  const auto& mb = this->getMB(target);
  for (const Var y : mb) {
    if (this->evaluateCandidatePC(target, y, mb, this->getMB(y))) {
      LOG_MESSAGE(info, "+ Adding %s to the PC of %s",
                  this->m_data.varName(y), this->m_data.varName(target));
      cpc.insert(y);
    }
  }
  return cpc;
}

template <typename Data, typename Var, typename Set>
/**
 * @brief Function that updates p-values for all the local pairs, given the current candidate MBs.
 *
 * @param myPV A list of the p-values corresponding to all the local pairs.
 * @param myBlankets A map with the candidate MBs of the primary variables on this processor.
 */
void
BlanketLearning<Data, Var, Set>::updatePValues(
  std::vector<std::tuple<Var, Var, double>>& myPV,
  const std::unordered_map<Var, Set>& myBlankets
) const
{
  for (auto& pv : myPV) {
    LOG_MESSAGE(debug, "Updating the p-value for the pair (%s, %s)",
                this->m_data.varName(std::get<0>(pv)), this->m_data.varName(std::get<1>(pv)));
    std::get<2>(pv) = this->m_data.pValue(std::get<0>(pv),
                                          std::get<1>(pv),
                                          myBlankets.at(std::get<0>(pv)));
  }
}

template <typename Data, typename Var, typename Set>
/**
 * @brief Function that grows candidate MBs of the primary variables on this processor.
 *
 * @param myPV A list of the p-values corresponding to all the local pairs.
 * @param myBlankets A map with the candidate MBs of the primary variables on this processor.
 *
 * @return Set of all the tuples which were added to the candidate MBs on this processor.
 */
std::set<std::tuple<Var, Var, double>>
BlanketLearning<Data, Var, Set>::growAll(
  const std::vector<std::tuple<Var, Var, double>>& myPV,
  std::unordered_map<Var, Set>& myBlankets
) const
{
  std::set<std::tuple<Var, Var, double>> added;
  std::vector<std::tuple<Var, Var, double>> minPV(myPV.size());
  auto comparePV = [] (const std::tuple<Var, Var, double>& a, const std::tuple<Var, Var, double>& b)
                      { return (std::get<0>(a) == std::get<0>(b)) ?
                               (std::islessequal(std::get<2>(a), std::get<2>(b)) ? a : b) : b; };
  // First, do a forward segmented parallel prefix with primary variable defining the segment boundaries
  // This will get the secondary variable with the minimum p-value to the corresponding primary variable boundary
  mxx::global_scan(myPV.begin(), myPV.end(), minPV.begin(), comparePV, false, this->m_comm);
  // Then, do a reverse segmented parallel prefix with the same segments as before
  // This will effectively broadcast the secondary variable with the minimum p-value within the segments
  mxx::global_scan_inplace(minPV.rbegin(), minPV.rend(), comparePV, false, this->m_comm.reverse());
  // There might be multiple local copies of the minimum p-value corresponding to every segment
  // Retain only one per segment
  auto comparePrimary = [] (const std::tuple<Var, Var, double>& a, const std::tuple<Var, Var, double>& b)
                           { return std::get<0>(a) == std::get<0>(b); };
  auto uniqueEnd = std::unique(minPV.begin(), minPV.end(), comparePrimary);
  for (auto it = minPV.begin(); it != uniqueEnd; ++it) {
    if (!this->m_data.isIndependent(std::get<2>(*it))) {
      // Add y to the blanket of x
      LOG_MESSAGE(info, "%d: + Adding %s to the MB of %s (p-value = %g)", this->m_comm.rank(),
                  this->m_data.varName(std::get<1>(*it)), this->m_data.varName(std::get<0>(*it)), std::get<2>(*it));
      myBlankets[std::get<0>(*it)].insert(std::get<1>(*it));
      added.insert(*it);
    }
  }
  return added;
}

template <typename Data, typename Var, typename Set>
/**
 * @brief Function that shrinks candidate MBs of the primary variables on this processor.
 *
 * @param myBlankets A map with the candidate MBs of the primary variables on this processor.
 *
 * @return Set of all the pairs which were removed from the candidate MBs on this processor.
 */
std::set<std::pair<Var, Var>>
BlanketLearning<Data, Var, Set>::shrinkAll(
  std::unordered_map<Var, Set>& myBlankets
) const
{
  std::set<std::pair<Var, Var>> removed;
  // Call shrink for all the locally stored blankets
  for (auto& mb : myBlankets) {
    auto shrunk = this->shrinkMB(mb.first, mb.second);
    for (const auto x : shrunk) {
      removed.insert(std::make_pair(mb.first, x));
    }
  }
  return removed;
}

template <typename Data, typename Var, typename Set>
/**
 * @brief Function that performs grow-shrink for multiple candidate MBs.
 *
 * @param myPV A list of the p-values corresponding to all the local pairs.
 * @param myBlankets A map with all the local candidate MBs.
 * @param myAdded A list of all the local candidate MB pairs.
 * @param imbalanceThreshold Specifies the amount of imbalance the algorithm should tolerate.
 *
 * @return The total size of all the blankets at the end.
 */
void
BlanketLearning<Data, Var, Set>::growShrink(
  std::vector<std::tuple<Var, Var, double>>&& myPV,
  std::unordered_map<Var, Set>& myBlankets,
  std::set<std::pair<Var, Var>>& myAdded,
  const double imbalanceThreshold
) const
{
  /* Grow Phase */
  TIMER_START(this->m_tGrow);
  bool changed = true;
  while (changed) {
    this->updatePValues(myPV, myBlankets);
    auto added = this->growAll(myPV, myBlankets);
    // Track blanket changes for all the primary variables across processors
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
          this->syncMissingSets(myPV, myBlankets);
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
  this->syncSets(myBlankets);
  TIMER_PAUSE(this->m_tSync);
  TIMER_PAUSE(this->m_tGrow);
  if (this->m_comm.is_first()) {
  }
  /* End of Grow Phase */

  /* Shrink Phase */
  TIMER_START(this->m_tShrink);
  this->shrinkAll(myBlankets);
  TIMER_PAUSE(this->m_tShrink);
  /* End of Shrink Phase */
}

template <typename Data, typename Var, typename Set>
/**
 * @brief Function for getting the undirected skeleton network in parallel.
 *
 * @param imbalanceThreshold Specifies the amount of imbalance the algorithm should tolerate.
 */
BayesianNetwork<Var>
BlanketLearning<Data, Var, Set>::getSkeleton_parallel(
  const double imbalanceThreshold,
  std::unordered_map<Var, Set>& allBlankets,
  std::unordered_map<Var, Set>& allNeighbors
) const
{
  TIMER_START(this->m_tBlankets);
  std::vector<std::tuple<Var, Var, double>> myPV;
  std::unordered_map<Var, Set> myBlankets;
  this->parallelInitialize(myPV, myBlankets);
  // Remember all the local MB ordered pairs
  std::set<std::pair<Var, Var>> myAdded;
  this->growShrink(std::move(myPV), myBlankets, myAdded, imbalanceThreshold);

  /* Symmetry correction */
  this->m_comm.barrier();
  TIMER_START(this->m_tSymmetry);
  auto myPairs = this->symmetryCorrect(std::move(myBlankets), std::move(myAdded));
  this->m_comm.barrier();
  TIMER_PAUSE(this->m_tSymmetry);
  TIMER_PAUSE(this->m_tBlankets);
  /* End of Symmetry Correction */

  TIMER_START(this->m_tNeighbors);
  for (const auto x : this->m_allVars) {
    allBlankets[x] = set_init(Set(), this->m_data.numVars());
  }
  for (const auto& p : myPairs) {
    allBlankets[p.first].insert(p.second);
    allBlankets[p.second].insert(p.first);
  }
  // Sync all the blankets across all the processors
  TIMER_START(this->m_tSync);
  this->syncSets(allBlankets);
  TIMER_PAUSE(this->m_tSync);

  // Get neighbors for the variables on this processor
  for (const auto& p : myPairs) {
    const auto& mbFirst = allBlankets.at(p.first);
    const auto& mbSecond = allBlankets.at(p.second);
    auto mbTest = mbFirst;
    // Pick the smaller of the two MBs
    if (mbSecond.size() > mbFirst.size()) {
      mbTest.erase(p.second);
    }
    else {
      mbTest = mbSecond;
      mbTest.erase(p.first);
    }
    if (!this->m_data.isIndependentAnySubset(p.first, p.second, mbTest, this->m_maxConditioning)) {
      if (allNeighbors.find(p.first) == allNeighbors.end()) {
        allNeighbors[p.first] = set_init(Set(), this->m_data.numVars());
      }
      allNeighbors[p.first].insert(p.second);
    }
  }
  // Now, sync all the neighbors across all the processors
  // XXX: Optimally, we can just gather all the data on process 0
  //      and create and write the network from that process for now.
  //      However, that will not be future proof, since we will need
  //      the network on all the processors when we want to use it.
  auto varNames = this->m_data.varNames(this->m_allVars);
  BayesianNetwork<Var> bn(varNames);
  for (const auto x : this->m_allVars) {
    if (allNeighbors.find(x) == allNeighbors.end()) {
      allNeighbors[x] = set_init(Set(), this->m_data.numVars());
    }
  }
  set_allunion_indexed(allNeighbors, this->m_allVars, this->m_data.numVars(), this->m_comm);
  // We can now create the Bayesian network independently on every processor
  for (const auto x : this->m_allVars) {
    for (const auto y : allNeighbors.at(x)) {
      if (x < y) {
        LOG_MESSAGE_IF(this->m_comm.is_first(), info, "+ Adding the edge %s <-> %s",
                       this->m_data.varName(x), this->m_data.varName(y));
        bn.addEdge(x, y);
        bn.addEdge(y, x);
        // The neighbor set of y has x in it only if x > y
        // Therefore, we need to make it symmetric here before returning
        allNeighbors[y].insert(x);
      }
    }
  }
  this->m_comm.barrier();
  TIMER_PAUSE(this->m_tNeighbors);
  return bn;
}

template <typename Data, typename Var, typename Set>
GS<Data, Var, Set>::GS(
  const mxx::comm& comm,
  const Data& data,
  const Var maxConditioning
) : BlanketLearning<Data, Var, Set>(comm, data, maxConditioning)
{
}

template <typename Data, typename Var, typename Set>
std::pair<Var, double>
GS<Data, Var, Set>::pickBestCandidate(
  const Var target,
  const Set& candidates,
  const Set& cmb
) const
{
  for (const Var y : candidates) {
    LOG_MESSAGE(debug, "Grow: Evaluating %s for addition to the MB", this->m_data.varName(y));
    double pv = this->m_data.pValue(target, y, cmb);
    if (!this->m_data.isIndependent(pv)) {
      LOG_MESSAGE(debug, "Grow: %s chosen as the best candidate", this->m_data.varName(y));
      return std::make_pair(y, pv);
    }
  }
  return std::make_pair(this->m_data.numVars(), 1.0);
}

template <typename Data, typename Var, typename Set>
Set
GS<Data, Var, Set>::getCandidateMB(
  const Var target,
  Set&& candidates
) const
{
  LOG_MESSAGE(info, "%s", std::string(60, '-'));
  LOG_MESSAGE(info, "GS: Getting MB for %s", this->m_data.varName(target));
  auto cmb = set_init(Set(), this->m_data.numVars());
  bool changed = true;
  Var x = this->m_data.numVars();
  double pvX = std::numeric_limits<double>::max();
  TIMER_START(this->m_tGrow);
  while ((candidates.size() > 0) && changed) {
    changed = false;
    std::tie(x, pvX) = this->pickBestCandidate(target, candidates, cmb);
    if (!this->m_data.isIndependent(pvX)) {
      LOG_MESSAGE(info, "+ Adding %s to the MB of %s (p-value = %g)",
                  this->m_data.varName(x), this->m_data.varName(target), pvX);
      cmb.insert(x);
      candidates.erase(x);
      changed = true;
    }
  }
  TIMER_PAUSE(this->m_tGrow);
  TIMER_START(this->m_tShrink);
  this->shrinkMB(target, cmb);
  TIMER_PAUSE(this->m_tShrink);
  LOG_MESSAGE(info, "%s", std::string(60, '-'));
  return cmb;
}

template <typename Data, typename Var, typename Set>
void
GS<Data, Var, Set>::updatePValues(
  std::vector<std::tuple<Var, Var, double>>& myPV,
  const std::unordered_map<Var, Set>& myBlankets
) const
{
  auto primary = std::numeric_limits<Var>::max();
  bool found = false;
  for (auto& mpv : myPV) {
    if (primary != std::get<0>(mpv)) {
      // Updating p-values for a new primary variable; reset
      found = false;
      primary = std::get<0>(mpv);
    }
    if (!found) {
      LOG_MESSAGE(debug, "Updating the p-value for the pair (%s, %s)",
                  this->m_data.varName(primary), this->m_data.varName(std::get<1>(mpv)));
      // A candidate to be added for this primary variable has not been found yet
      // Conduct CI test for this primary-secondary variable pair
      found = !this->m_data.isIndependent(primary, std::get<1>(mpv),
                                          myBlankets.at(std::get<0>(mpv)));
    }
    std::get<2>(mpv) = 1.0 - static_cast<double>(found);
  }
}

template <typename Data, typename Var, typename Set>
IAMB<Data, Var, Set>::IAMB(
  const mxx::comm& comm,
  const Data& data,
  const Var maxConditioning
) : BlanketLearning<Data, Var, Set>(comm, data, maxConditioning)
{
}

template <typename Data, typename Var, typename Set>
Set
IAMB<Data, Var, Set>::getCandidateMB(
  const Var target,
  Set&& candidates
) const
{
  LOG_MESSAGE(info, "%s", std::string(60, '-'));
  LOG_MESSAGE(info, "IAMB: Getting MB for %s", this->m_data.varName(target));
  auto cmb = set_init(Set(), this->m_data.numVars());
  bool changed = true;
  Var x = this->m_data.numVars();
  double pvX = std::numeric_limits<double>::max();
  TIMER_START(this->m_tGrow);
  while ((candidates.size() > 0) && changed) {
    changed = false;
    std::tie(x, pvX) = this->pickBestCandidate(target, candidates, cmb);
    // Add the variable to the candidate MB if it is not
    // independedent of the target
    if (!this->m_data.isIndependent(pvX)) {
      LOG_MESSAGE(info, "+ Adding %s to the MB of %s (p-value = %g)",
                  this->m_data.varName(x), this->m_data.varName(target), pvX);
      cmb.insert(x);
      candidates.erase(x);
      changed = true;
    }
  }
  TIMER_PAUSE(this->m_tGrow);
  TIMER_START(this->m_tShrink);
  this->shrinkMB(target, cmb);
  TIMER_PAUSE(this->m_tShrink);
  LOG_MESSAGE(info, "%s", std::string(60, '-'));
  return cmb;
}

template <typename Data, typename Var, typename Set>
InterIAMB<Data, Var, Set>::InterIAMB(
  const mxx::comm& comm,
  const Data& data,
  const Var maxConditioning
) : BlanketLearning<Data, Var, Set>(comm, data, maxConditioning)
{
}

template <typename Data, typename Var, typename Set>
Set
InterIAMB<Data, Var, Set>::getCandidateMB(
  const Var target,
  Set&& candidates
) const
{
  LOG_MESSAGE(info, "%s", std::string(60, '-'));
  LOG_MESSAGE(info, "InterIAMB: Getting MB for %s", this->m_data.varName(target));
  auto cmb = set_init(Set(), this->m_data.numVars());
  bool changed = true;
  Var x = this->m_data.numVars();
  double pvX = std::numeric_limits<double>::max();
  while ((candidates.size() > 0) && changed) {
    TIMER_START(this->m_tGrow);
    changed = false;
    std::tie(x, pvX) = this->pickBestCandidate(target, candidates, cmb);
    // Add the variable to the candidate MB if it is not
    // independedent of the target
    if (!this->m_data.isIndependent(pvX)) {
      LOG_MESSAGE(info, "+ Adding %s to the MB of %s (p-value = %g)",
                  this->m_data.varName(x), this->m_data.varName(target), pvX);
      cmb.insert(x);
      candidates.erase(x);
      changed = true;
    }
    TIMER_PAUSE(this->m_tGrow);
    if (changed) {
      TIMER_START(this->m_tShrink);
      auto removed = this->shrinkMB(target, cmb);
      TIMER_PAUSE(this->m_tShrink);
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

template <typename Data, typename Var, typename Set>
void
InterIAMB<Data, Var, Set>::growShrink(
  std::vector<std::tuple<Var, Var, double>>&& myPV,
  std::unordered_map<Var, Set>& myBlankets,
  std::set<std::pair<Var, Var>>& myAdded,
  const double imbalanceThreshold
) const
{
  bool changed = true;
  bool redistributed = false;
  while (changed) {
    /* Grow Phase */
    TIMER_START(this->m_tGrow);
    this->updatePValues(myPV, myBlankets);
    auto prevBlankets = myBlankets;
    auto added = this->growAll(myPV, myBlankets);
    TIMER_PAUSE(this->m_tGrow);
    /* End of Grow Phase */
    TIMER_START(this->m_tSync);
    this->syncSets(myBlankets);
    TIMER_PAUSE(this->m_tSync);
    /* Shrink Phase */
    TIMER_START(this->m_tShrink);
    auto removed = this->shrinkAll(myBlankets);
    TIMER_PAUSE(this->m_tShrink);
    /* End of Shrink Phase */
    // Track changes for all the variables across processors
    auto changes = set_init(Set(), this->m_data.numVars());
    std::set<std::pair<Var, Var>> newAdded;
    // Only the candidate blankets which grew can change
    for (const auto& as : added) {
      if (prevBlankets.at(std::get<0>(as)) != myBlankets.at(std::get<0>(as))) {
        // The blanket for this variable changed
        changes.insert(std::get<0>(as));
      }
    }
    set_allunion(changes, this->m_comm);
    if (!changes.empty()) {
      if (!added.empty()) {
        for (const auto& mpv : myPV) {
          if (added.find(mpv) != added.end()) {
            newAdded.insert(std::make_pair(std::get<0>(mpv), std::get<1>(mpv)));
          }
        }
      }
      bool sortPV = false;
      for (const auto& p : removed) {
        if (newAdded.find(p) != newAdded.end()) {
          // If this was a newly added pair during this iteration's grow phase,
          // then we can just remove it from the set
          newAdded.erase(p);
        }
        else if (myAdded.find(p) != myAdded.end()) {
          // Otherwise, if it was a pair belonging to this processor then we
          // need to add it back and sort the list later
          myAdded.erase(p);
          sortPV = true;
          myPV.push_back(std::make_tuple(std::get<0>(p), std::get<1>(p), 0.0));
        }
      }
      // Remove added tuples from future consideration, if they were on this processor
      // Also remove the tuples corresponding to primary variables whose MB did not change
      auto addedOrUnchanged = [&newAdded, &changes] (const std::tuple<Var, Var, double>& mpv)
                                                    { return (newAdded.find(std::make_pair(std::get<0>(mpv),
                                                                                           std::get<1>(mpv))) != newAdded.end()) ||
                                                             !changes.contains(std::get<0>(mpv)); };
      auto newEnd = std::remove_if(myPV.begin(), myPV.end(), addedOrUnchanged);
      myPV.resize(std::distance(myPV.begin(), newEnd));
      if (imbalanceThreshold > 1.0) {
        TIMER_START(this->m_tDist);
        bool fixed = this->fixImbalance(myPV, imbalanceThreshold);
        bool sorted = false;
        if ((redistributed || fixed) && mxx::any_of(sortPV, this->m_comm)) {
          // If the p-value list was redistributed in any of the previous iterations AND
          // there are newly added elements in the list, then we need to globally sort the list
          mxx::comm nonzero_comm(static_cast<MPI_Comm>(this->m_comm));
          if (mxx::any_of(myPV.size() == 0, this->m_comm)) {
            nonzero_comm = this->m_comm.split(myPV.size() > 0);
          }
          if (myPV.size() > 0) {
            mxx::sort(myPV.begin(), myPV.end(), nonzero_comm);
          }
          sorted = true;
        }
        if (fixed || sorted) {
          // We need to get the missing blankets if the p-value list was redistributed in this iteration
          // This can happen if the imbalance was fixed OR if the list was sorted because of an addition
          TIMER_START(this->m_tSync);
          this->syncMissingSets(myPV, myBlankets);
          TIMER_PAUSE(this->m_tSync);
          redistributed = true;
        }
        TIMER_PAUSE(this->m_tDist);
      }
      if (!redistributed && sortPV) {
        // If the p-values have never been redistributed, then we only need to sort locally
        std::sort(myPV.begin(), myPV.end());
      }
      // Finally, include the newly added pairs to the list of this processor's added pairs
      myAdded.insert(newAdded.begin(), newAdded.end());
    }
    else {
      changed = false;
    }
  }
}

#endif // DETAIL_BLANKETLEARNING_HPP_
