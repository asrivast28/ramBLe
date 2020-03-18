/**
 * @file DirectDiscovery.hpp
 * @brief Implementation of the DirectDiscovery as well as
 *        the derived class functions.
 */
#ifndef DETAIL_DIRECTDISCOVERY_HPP_
#define DETAIL_DIRECTDISCOVERY_HPP_

#include "../SetUtils.hpp"

#include "mxx/distribution.hpp"
#include "mxx/reduction.hpp"
#include "mxx/sort.hpp"
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
DirectDiscovery<Data, Var, Set>::getCandidatePC(
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
DirectDiscovery<Data, Var, Set>::updatePValues(
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
DirectDiscovery<Data, Var, Set>::growAll(
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
DirectDiscovery<Data, Var, Set>::shrinkAll(
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
 * @brief Function that synchronizes the candidate blankets by taking a
 *        union of the blankets across all the processors.
 *
 * @param myBlankets A map with the candidate MBs of the primary variables on this processor.
 */
void
DirectDiscovery<Data, Var, Set>::syncBlankets(
  std::unordered_map<Var, Set>& myBlankets
) const
{
  set_allunion_indexed(myBlankets, this->m_allVars, this->m_data.numVars(), this->m_comm);
}

template <typename Data, typename Var, typename Set>
/**
 * @brief Function that fixes the imbalance in the p-value list across processors.
 *
 * @param myPV A list of the p-values corresponding to all the local pairs.
 * @param imbalanceThreshold The amount of imbalance that can be tolerated.
 *
 * @return true if the list was redistributed to fix the imbalance.
 */
bool
DirectDiscovery<Data, Var, Set>::fixImbalance(
  std::vector<std::tuple<Var, Var, double>>& myPV,
  const double imbalanceThreshold
) const
{
  bool fixed = false;
  const auto minSize = mxx::allreduce(myPV.size(), mxx::min<size_t>(), this->m_comm);
  const auto maxSize = mxx::allreduce(myPV.size(), mxx::max<size_t>(), this->m_comm);
  if (minSize * imbalanceThreshold < static_cast<double>(maxSize)) {
    // Redistribute the pairs to fix the imbalance
    mxx::stable_distribute_inplace(myPV, this->m_comm);
    fixed = true;
  }
  return fixed;
}

template <typename Data, typename Var, typename Set>
/**
 * @brief Function that gets the missing candidate blankets for the
 *        primary variables on this processor.
 *
 * @param myPV A list of the p-values corresponding to all the local pairs.
 * @param myBlankets A map with the candidate MBs of the primary variables on this processor.
 */
void
DirectDiscovery<Data, Var, Set>::syncMissingBlankets(
  const std::vector<std::tuple<Var, Var, double>>& myPV,
  std::unordered_map<Var, Set>& myBlankets
) const
{
  for (const auto& mpv : myPV) {
    const auto primary = std::get<0>(mpv);
    if (myBlankets.find(primary) == myBlankets.end()) {
      // This primary variable's blanket was not available on this processor
      // Initialize the blanket for the new primary variable
      myBlankets.insert(std::make_pair(primary, set_init(Set(), this->m_data.numVars())));
    }
  }
  // Synchronize all the blankets
  // XXX: We do not need to synchronize all the blankets. However, some performance testing
  // shows that tracking which blankets should be synced does not result in better performance.
  // Therefore, going with the easier way for now
  this->syncBlankets(myBlankets);
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
DirectDiscovery<Data, Var, Set>::growShrink(
  std::vector<std::tuple<Var, Var, double>>&& myPV,
  std::unordered_map<Var, Set>& myBlankets,
  std::set<std::pair<Var, Var>>& myAdded,
  const double imbalanceThreshold
) const
{
  /* Grow Phase */
  TIMER_DECLARE(tGrow);
  TIMER_DECLARE(tDist);
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
      myPV.resize(std::distance(myPV.begin(), newEnd));
      if (imbalanceThreshold > 1.0) {
        TIMER_START(tDist);
        if (this->fixImbalance(myPV, imbalanceThreshold)) {
          this->syncMissingBlankets(myPV, myBlankets);
        }
        TIMER_PAUSE(tDist);
      }
    }
    else {
      changed = false;
    }
  }
  this->syncBlankets(myBlankets);
  if (this->m_comm.is_first()) {
    TIMER_ELAPSED_NONZERO("Time taken in redistributing: ", tDist);
    TIMER_ELAPSED("Time taken in growing the candidate blankets: ", tGrow);
  }
  /* End of Grow Phase */

  /* Shrink Phase */
  TIMER_DECLARE(tShrink);
  this->shrinkAll(myBlankets);
  if (this->m_comm.is_first()) {
    TIMER_ELAPSED("Time taken in shrinking the candidate blankets: ", tShrink);
  }
  /* End of Shrink Phase */
}

template <typename Data, typename Var, typename Set>
/**
 * @brief Function that performs symmetry correction for all the candidate MBs.
 *
 * @param myBlankets A map with the candidate MBs of the primary variables on this processor.
 * @param myAdded A list of tuples corresponding to the local candidate MB pairs.
 *
 * @return The pairs assigned to this processor after the symmetry correction.
 */
std::vector<std::pair<Var, Var>>
DirectDiscovery<Data, Var, Set>::symmetryCorrect(
  const std::unordered_map<Var, Set>&& myBlankets,
  const std::set<std::pair<Var, Var>>&& myAdded
) const
{
  std::vector<std::pair<Var, Var>> myPairs(myAdded.size());
  auto i = 0u;
  for (const auto& mb : myBlankets) {
    for (const auto secondary : mb.second) {
      // Create an ordered pair corresponding to this pair only
      // if it was a pair which was added on this processor
      if (myAdded.find(std::make_pair(mb.first, secondary)) != myAdded.end()) {
        if (mb.first < secondary) {
          myPairs[i] = std::make_pair(mb.first, secondary);
        }
        else {
          myPairs[i] = std::make_pair(secondary, mb.first);
        }
        ++i;
      }
    }
  }
  myPairs.resize(i);
  // Redistribute and sort the list across all the processors
  mxx::stable_distribute_inplace(myPairs, this->m_comm);
  mxx::sort(myPairs.begin(), myPairs.end(), this->m_comm);
  // There should be a duplicate for every ordered pair,
  // otherwise symmetry correction is required
  std::vector<bool> remove(myPairs.size(), true);
  // First check local pairs
  auto it = myPairs.begin();
  while (it != myPairs.end()) {
    it = std::adjacent_find(it, myPairs.end());
    if (it != myPairs.end()) {
      auto d = std::distance(myPairs.begin(), it);
      remove[d] = false;
      remove[d+1] = false;
      it += 2;
    }
  }
  // Now, check the boundary pairs
  auto elem = std::make_pair(static_cast<Var>(0), static_cast<Var>(0));
  // First, check if the first pair on this processor matches the
  // last pair on the previous processor
  if (myPairs.size() > 0) {
    elem = *myPairs.rbegin();
  }
  auto left = mxx::right_shift(elem, this->m_comm);
  if (*remove.begin() && (left == *myPairs.begin())) {
    *remove.begin() = false;
  }
  // Then, check if the last pair on this processor matches the
  // first pair on the next processor
  if (myPairs.size() > 0) {
    elem = *myPairs.begin();
  }
  auto right = mxx::left_shift(elem, this->m_comm);
  if (*remove.rbegin() && (right == *myPairs.rbegin())) {
    *remove.rbegin() = false;
  }
  // Remove the pairs which did not have duplicates
  it = myPairs.begin();
  for (const auto r : remove) {
    if (r) {
      LOG_MESSAGE(info, "- Removing %s from the MB of %s (asymmetry)",
                        this->m_data.varName(it->second), this->m_data.varName(it->first));
      it = myPairs.erase(it);
    }
    else {
      ++it;
    }
  }
  // Only retain unique elements now
  it = std::unique(myPairs.begin(), it);
  myPairs.resize(std::distance(myPairs.begin(), it));
  mxx::stable_distribute_inplace(myPairs, this->m_comm);

  return myPairs;
}

template <typename Data, typename Var, typename Set>
/**
 * @brief Function for getting the undirected skeleton network in parallel.
 *
 * @param imbalanceThreshold Specifies the amount of imbalance the algorithm should tolerate.
 */
BayesianNetwork<Var>
DirectDiscovery<Data, Var, Set>::getSkeleton_parallel(
  const double imbalanceThreshold
) const
{
  TIMER_DECLARE(tBlankets);
  // First, block decompose all the variable pairs on all the processors
  auto n = this->m_allVars.size();
  auto totalPairs = n * (n - 1);
  auto batchSize = (totalPairs / this->m_comm.size()) + ((totalPairs % this->m_comm.size() != 0) ? 1 : 0);
  auto myOffset = std::min(this->m_comm.rank() * batchSize, totalPairs);
  auto mySize = std::min(batchSize, totalPairs - myOffset);

  auto vars = std::vector<Var>(this->m_allVars.begin(), this->m_allVars.end());
  auto primary = vars.begin() + (myOffset / (n - 1));
  auto secondary = vars.begin() + (myOffset % (n - 1));
  if (primary == secondary) {
    ++secondary;
  }

  std::vector<std::tuple<Var, Var, double>> myPV(mySize);
  std::unordered_map<Var, Set> myBlankets;
  myBlankets.insert(std::make_pair(*primary, set_init(Set(), this->m_data.numVars())));
  for (auto i = 0u; i < mySize; ++i, ++secondary) {
    if (secondary == vars.end()) {
      // Start a new primary variable if the secondary variable is exhausted
      secondary = vars.begin();
      ++primary;
      myBlankets.insert(std::make_pair(*primary, set_init(Set(), this->m_data.numVars())));
    }
    if (secondary == primary) {
      // Increment secondary variable once more if it is the same as the primary variable
      ++secondary;
    }
    // Initialize the p-values
    myPV[i] = std::make_tuple(*primary, *secondary, std::numeric_limits<double>::max());
  }

  // Remember all the local MB ordered pairs
  std::set<std::pair<Var, Var>> myAdded;
  this->growShrink(std::move(myPV), myBlankets, myAdded, imbalanceThreshold);

  /* Symmetry correction */
  this->m_comm.barrier();
  TIMER_DECLARE(tSymmetry);
  auto myPairs = this->symmetryCorrect(std::move(myBlankets), std::move(myAdded));
  this->m_comm.barrier();
  if (this->m_comm.is_first()) {
    TIMER_ELAPSED("Time taken in symmetry correcting the blankets: ", tSymmetry);
    TIMER_ELAPSED("Time taken in getting the blankets: ", tBlankets);
  }
  /* End of Symmetry Correction */

  TIMER_DECLARE(tNeighbors);
  std::unordered_map<Var, Set> allBlankets;
  for (const auto x : this->m_allVars) {
    allBlankets[x] = set_init(Set(), this->m_data.numVars());
  }
  for (const auto& p : myPairs) {
    allBlankets[p.first].insert(p.second);
    allBlankets[p.second].insert(p.first);
  }
  // Sync all the blankets across all the processors
  this->syncBlankets(allBlankets);

  // Get neighbors for the variables on this processor
  std::unordered_map<Var, Set> myNeighbors;
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
      if (myNeighbors.find(p.first) == myNeighbors.end()) {
        myNeighbors[p.first] = set_init(Set(), this->m_data.numVars());
      }
      myNeighbors[p.first].insert(p.second);
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
    if (myNeighbors.find(x) == myNeighbors.end()) {
      myNeighbors.insert(std::make_pair(x, set_init(Set(), this->m_data.numVars())));
    }
  }
  set_allunion_indexed(myNeighbors, this->m_allVars, this->m_data.numVars(), this->m_comm);
  // We can now create the Bayesian network independently on every processor
  for (const auto x : this->m_allVars) {
    for (const auto y : myNeighbors.at(x)) {
      LOG_MESSAGE_IF(this->m_comm.is_first(), info, "+ Adding the edge %s <-> %s",
                     this->m_data.varName(x), this->m_data.varName(y));
      bn.addEdge(x, y);
      bn.addEdge(y, x);
    }
  }
  this->m_comm.barrier();
  if (this->m_comm.is_first()) {
    TIMER_ELAPSED("Time taken in getting the neighbors: ", tNeighbors);
  }
  return bn;
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
GSMB<Data, Var, Set>::getCandidateMB(
  const Var target,
  Set&& candidates
) const
{
  LOG_MESSAGE(info, "%s", std::string(60, '-'));
  LOG_MESSAGE(info, "GSMB: Getting MB for %s", this->m_data.varName(target));
  auto cmb = set_init(Set(), this->m_data.numVars());
  bool changed = true;
  Var x = this->m_data.numVars();
  double pvX = std::numeric_limits<double>::max();
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
  this->shrinkMB(target, cmb);
  LOG_MESSAGE(info, "%s", std::string(60, '-'));
  return cmb;
}

template <typename Data, typename Var, typename Set>
void
GSMB<Data, Var, Set>::updatePValues(
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
) : DirectDiscovery<Data, Var, Set>(comm, data, maxConditioning)
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

template <typename Data, typename Var, typename Set>
void
InterIAMB<Data, Var, Set>::growShrink(
  std::vector<std::tuple<Var, Var, double>>&& myPV,
  std::unordered_map<Var, Set>& myBlankets,
  std::set<std::pair<Var, Var>>& myAdded,
  const double imbalanceThreshold
) const
{
  TIMER_DECLARE(tGrow);
  TIMER_DECLARE(tShrink);
  TIMER_DECLARE(tDist);
  TIMER_RESET(tDist);
  TIMER_DECLARE(tSync);
  bool changed = true;
  bool redistributed = false;
  while (changed) {
    /* Grow Phase */
    TIMER_START(tGrow);
    this->updatePValues(myPV, myBlankets);
    auto prevBlankets = myBlankets;
    auto added = this->growAll(myPV, myBlankets);
    TIMER_PAUSE(tGrow);
    /* End of Grow Phase */
    TIMER_START(tSync);
    this->syncBlankets(myBlankets);
    TIMER_PAUSE(tSync);
    /* Shrink Phase */
    TIMER_START(tShrink);
    auto removed = this->shrinkAll(myBlankets);
    TIMER_PAUSE(tShrink);
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
        TIMER_START(tDist);
        bool fixed = this->fixImbalance(myPV, imbalanceThreshold);
        bool sorted = false;
        if ((redistributed || fixed) && mxx::any_of(sortPV, this->m_comm)) {
          // If the p-value list was redistributed in any of the previous iterations AND
          // there are newly added elements in the list, then we need to globally sort the list
          mxx::sort(myPV.begin(), myPV.end(), this->m_comm);
          sorted = true;
        }
        if (fixed || sorted) {
          // We need to get the missing blankets if the p-value list was redistributed in this iteration
          // This can happen if the imbalance was fixed OR if the list was sorted because of an addition
          TIMER_START(tSync);
          this->syncMissingBlankets(myPV, myBlankets);
          TIMER_PAUSE(tSync);
          redistributed = true;
        }
        TIMER_PAUSE(tDist);
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
  if (this->m_comm.is_first()) {
    TIMER_ELAPSED_NONZERO("Time taken in redistributing: ", tDist);
    TIMER_ELAPSED("Time taken in synchronizing the candidate blankets: ", tSync);
    TIMER_ELAPSED("Time taken in growing the candidate blankets: ", tGrow);
    TIMER_ELAPSED("Time taken in shrinking the candidate blankets: ", tShrink);
  }
}

#endif // DETAIL_DIRECTDISCOVERY_HPP_
