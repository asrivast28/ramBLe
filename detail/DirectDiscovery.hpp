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

  return !this->m_data.isIndependentAnySubset(x, y, mbTest, this->m_maxConditioning, this->m_comm);
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
  LOG_MESSAGE(info, "Neighbors: Getting PC from MB for %s",
              this->m_data.varName(target));
  auto cpc = set_init(Set(), this->m_data.numVars());
  auto mb = this->getMB(target);
  for (const Var y: mb) {
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
 * @brief Function that updates scores for all the local pairs, given the current candidate MBs.
 *
 * @param myScores A list of the scores corresponding to all the local pairs.
 * @param myBlankets A map with the candidate MBs of the primary variables on this processor.
 */
void
DirectDiscovery<Data, Var, Set>::updateScores(
  std::vector<std::tuple<Var, Var, double>>& myScores,
  const std::unordered_map<Var, Set>& myBlankets
) const
{
  for (auto& score : myScores) {
    LOG_MESSAGE(debug, "Updating the score for the pair (%s, %s)",
                this->m_data.varName(std::get<0>(score)), this->m_data.varName(std::get<1>(score)));
    std::get<2>(score) = this->m_data.assocScore(std::get<0>(score),
                                                 std::get<1>(score),
                                                 myBlankets.at(std::get<0>(score)));
  }
}

template <typename Data, typename Var, typename Set>
/**
 * @brief Function that grows candidate MBs of the primary variables on this processor.
 *
 * @param myScores A list of the scores corresponding to all the local pairs.
 * @param myBlankets A map with the candidate MBs of the primary variables on this processor.
 *
 * @return Set of all the scores which were added to the candidate MBs on this processor.
 */
std::set<std::tuple<Var, Var, double>>
DirectDiscovery<Data, Var, Set>::growAll(
  const std::vector<std::tuple<Var, Var, double>>& myScores,
  std::unordered_map<Var, Set>& myBlankets
) const
{
  std::set<std::tuple<Var, Var, double>> added;
  std::vector<std::tuple<Var, Var, double>> maxScores(myScores.size());
  auto compareScores = [] (const std::tuple<Var, Var, double>& a, const std::tuple<Var, Var, double>& b)
                          { return (std::get<0>(a) == std::get<0>(b)) ?
                                   (std::isgreaterequal(std::get<2>(a), std::get<2>(b)) ? a : b) : b; };
  // First, do a forward segmented parallel prefix with primary variable defining the segment boundaries
  // This will get the secondary variable with the max score to the corresponding primary variable boundary
  mxx::global_scan(myScores.begin(), myScores.end(), maxScores.begin(), compareScores, this->m_comm, false);
  // Then, do a reverse segmented parallel prefix with the same segments as before
  // This will effectively broadcast the secondary variable with the maximum score within the segments
  mxx::global_scan_inplace(maxScores.rbegin(), maxScores.rend(), compareScores, this->m_comm.reverse(), false);
  // There might be multiple local copies of the max score corresponding to every segment
  // Retain only one per segment
  auto comparePrimary = [] (const std::tuple<Var, Var, double>& a, const std::tuple<Var, Var, double>& b)
                           { return std::get<0>(a) == std::get<0>(b); };
  auto uniqueEnd = std::unique(maxScores.begin(), maxScores.end(), comparePrimary);
  auto changes = set_init(Set(), this->m_data.numVars());
  for (auto it = maxScores.begin(); it != uniqueEnd; ++it) {
    if (!this->m_data.isIndependent(std::get<2>(*it))) {
      // Need to keep looping as long as even one variable's MB keeps changing
      changes.insert(std::get<0>(*it));
      // Add y to the blanket of x
      LOG_MESSAGE(info, "%d: + Adding %s to the MB of %s (score = %g)", this->m_comm.rank(),
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
 * @brief Function that performs grow-shrink for multiple candidate MBs.
 *
 * @param myScores A list of the scores corresponding to all the local pairs.
 * @param myBlankets A map with all the local candidate MBs.
 * @param myAdded A list of all the local candidate MB pairs.
 *
 * @return The total size of all the blankets at the end.
 */
void
DirectDiscovery<Data, Var, Set>::growShrink(
  std::vector<std::tuple<Var, Var, double>>& myScores,
  std::unordered_map<Var, Set>& myBlankets,
  std::set<std::pair<Var, Var>>& myAdded
) const
{
  /* Grow Phase */
  TIMER_DECLARE(tGrow);
  bool changed = true;
  while (changed) {
    this->updateScores(myScores, myBlankets);
    auto added = this->growAll(myScores, myBlankets);
    // Track blanket changes for all the primary variables across processors
    auto changes = set_init(Set(), this->m_data.numVars());
    for (const auto& as : added) {
      changes.insert(std::get<0>(as));
    }
    set_allunion(changes, this->m_comm);
    if (!changes.empty()) {
      if (!added.empty()) {
        // Record the added scores belonging to this processor
        for (const auto& ms : myScores) {
          if (added.find(ms) != added.end()) {
            myAdded.insert(std::make_pair(std::get<0>(ms), std::get<1>(ms)));
          }
        }
      }
      // Remove added scores from future consideration, if they were on this processor
      // Also remove the scores corresponding to primary variables whose MB did not change
      auto addedOrUnchanged = [&added, &changes] (const std::tuple<Var, Var, double>& ms)
                                                 { return (added.find(ms) != added.end()) ||
                                                          !changes.contains(std::get<0>(ms)); };
      auto newEnd = std::remove_if(myScores.begin(), myScores.end(), addedOrUnchanged);
      myScores.resize(std::distance(myScores.begin(), newEnd));
    }
    else {
      changed = false;
    }
  }
  if (this->m_comm.is_first()) {
    TIMER_ELAPSED("Time taken in growing the candidate blankets: ", tGrow);
  }
  /* End of Grow Phase */

  /* Shrink Phase */
  TIMER_DECLARE(tShrink);
  auto removed = this->shrinkAll(myBlankets);
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
 * @param myAdded A list of scores corresponding to the local candidate MB pairs.
 *
 * @return The pairs assigned to this processor after the symmetry correction.
 */
std::vector<std::pair<Var, Var>>
DirectDiscovery<Data, Var, Set>::symmetryCorrect(
  const std::unordered_map<Var, Set>& myBlankets,
  const std::set<std::pair<Var, Var>>& myAdded
) const
{
  std::vector<std::pair<Var, Var>> myPairs(myAdded.size());
  auto i = 0u;
  for (const auto& mb : myBlankets) {
    for (const auto secondary : mb.second) {
      // Create an ordered pair corresponding to this pair only if
      // it was a pair originally assigned to this processor
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
 */
BayesianNetwork<Var>
DirectDiscovery<Data, Var, Set>::getSkeleton_parallel(
) const
{
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

  std::vector<std::tuple<Var, Var, double>> myScores(mySize);
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
    // Compute the initial marginal associations
    myScores[i] = std::make_tuple(*primary, *secondary, 0.0);
  }

  // Remember all the local MB ordered pairs
  std::set<std::pair<Var, Var>> myAdded;
  this->growShrink(myScores, myBlankets, myAdded);

  /* Symmetry correction */
  TIMER_DECLARE(tSymmetry);
  auto myPairs = this->symmetryCorrect(myBlankets, myAdded);
  if (this->m_comm.is_first()) {
    TIMER_ELAPSED("Time taken in symmetry correcting the blankets: ", tSymmetry);
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
  for (const auto x : this->m_allVars) {
    set_allunion(allBlankets[x], this->m_comm);
  }

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
  // and create the bayesian network
  auto varNames = this->m_data.varNames(this->m_allVars);
  BayesianNetwork<Var> bn(varNames);
  for (const auto x : this->m_allVars) {
    auto pcX = set_init(Set(), this->m_data.numVars());
    if (myNeighbors.find(x) != myNeighbors.end()) {
      pcX = myNeighbors.at(x);
    }
    set_allunion(pcX, this->m_comm);
    for (const auto y : pcX) {
      LOG_MESSAGE_IF(this->m_comm.is_first(), info, "+ Adding the edge %s <-> %s",
                     this->m_data.varName(x), this->m_data.varName(y));
      bn.addEdge(x, y);
      bn.addEdge(y, x);
    }
  }
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
      LOG_MESSAGE(info, "+ Adding %s to the MB of %s (score = %g)",
                  this->m_data.varName(x), this->m_data.varName(target), scoreX);
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
GSMB<Data, Var, Set>::updateScores(
  std::vector<std::tuple<Var, Var, double>>& myScores,
  const std::unordered_map<Var, Set>& myBlankets
) const
{
  auto score = myScores.begin();
  auto primary = std::numeric_limits<Var>::max();
  bool found = false;
  while (score != myScores.end()) {
    if (primary != std::get<0>(*score)) {
      // Updating scores for a new primary variable; reset
      found = false;
      primary = std::get<0>(*score);
    }
    if (!found) {
      LOG_MESSAGE(debug, "Updating the score for the pair (%s, %s)",
                  this->m_data.varName(primary), this->m_data.varName(std::get<1>(*score)));
      // A candidate to be added for this primary variable has not been found yet
      // Conduct CI test for this primary-secondary variable pair
      found = !this->m_data.isIndependent(primary, std::get<1>(*score),
                                          myBlankets.at(std::get<0>(*score)));
    }
    std::get<2>(*score) = static_cast<double>(found);
    ++score;
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
      LOG_MESSAGE(info, "+ Adding %s to the MB of %s",
                  this->m_data.varName(x), this->m_data.varName(target));
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

template <typename Data, typename Var, typename Set>
void
InterIAMB<Data, Var, Set>::growShrink(
  std::vector<std::tuple<Var, Var, double>>& myScores,
  std::unordered_map<Var, Set>& myBlankets,
  std::set<std::pair<Var, Var>>& myAdded
) const
{
  TIMER_DECLARE(tGrow);
  TIMER_DECLARE(tShrink);
  bool changed = true;
  while (changed) {
    /* Grow Phase */
    TIMER_START(tGrow);
    this->updateScores(myScores, myBlankets);
    auto prevBlankets = myBlankets;
    auto added = this->growAll(myScores, myBlankets);
    TIMER_PAUSE(tGrow);
    /* End of Grow Phase */
    /* Shrink Phase */
    TIMER_START(tShrink);
    auto removed = this->shrinkAll(myBlankets);
    TIMER_PAUSE(tShrink);
    /* End of Shrink Phase */
    // Track changes for all the variables across processors
    auto changes = set_init(Set(), this->m_data.numVars());
    std::set<std::pair<Var, Var>> newAdded;
    for (const auto& as : added) {
      if (prevBlankets.at(std::get<0>(as)) != myBlankets.at(std::get<0>(as))) {
        // The blanket for this variable changed
        changes.insert(std::get<0>(as));
      }
    }
    set_allunion(changes, this->m_comm);
    if (!changes.empty()) {
      if (!added.empty()) {
        for (const auto& ms : myScores) {
          if (added.find(ms) != added.end()) {
            newAdded.insert(std::make_pair(std::get<0>(ms), std::get<1>(ms)));
          }
        }
      }
      bool sortScores = false;
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
          sortScores = true;
          myScores.push_back(std::make_tuple(std::get<0>(p), std::get<1>(p), 0.0));
        }
      }
      // Remove added scores from future consideration, if they were on this processor
      // Also remove the scores corresponding to primary variables whose MB did not change
      auto addedOrUnchanged = [&newAdded, &changes] (const std::tuple<Var, Var, double>& ms)
                                                    { return (newAdded.find(std::make_pair(std::get<0>(ms),
                                                                                           std::get<1>(ms))) != newAdded.end()) ||
                                                             !changes.contains(std::get<0>(ms)); };
      auto newEnd = std::remove_if(myScores.begin(), myScores.end(), addedOrUnchanged);
      myScores.resize(std::distance(myScores.begin(), newEnd));
      if (sortScores) {
        std::sort(myScores.begin(), myScores.end());
      }
      // Finally, include the newly added pairs to the list of this processor's added pairs
      myAdded.insert(newAdded.begin(), newAdded.end());
    }
    else {
      changed = false;
    }
  }
  if (this->m_comm.is_first()) {
    TIMER_ELAPSED("Time taken in growing the candidate blankets: ", tGrow);
    TIMER_ELAPSED("Time taken in shrinking the candidate blankets: ", tShrink);
  }
}

#endif // DETAIL_DIRECTDISCOVERY_HPP_
