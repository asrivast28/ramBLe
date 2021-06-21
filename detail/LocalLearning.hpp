/**
 * @file LocalLearning.hpp
 * @brief Implementation of the classes for local-to-global causal discovery.
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
#ifndef DETAIL_LOCALLEARNING_HPP_
#define DETAIL_LOCALLEARNING_HPP_

#include "mxx/sort.hpp"


template <typename Data, typename Var, typename Set>
/**
 * @brief Constructs the object with the given data.
 *
 * @param data Reference to an object of the Data.
 */
LocalLearning<Data, Var, Set>::LocalLearning(
  const mxx::comm& comm,
  const Data& data,
  const double alpha,
  const Var maxConditioning
) : ConstraintBasedLearning<Data, Var, Set>(comm, data, alpha, maxConditioning),
    m_cachedPC(),
    m_cachedMB(),
    m_cachedPCSymmetric(),
    m_cachedMBSymmetric()
{
}

template <typename Data, typename Var, typename Set>
/**
 * @brief Finds the candidate PC for a variable, and caches the result.
 *
 * @param target The index of the target variable.
 * @param candidates The indices of all the candidate variables.
 *
 * @return A pair with the first element being a set containing the indices of
 *         all the variables in the candidate PC of the given target variable,
 *         and the second element specifying if the set has been symmetry corrected.
 */
Set&
LocalLearning<Data, Var, Set>::getCandidatePC_cache(
  const Var target,
  Set&& candidates
) const
{
  auto cacheIt = m_cachedPC.find(target);
  LOG_MESSAGE_IF(cacheIt != m_cachedPC.end(), trace, "* Found candidate PC for %s in the cache",
                                                     this->m_data.varName(target));
  if (cacheIt == m_cachedPC.end()) {
    auto cpc = this->getCandidatePC(target, std::move(candidates));
    cacheIt = m_cachedPC.insert(cacheIt, std::make_pair(target, cpc));
    m_cachedPCSymmetric[target] = false;
  }
  return cacheIt->second;
}

template <typename Data, typename Var, typename Set>
/**
 * @brief Performs symmetry correction on the given PC set for the given variable.
 *
 * @param target The index of the target variable.
 * @param cpc The set containing the indices of all the variables in
 *            the candidate PC set.
 */
void
LocalLearning<Data, Var, Set>::symmetryCorrectPC(
  const Var target,
  Set& cpc
) const
{
  auto initial = cpc;
  for (const Var x : initial) {
    auto candidatesX = this->getCandidates(x);
    const auto& cpcX = this->getCandidatePC_cache(x, std::move(candidatesX));
    if (!set_contains(cpcX, target)) {
      LOG_MESSAGE(info, "- Removing %s from the PC of %s (asymmetry)", this->m_data.varName(x), this->m_data.varName(target));
      cpc.erase(x);
    }
  }
}

template <typename Data, typename Var, typename Set>
/**
 * @brief The top level function for getting the correct PC set
 *        for the given target variable.
 *
 * @param target The index of the target variable.
 *
 * @return A set containing the indices of all the variables
 *         in the correct PC set of the target variable.
 */
const Set&
LocalLearning<Data, Var, Set>::getPC(
  const Var target
) const
{
  auto candidates = this->getCandidates(target);
  auto& cpc = this->getCandidatePC_cache(target, std::move(candidates));
  if (!m_cachedPCSymmetric.at(target)) {
    this->symmetryCorrectPC(target, cpc);
    m_cachedPCSymmetric[target] = true;
  }
  return cpc;
}

template <typename Data, typename Var, typename Set>
/**
 * @brief Finds the candidate MB for a variable, and caches the result.
 *
 * @param target The index of the target variable.
 * @param candidates The indices of all the candidate variables.
 *
 * @return A pair with the first element being a set containing the indices of
 *         all the variables in the candidate MB of the given target variable,
 *         and the second element specifying if the set has been symmetry corrected.
 */
Set&
LocalLearning<Data, Var, Set>::getCandidateMB_cache(
  const Var target,
  Set&& candidates
) const
{
  auto cacheIt = m_cachedMB.find(target);
  LOG_MESSAGE_IF(cacheIt != m_cachedMB.end(), trace, "* Found candidate MB for %s in the cache",
                                                     this->m_data.varName(target));
  if (cacheIt == m_cachedMB.end()) {
    auto cmb = this->getCandidateMB(target, std::move(candidates));
    cacheIt = m_cachedMB.insert(cacheIt, std::make_pair(target, cmb));
    m_cachedMBSymmetric[target] = false;
  }
  return cacheIt->second;
}

template <typename Data, typename Var, typename Set>
/**
 * @brief Symmetry corrects the candidate MB of the target variable.
 *
 * @param target The index of the target variable.
 * @param cmb The indices of the variables in the candidate MB of the target variable.
 *                The function removes the indices from the candidate set.
 */
void
LocalLearning<Data, Var, Set>::symmetryCorrectMB(
  const Var target,
  Set& cmb
) const
{
  auto initial = cmb;
  for (const Var x : initial) {
    auto candidatesX = this->getCandidates(x);
    const auto& cmbX = this->getCandidateMB_cache(x, std::move(candidatesX));
    if (!set_contains(cmbX, target)) {
      LOG_MESSAGE(info, "- Removing %s from the MB of %s (asymmetry)", this->m_data.varName(x), this->m_data.varName(target));
      cmb.erase(x);
    }
  }
}

template <typename Data, typename Var, typename Set>
/**
 * @brief Top level function for getting the MB of the target variable.
 *
 * @param target The index of the target variable.
 */
const Set&
LocalLearning<Data, Var, Set>::getMB(
  const Var target
) const
{
  auto candidates = this->getCandidates(target);
  auto& cmb = this->getCandidateMB_cache(target, std::move(candidates));
  if (!m_cachedMBSymmetric.at(target)) {
    this->symmetryCorrectMB(target, cmb);
    m_cachedMBSymmetric[target] = true;
  }
  return cmb;
}

template <typename Data, typename Var, typename Set>
/**
 * @brief Function for getting the undirected skeleton network sequentially.
 */
BayesianNetwork<Var>
LocalLearning<Data, Var, Set>::getSkeleton_sequential(
  const bool
) const
{
  BayesianNetwork<Var> bn(this->m_data.varNames(this->m_allVars));
  for (const auto x : this->m_allVars) {
    const auto& pcX = this->getPC(x);
    for (const auto y : pcX) {
      if (x < y) {
        LOG_MESSAGE(info, "+ Adding the edge %s <-> %s", this->m_data.varName(x), this->m_data.varName(y));
        bn.addEdge(x, y, true);
      }
    }
  }
  return bn;
}

template <typename Data, typename Var, typename Set>
/**
 * @brief Function for initializing the data structures used for
 *        learning the skeleton in parallel.
 *
 * @param myPV The array to be initialized with the variable pairs
 *             to be handled by this processor.
 * @param mySets The map containing the MB or PC sets
 *               for all the primary variables on this processor.
 */
void
LocalLearning<Data, Var, Set>::parallelInitialize(
  std::vector<std::tuple<Var, Var, double>>& myPV,
  std::unordered_map<Var, Set>& mySets
) const
{
  // First, block decompose all the variable pairs on all the processors
  auto n = this->m_allVars.size();
  mxx::blk_dist dist(n * (n - 1), this->m_comm);
  auto myOffset = dist.eprefix_size();
  auto mySize = dist.local_size();

  auto vars = std::vector<Var>(this->m_allVars.begin(), this->m_allVars.end());
  auto primary = vars.begin() + (myOffset / (n - 1));
  auto secondary = vars.begin() + (myOffset % (n - 1));
  if (secondary >= primary) {
    ++secondary;
  }

  myPV.resize(mySize);
  mySets.insert(std::make_pair(*primary, set_init(Set(), this->m_data.numVars())));
  for (auto i = 0u; i < mySize; ++i, ++secondary) {
    if (secondary == vars.end()) {
      // Start a new primary variable if the secondary variable is exhausted
      secondary = vars.begin();
      ++primary;
      mySets.insert(std::make_pair(*primary, set_init(Set(), this->m_data.numVars())));
    }
    if (secondary == primary) {
      // Increment secondary variable once more if it is the same as the primary variable
      ++secondary;
    }
    // Initialize the p-values
    myPV[i] = std::make_tuple(*primary, *secondary, 0.0);
  }
}

template <typename Data, typename Var, typename Set>
/**
 * @brief Function that gets the missing candidate sets for the
 *        primary variables on this processor.
 *
 * @param myPV A list of the p-values corresponding to all the local pairs.
 * @param mySets A map with the candidate MB or PC sets of the
 *               primary variables on this processor.
 */
void
LocalLearning<Data, Var, Set>::syncMissingSets(
  const std::vector<std::tuple<Var, Var, double>>& myPV,
  std::unordered_map<Var, Set>& mySets
) const
{
  for (const auto& mpv : myPV) {
    const auto primary = std::get<0>(mpv);
    if (mySets.find(primary) == mySets.end()) {
      // This primary variable's candidate set was not available on this processor
      // Initialize the set for the new primary variable
      mySets.insert(std::make_pair(primary, set_init(Set(), this->m_data.numVars())));
    }
  }
  // Synchronize all the sets
  // XXX: We do not need to synchronize all the sets. However, some performance testing
  // shows that tracking which sets should be synced does not result in better performance.
  // Therefore, going with the easier way for now
  this->syncSets(mySets);
}

template <typename Data, typename Var, typename Set>
/**
 * @brief Function that performs symmetry correction for all the candidate MBs or PCs.
 *
 * @param mySets The map containing the MB or PC sets
 *               for all the primary variables on this processor.
 * @param myAdded A list of tuples corresponding to the local candidate MB pairs.
 *
 * @return The pairs assigned to this processor after the symmetry correction.
 */
std::vector<std::pair<Var, Var>>
LocalLearning<Data, Var, Set>::symmetryCorrect(
  const std::unordered_map<Var, Set>&& mySets,
  const std::set<std::pair<Var, Var>>&& myAdded
) const
{
  std::vector<std::pair<Var, Var>> myPairs(myAdded.size());
  auto i = 0u;
  for (const auto& ms : mySets) {
    for (const auto secondary : ms.second) {
      // Create an ordered pair corresponding to this pair only
      // if it was a pair which was added on this processor
      if (myAdded.find(std::make_pair(ms.first, secondary)) != myAdded.end()) {
        if (ms.first < secondary) {
          myPairs[i] = std::make_pair(ms.first, secondary);
        }
        else {
          myPairs[i] = std::make_pair(secondary, ms.first);
        }
        ++i;
      }
    }
  }
  myPairs.resize(i);
  // Redistribute and sort the list across all the processors
  TIMER_START(this->m_tMxx);
  mxx::stable_distribute_inplace(myPairs, this->m_comm);
  mxx::comm nonzero_comm(static_cast<MPI_Comm>(this->m_comm));
  if (mxx::any_of(myPairs.size() == 0, this->m_comm)) {
    nonzero_comm = this->m_comm.split(myPairs.size() > 0);
  }
  TIMER_PAUSE(this->m_tMxx);
  if (myPairs.size() > 0) {
    TIMER_START(this->m_tMxx);
    mxx::sort(myPairs.begin(), myPairs.end(), nonzero_comm);
    TIMER_PAUSE(this->m_tMxx);
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
    // First, check if the first pair on this processor matches the
    // last pair on the previous processor
    auto elem = *myPairs.rbegin();
    TIMER_START(this->m_tMxx);
    auto left = mxx::right_shift(elem, nonzero_comm);
    TIMER_PAUSE(this->m_tMxx);
    if (*remove.begin() && (left == *myPairs.begin())) {
      *remove.begin() = false;
    }
    // Then, check if the last pair on this processor matches the
    // first pair on the next processor
    elem = *myPairs.begin();
    TIMER_START(this->m_tMxx);
    auto right = mxx::left_shift(elem, nonzero_comm);
    TIMER_PAUSE(this->m_tMxx);
    if (*remove.rbegin() && (right == *myPairs.rbegin())) {
      *remove.rbegin() = false;
    }
    // Remove the pairs which did not have duplicates
    it = myPairs.begin();
    for (const auto r : remove) {
      if (r) {
        LOG_MESSAGE(info, "- Removing %s from the candidate set of %s (asymmetry)",
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
  }
  TIMER_START(this->m_tMxx);
  mxx::stable_distribute_inplace(myPairs, this->m_comm);
  TIMER_PAUSE(this->m_tMxx);

  return myPairs;
}

#endif // DETAIL_LOCALLEARNING_HPP_
