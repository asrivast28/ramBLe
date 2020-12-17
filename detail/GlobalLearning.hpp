/**
 * @file GlobalLearning.hpp
 * @brief Implementation of the classes for global learning algorithms.
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
#ifndef DETAIL_GLOBALLEARNING_HPP_
#define DETAIL_GLOBALLEARNING_HPP_

#include "SetUtils.hpp"

#include "utils/Logging.hpp"


template <typename Data, typename Var, typename Set>
/**
 * @brief Constructs the object with the given data.
 */
GlobalLearning<Data, Var, Set>::GlobalLearning(
  const mxx::comm& comm,
  const Data& data,
  const double alpha,
  const Var maxConditioning
) : ConstraintBasedLearning<Data, Var, Set>(comm, data, alpha, maxConditioning),
    m_cachedNeighbors(),
    m_removedEdges()
{
  TIMER_RESET(m_tSync);
  TIMER_RESET(m_tDist);
  TIMER_RESET(m_tRemoved);
}

template <typename Data, typename Var, typename Set>
/**
 * @brief Default destructor that prints out timing information.
 */
GlobalLearning<Data, Var, Set>::~GlobalLearning(
)
{
  if (this->m_comm.is_first()) {
    TIMER_ELAPSED_NONZERO("Time taken in redistributing: ", m_tDist);
    TIMER_ELAPSED_NONZERO("Time taken in synchronizing neighbors: ", m_tSync);
    TIMER_ELAPSED_NONZERO("Time taken in tracking removed edges: ", m_tRemoved);
  }
}

template <typename Data, typename Var, typename Set>
const Set&
GlobalLearning<Data, Var, Set>::getPC(
  const Var target
) const
{
  return m_cachedNeighbors.at(target);
}

template <typename Data, typename Var, typename Set>
const Set&
GlobalLearning<Data, Var, Set>::getMB(
  const Var target
) const
{
  return m_cachedNeighbors.at(target);
}

template <typename Data, typename Var, typename Set>
/**
 * @brief Function for initializing the data structures used for
 *        learning the skeleton in parallel.
 *
 * @param myEdges The array to be initialized with the edges
 *                to be handled by this processor.
 * @param allNeighbors The map containing the neighborhood set
 *                     for all the variables in the network.
 */
void
GlobalLearning<Data, Var, Set>::initializeLearning(
  std::vector<std::tuple<Var, Var, double>>& myEdges,
  std::unordered_map<Var, Set>& allNeighbors
) const
{
  // First, block decompose all the edges on all the processors
  auto n = this->m_allVars.size();
  auto totalPairs = (n * (n - 1)) / 2;
  auto batchSize = (totalPairs / this->m_comm.size()) + ((totalPairs % this->m_comm.size() != 0) ? 1 : 0);
  auto myOffset = std::min(this->m_comm.rank() * batchSize, totalPairs);
  auto mySize = static_cast<uint32_t>(std::min(batchSize, totalPairs - myOffset));

  auto vars = std::vector<Var>(this->m_allVars.begin(), this->m_allVars.end());
  auto k = 0u;
  auto currOffset = 0u;
  while (currOffset + (n - (k + 1)) <= static_cast<uint32_t>(myOffset)) {
    currOffset += (n - (k + 1));
    ++k;
  }
  auto primary = vars.begin() + k;
  auto secondary = primary + 1 + (myOffset - currOffset);

  myEdges.resize(mySize);
  for (auto i = 0u; i < mySize; ++i, ++secondary) {
    if (secondary == vars.end()) {
      // Start a new primary variable if the secondary variable is exhausted
      ++primary;
      secondary = primary + 1;
    }
    // Initialize the p-values
    myEdges[i] = std::make_tuple(*primary, *secondary, 0.0);
  }
  // Then, initialize all the neighbor sets
  for (const auto v : vars) {
    allNeighbors.insert(std::make_pair(v, set_init(this->getCandidates(v), this->m_data.numVars())));
  }
}

template <typename Data, typename Var, typename Set>
/**
 * @brief Function that synchronizes the candidate sets by taking an
 *        intersection of the sets across all the processors.
 *
 * @param mySets A map with the neighbors of the variables
 *               remaining on this processor.
 */
void
GlobalLearning<Data, Var, Set>::syncSets(
  std::unordered_map<Var, Set>& mySets
) const
{
  set_allintersect_indexed(mySets, this->m_allVars, this->m_data.numVars(), this->m_comm);
}

template <typename Data, typename Var, typename Set>
/**
 * @brief Function for storing removed edges if they might be relevant
 *        for directing edges later.
 *
 * @param myRemoved A vector of all the removed edges.
 * @param myDSepSets Vector containing d-separating sets for all the removed edges.
 * @param allNeighbors Neighbor sets for all the variables.
 */
void
GlobalLearning<Data, Var, Set>::storeRemovedEdges(
  std::vector<std::tuple<Var, Var, double>>&& myRemoved,
  std::vector<Set>&& myDSepSets,
  const std::unordered_map<Var, Set>& allNeighbors
) const
{
  LOG_MESSAGE_IF(myRemoved.size() != myDSepSets.size(),
                 error, "Mismatch between number of edges and d-separating sets.");
  // Only retain information which can be used for directing edges later
  auto e = myRemoved.begin();
  auto s = myDSepSets.begin();
  while (e != myRemoved.end()) {
    Var x, y;
    std::tie(x, y, std::ignore) = *e;
    if (set_intersection(allNeighbors.at(x), allNeighbors.at(y)).empty()) {
      e = myRemoved.erase(e);
      s = myDSepSets.erase(s);
    }
    else {
      ++e;
      ++s;
    }
  }
  LOG_MESSAGE_IF(myRemoved.size() != myDSepSets.size(),
                 error, "Mismatch between number of edges and d-separating sets.");
  // Gather this information from all the processes
  auto removedSizes = mxx::allgather(myRemoved.size(), this->m_comm);
  auto allRemoved = mxx::allgatherv(myRemoved, removedSizes, this->m_comm);
  auto allDSepSets = set_allgatherv(myDSepSets, removedSizes, this->m_data.numVars(), this->m_comm);
  // Store this information in expected format
  Var x, y;
  double pv;
  for (auto e = 0u; e < allRemoved.size(); ++e) {
    std::tie(x, y, pv) = allRemoved.at(e);
    m_removedEdges.insert(std::make_pair(std::make_pair(x, y),
                          std::make_pair(pv, allDSepSets.at(e))));
  }
}

template <typename Data, typename Var, typename Set>
/**
 * @brief Function for constructing skeleton from the neighbor sets.
 *
 * @param allNeighbors Neighborhood sets for all the variables.
 */
BayesianNetwork<Var>
GlobalLearning<Data, Var, Set>::constructSkeleton(
  const std::unordered_map<Var, Set>&& allNeighbors
) const
{
  BayesianNetwork<Var> bn(this->m_data.varNames(this->m_allVars));
  for (const auto& neighbors : allNeighbors) {
    const auto x = neighbors.first;
    for (const auto y : neighbors.second) {
      if (x < y) {
        LOG_MESSAGE(info, "+ Adding the edge %s <-> %s",
                          this->m_data.varName(x), this->m_data.varName(y));
        bn.addEdge(x, y, true);
      }
    }
  }
  m_cachedNeighbors = std::move(allNeighbors);
  return bn;
}

template <typename Data, typename Var, typename Set>
/**
 * @brief Checks if the given v-structure (y-x-z) is a collider.
 *        Also returns the p-value for the collider.
 *
 * This function checks if x is part of the set that rendered
 * y and z independent. If not, y-x-z is a collider.
 */
std::pair<bool, double>
GlobalLearning<Data, Var, Set>::checkCollider(
  const Var y,
  const Var x,
  const Var z
) const
{
  LOG_MESSAGE_IF(m_removedEdges.find(std::make_pair(y, z)) == m_removedEdges.end(),
                 error, "Unable to check collider %s - %s - %s",
                        this->m_data.varName(y), this->m_data.varName(x), this->m_data.varName(z));
  auto p = m_removedEdges.at(std::make_pair(y, z));
  auto collider = false;
  if (!p.second.contains(x)) {
    collider = true;
  }
  return std::make_pair(collider, p.first);
}

template <typename Data, typename Var, typename Set>
PCStable<Data, Var, Set>::PCStable(
  const mxx::comm& comm,
  const Data& data,
  const double alpha,
  const Var maxConditioning
) : GlobalLearning<Data, Var, Set>(comm, data, alpha, maxConditioning)
{
}

template <typename Data, typename Var, typename Set>
std::pair<double, Set>
PCStable<Data, Var, Set>::checkEdge(
  const std::tuple<Var, Var, double>& edge,
  std::unordered_map<Var, Set>& allNeighbors,
  std::unordered_map<Var, Set>& removedNeighbors,
  const uint32_t setSize
) const
{
  auto pv = 0.0;
  auto dsep = set_init(Set(), this->m_data.numVars());
  const auto x = std::get<0>(edge);
  const auto y = std::get<1>(edge);
  auto& xNeighbors = allNeighbors.at(x);
  auto& yNeighbors = allNeighbors.at(y);
  LOG_MESSAGE(debug, "Checking the edge %s <-> %s, d-separating set of size %u",
                     this->m_data.varName(x), this->m_data.varName(y), setSize);
  if ((xNeighbors.size() <= setSize) && (yNeighbors.size() <= setSize)) {
    // The neighborhoods of both x and y are too small to remove the edge now
    LOG_MESSAGE(debug, "+ Fixing the edge %s <-> %s",
                       this->m_data.varName(x), this->m_data.varName(y));
    // Mark it for removal so that it does not get tested in subsequent iterations
    pv = std::numeric_limits<double>::max();
  }
  else {
    auto remove = false;
    // First, check the edge using neighbors of x
    if (xNeighbors.size() > setSize) {
      xNeighbors.erase(y);
      std::tie(pv, dsep) = this->m_data.maxPValueSubset(this->m_alpha, x, y, xNeighbors, setSize, setSize);
      xNeighbors.insert(y);
      remove = this->m_data.isIndependent(this->m_alpha, pv);
    }
    // Then, check the edge using neighbors of y if all of the following hold:
    // 1. The neighbors of x did not remove it already
    // 2. The conditioning set used for checking is not empty
    // 3. The size of the neighborhood of y is greater than the conditioning set size
    if (!remove && (setSize > 0) && (yNeighbors.size() > setSize)) {
      yNeighbors.erase(x);
      // Further, only check if neighborhood of y has some elements which are
      // not present in the neighborhood of x
      if (!set_difference(yNeighbors, xNeighbors).empty()) {
        std::tie(pv, dsep) = this->m_data.maxPValueSubset(this->m_alpha, x, y, yNeighbors, setSize, setSize);
      }
      yNeighbors.insert(x);
      remove = this->m_data.isIndependent(this->m_alpha, pv);
    }
    LOG_MESSAGE(debug, "%s and %s are " + std::string(this->m_data.isIndependent(this->m_alpha, pv) ? "independent" : "dependent") +
                       " (p-value = %g)",
                       this->m_data.varName(x), this->m_data.varName(y), pv);
    LOG_MESSAGE_IF(remove, debug, "- Removing the edge %s <-> %s",
                                  this->m_data.varName(x), this->m_data.varName(y));
    if (remove) {
      // Mark y for removal from neighborhood of x
      auto it = removedNeighbors.find(x);
      if (it == removedNeighbors.end()) {
        it = removedNeighbors.insert(it, std::make_pair(x, set_init(Set(), this->m_data.numVars())));
      }
      it->second.insert(y);
      // Mark x for removal from neighborhood of y
      it = removedNeighbors.find(y);
      if (it == removedNeighbors.end()) {
        it = removedNeighbors.insert(it, std::make_pair(y, set_init(Set(), this->m_data.numVars())));
      }
      it->second.insert(x);
    }
  }
  return std::make_pair(pv, dsep);
}

template <typename Data, typename Var, typename Set>
BayesianNetwork<Var>
PCStable<Data, Var, Set>::getSkeleton_sequential(
  const bool directEdges
) const
{
  std::unordered_map<Var, Set> allNeighbors;
  std::vector<std::tuple<Var, Var, double>> allEdges;
  this->initializeLearning(allEdges, allNeighbors);
  auto maxSize = std::min(this->m_maxConditioning, static_cast<Var>(this->m_allVars.size() - 2));
  for (auto s = 0u; (s <= maxSize) && !allEdges.empty(); ++s) {
    LOG_MESSAGE(debug, "Testing %u edges using sets of size %u", allEdges.size(), s);
    TIMER_DECLARE(tIter);
    std::unordered_map<Var, Set> removedNeighbors;
    for (auto& e : allEdges) {
      auto result = this->checkEdge(e, allNeighbors, removedNeighbors, s);
      std::get<2>(e) = result.first;
      // We need to store the removed edges since the d-separating
      // set may be required for directing edges later
      if (directEdges &&
          std::isless(result.first, std::numeric_limits<double>::max()) &&
          this->m_data.isIndependent(this->m_alpha, result.first)) {
        Var x, y;
        std::tie(x, y, std::ignore) = e;
        TIMER_START(this->m_tRemoved);
        // We only use this information for directing colliders of type x-z-y
        // However, there is no need to store it if there are no candidates for z
        // i.e., no common neighbors between x and y
        if (!set_intersection(allNeighbors.at(x), allNeighbors.at(y)).empty()) {
          this->m_removedEdges.insert(std::make_pair(std::make_pair(x, y), result));
        }
        TIMER_PAUSE(this->m_tRemoved);
      }
    }
    auto newEnd = std::remove_if(allEdges.begin(), allEdges.end(),
                                 [this] (const std::tuple<Var, Var, double>& e)
                                        { return this->m_data.isIndependent(this->m_alpha, std::get<2>(e)); });
    allEdges.resize(std::distance(allEdges.begin(), newEnd));
    for (const auto& rn : removedNeighbors) {
      allNeighbors[rn.first] = set_difference(allNeighbors.at(rn.first), rn.second);
    }
    TIMER_ELAPSED("Time taken in testing all sets of size " + std::to_string(s) + ": ", tIter);
  }
  return this->constructSkeleton(std::move(allNeighbors));
}

template <typename Data, typename Var, typename Set>
BayesianNetwork<Var>
PCStable<Data, Var, Set>::getSkeleton_parallel(
  const bool directEdges,
  const double imbalanceThreshold
) const
{
  std::vector<std::tuple<Var, Var, double>> myEdges;
  std::unordered_map<Var, Set> allNeighbors;
  this->initializeLearning(myEdges, allNeighbors);
  // The following two data structures are used for getting
  // d-separating set and p-value for all the removed edges
  std::vector<std::tuple<Var, Var, double>> myRemoved;
  std::vector<Set> myDSepSets;
  auto maxSize = std::min(this->m_maxConditioning, static_cast<Var>(this->m_allVars.size() - 2));
  for (auto s = 0u; (s <= maxSize) && mxx::any_of(myEdges.size() > 0, this->m_comm); ++s) {
    LOG_MESSAGE(debug, "Testing %u edges using sets of size %u", myEdges.size(), s);
    TIMER_DECLARE(tIter);
    std::unordered_map<Var, Set> removedNeighbors;
    for (auto& e : myEdges) {
      auto result = this->checkEdge(e, allNeighbors, removedNeighbors, s);
      std::get<2>(e) = result.first;
      // We need to store the removed edges since the d-separating
      // set may be required for directing edges later
      if (directEdges &&
          std::isless(result.first, std::numeric_limits<double>::max()) &&
          this->m_data.isIndependent(this->m_alpha, result.first)) {
        TIMER_START(this->m_tRemoved);
        Var x, y;
        std::tie(x, y, std::ignore) = e;
        // We only use this information for directing colliders of type x-z-y
        // However, there is no need to store it if there are no candidates for z
        // i.e., no common neighbors between x and y
        if (!set_intersection(allNeighbors.at(x), allNeighbors.at(y)).empty()) {
          myRemoved.push_back(e);
          myDSepSets.push_back(result.second);
        }
        TIMER_PAUSE(this->m_tRemoved);
      }
    }
    auto newEnd = std::remove_if(myEdges.begin(), myEdges.end(),
                                 [this] (const std::tuple<Var, Var, double>& e)
                                        { return this->m_data.isIndependent(this->m_alpha, std::get<2>(e)); });
    myEdges.resize(std::distance(myEdges.begin(), newEnd));
    TIMER_START(this->m_tDist);
    this->fixImbalance(myEdges, imbalanceThreshold);
    TIMER_PAUSE(this->m_tDist);
    for (const auto& rn : removedNeighbors) {
      allNeighbors[rn.first] = set_difference(allNeighbors.at(rn.first), rn.second);
    }
    TIMER_START(this->m_tSync);
    this->syncSets(allNeighbors);
    TIMER_PAUSE(this->m_tSync);
    TIMER_ELAPSED("Time taken in testing all sets of size " + std::to_string(s) + ": ", tIter);
  }
  if (directEdges) {
    TIMER_START(this->m_tRemoved);
    this->storeRemovedEdges(std::move(myRemoved), std::move(myDSepSets), allNeighbors);
    TIMER_PAUSE(this->m_tRemoved);
  }
  return this->constructSkeleton(std::move(allNeighbors));
}

#endif // DETAIL_GLOBALLEARNING_HPP_
