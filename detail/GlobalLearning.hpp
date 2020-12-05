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

#include "../SetUtils.hpp"

#include "utils/Logging.hpp"
#include "../prettyprint.hpp"

#include <algorithm>


template <typename Data, typename Var, typename Set>
/**
 * @brief Constructs the object with the given data.
 */
GlobalLearning<Data, Var, Set>::GlobalLearning(
  const mxx::comm& comm,
  const Data& data,
  const Var maxConditioning
) : ConstraintBasedLearning<Data, Var, Set>(comm, data, maxConditioning),
    m_cachedNeighbors(),
    m_removedEdges()
{
}

template <typename Data, typename Var, typename Set>
/**
 * @brief Default destructor that prints out timing information.
 */
GlobalLearning<Data, Var, Set>::~GlobalLearning(
)
{
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
  const Var maxConditioning
) : GlobalLearning<Data, Var, Set>(comm, data, maxConditioning)
{
}

template <typename Data, typename Var, typename Set>
BayesianNetwork<Var>
PCStable<Data, Var, Set>::getSkeleton_sequential(
) const
{
  std::unordered_map<Var, Set> allNeighbors;
  auto n = this->m_allVars.size();
  // We start with nC2 edges, one between each pair of vertices
  std::vector<std::tuple<Var, Var, bool>> candidateEdges((n * (n - 1)) / 2);
  auto e = 0u;
  for (const auto x : this->m_allVars) {
    allNeighbors[x] = this->getCandidates(x);
    for (const auto y : allNeighbors.at(x)) {
      if (x < y) {
        candidateEdges[e++] = std::make_tuple(x, y, false);
      }
    }
  }
  auto maxSize = std::min(this->m_maxConditioning, static_cast<Var>(this->m_allVars.size() - 2));
  BayesianNetwork<Var> bn(this->m_data.varNames(this->m_allVars));
  for (auto s = 0u; (s < maxSize) && !candidateEdges.empty(); ++s) {
    std::unordered_map<Var, Set> removedNeighbors;
    for (auto& e : candidateEdges) {
      const auto x = std::get<0>(e);
      const auto y = std::get<1>(e);
      auto& xNeighbors = allNeighbors.at(x);
      auto& yNeighbors = allNeighbors.at(y);
      if ((xNeighbors.size() <= s) && (yNeighbors.size() <= s)) {
        // The neighborhoods of both x and y are too small to remove the edge now
        LOG_MESSAGE(info, "+ Fixing the edge %s <-> %s", this->m_data.varName(x), this->m_data.varName(y));
        bn.addEdge(x, y, true);
        // Mark it for removal so that it does not get tested in subsequent iterations
        std::get<2>(e) = true;
        continue;
      }
      LOG_MESSAGE(debug, "Investigating the edge %s <-> %s, d-separating set of size %u",
                         this->m_data.varName(x), this->m_data.varName(y), s);
      bool remove = false;
      auto pv = std::numeric_limits<double>::lowest();
      auto dsep = set_init(Set(), this->m_data.numVars());
      // First, check the edge using neighbors of x
      if (xNeighbors.size() > s) {
        xNeighbors.erase(y);
        std::tie(pv, dsep) = this->m_data.maxPValueSubset(x, y, xNeighbors, s, s);
        xNeighbors.insert(y);
        remove = this->m_data.isIndependent(pv);
      }
      // Then, check the edge using neighbors of y if all of the following hold:
      // 1. The neighbors of x did not remove it already
      // 2. The conditioning set used for checking is not empty
      // 3. The size of the neighborhood of y is greater than the conditioning set size
      if (!remove && (s > 0) && (yNeighbors.size() > s)) {
        yNeighbors.erase(x);
        // Further, only check if neighborhood of y has some elements which are
        // not present in the neighborhood of x
        if (!set_difference(yNeighbors, xNeighbors).empty()) {
          std::tie(pv, dsep) = this->m_data.maxPValueSubset(x, y, yNeighbors, s, s);
        }
        yNeighbors.insert(x);
        remove = this->m_data.isIndependent(pv);
      }
      LOG_MESSAGE(debug, "%s and %s are " + std::string(this->m_data.isIndependent(pv) ? "independent" : "dependent") +
                         " (p-value = %g)",
                         this->m_data.varName(x), this->m_data.varName(y), pv);
      LOG_MESSAGE_IF(remove, debug, "- Removing the edge %s <-> %s", this->m_data.varName(x), this->m_data.varName(y));
      if (remove) {
        this->m_removedEdges.insert(std::make_pair(std::make_pair(x, y), std::make_pair(pv, dsep)));
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
      std::get<2>(e) = remove;
    }
    auto newEnd = std::remove_if(candidateEdges.begin(), candidateEdges.end(),
                                 [] (std::tuple<Var, Var, bool>& e) { return std::get<2>(e); });
    candidateEdges.resize(std::distance(candidateEdges.begin(), newEnd));
    for (const auto& rn : removedNeighbors) {
      allNeighbors[rn.first] = set_difference(allNeighbors.at(rn.first), rn.second);
    }
  }
  for (const auto& e : candidateEdges) {
    LOG_MESSAGE(info, "+ Fixing the edge %s <-> %s", this->m_data.varName(std::get<0>(e)), this->m_data.varName(std::get<1>(e)));
    bn.addEdge(std::get<0>(e), std::get<1>(e), true);
  }
  this->m_cachedNeighbors = std::move(allNeighbors);
  return bn;
}

template <typename Data, typename Var, typename Set>
BayesianNetwork<Var>
PCStable<Data, Var, Set>::getSkeleton_parallel(
  const double
) const
{
  throw NotImplementedError("PCStable: Parallel algorithm is not implemented yet");
  return BayesianNetwork<Var>(this->m_data.varNames(this->m_allVars));
}

#endif // DETAIL_GLOBALLEARNING_HPP_
