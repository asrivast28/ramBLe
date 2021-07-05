/**
 * @file ConstraintBasedLearning.hpp
 * @brief Implementation of the classes for constraint-based learning.
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
#ifndef DETAIL_CONSTRAINTBASEDLEARNING_HPP_
#define DETAIL_CONSTRAINTBASEDLEARNING_HPP_

#include "common/SetUtils.hpp"
#include "mxx/distribution.hpp"
#include "mxx/reduction.hpp"
#include "utils/Logging.hpp"
#include "utils/Timer.hpp"


NotImplementedError::NotImplementedError(
  const std::string& message
) : std::logic_error(message)
{
}

template <typename Data, typename Var, typename Set>
/**
 * @brief Constructs the object with the given data.
 *
 * @param data Reference to an object of the Data.
 */
ConstraintBasedLearning<Data, Var, Set>::ConstraintBasedLearning(
  const mxx::comm& comm,
  const Data& data,
  const double alpha,
  const Var maxConditioning
) : m_comm(comm),
    m_data(data),
    m_allVars(set_init(Set(), data.numVars())),
    m_alpha(alpha),
    m_maxConditioning(maxConditioning)
{
  for (auto i = 0u; i < data.numVars(); ++i) {
    m_allVars.insert(m_allVars.end(), i);
  }
  TIMER_RESET(m_tMxx);
  TIMER_RESET(m_tDirect);
}

template <typename Data, typename Var, typename Set>
/**
 * @brief Default destructor.
 */
ConstraintBasedLearning<Data, Var, Set>::~ConstraintBasedLearning(
)
{
#if TIMER
  auto myMxx = m_tMxx.elapsed();
  auto maxMxx= mxx::reduce(myMxx, 0, mxx::max<float>(), m_comm);
  if (m_comm.is_first()) {
    std::cout << "Time taken in mxx calls: " << maxMxx << std::endl;
    TIMER_ELAPSED_NONZERO("Time taken in directing the edges: ", m_tDirect);
  }
#endif
}

template <typename Data, typename Var, typename Set>
/**
 * @brief Function for getting all the candidates for the given target variable.
 *
 * @param target The index of the target variable.
 *
 * @return The indices of all the variables except the target.
 */
Set
ConstraintBasedLearning<Data, Var, Set>::getCandidates(
  const Var target
) const
{
  auto candidates = m_allVars;
  candidates.erase(target);
  return candidates;
}

template <typename Data, typename Var, typename Set>
/**
 * @brief Function that synchronizes the candidate sets by taking a
 *        union of the sets across all the processors.
 *
 * @param mySets A map with the candidate MB or PC sets of the
 *               primary variables on this processor.
 */
void
ConstraintBasedLearning<Data, Var, Set>::syncSets(
  std::unordered_map<Var, Set>& mySets
) const
{
  TIMER_START(this->m_tMxx);
  set_allunion_indexed(mySets, this->m_allVars, this->m_data.numVars(), this->m_comm);
  TIMER_PAUSE(this->m_tMxx);
}

template <typename Data, typename Var, typename Set>
/**
 * @brief Function that fixes the imbalance in the p-value list across processors.
 *
 * @param myPairs A list of tuples corresponding to all the local pairs.
 * @param imbalanceThreshold The amount of imbalance that can be tolerated.
 *
 * @return true if the list was redistributed to fix the imbalance.
 */
bool
ConstraintBasedLearning<Data, Var, Set>::fixImbalance(
  std::vector<std::tuple<Var, Var, double>>& myPairs,
  const double imbalanceThreshold
) const
{
  bool fixed = false;
  TIMER_START(m_tMxx);
  auto totalSize = mxx::allreduce(myPairs.size(), m_comm);
  TIMER_PAUSE(m_tMxx);
  auto imbalance = 0.0;
  if (totalSize > 0) {
    TIMER_START(m_tMxx);
    auto maxSize = mxx::allreduce(myPairs.size(), mxx::max<size_t>(), m_comm);
    TIMER_PAUSE(m_tMxx);
    auto avgSize = totalSize / static_cast<double>(m_comm.size());
    imbalance = (maxSize / avgSize) - 1.0;
    //if (m_comm.is_first()) {
      //std::cout << "Imbalance: " << imbalance << std::endl;
    //}
  }
  if (std::isgreater(imbalance, imbalanceThreshold)) {
    // Redistribute the pairs to fix the imbalance
    TIMER_START(m_tMxx);
    mxx::stable_distribute_inplace(myPairs, this->m_comm);
    TIMER_PAUSE(m_tMxx);
    fixed = true;
  }
  return fixed;
}

template <typename Data, typename Var, typename Set>
/**
 * @brief Finds all the potential v-structures centered
 *        on the given variable.
 */
std::vector<std::tuple<double, Var, Var, Var>>
ConstraintBasedLearning<Data, Var, Set>::findVStructures(
  const Var target
) const
{
  std::vector<std::tuple<double, Var, Var, Var>> vStructures;
  std::set<std::pair<Var, Var>> checked;
  const auto& pcTarget = this->getPC(target);
  for (const auto y : pcTarget) {
    // Candidate parents of target, which are not connected to y
    const auto& pcY = this->getPC(y);
    auto cpaTarget = set_difference(pcTarget, pcY);
    cpaTarget.erase(y);
    for (const auto z : cpaTarget) {
      auto curr = (y < z) ? std::make_pair(y, z) : std::make_pair(z, y);
      if (checked.find(curr) == checked.end()) {
        checked.insert(curr);
        auto res = this->checkCollider(curr.first, target, curr.second);
        if (res.first) {
          LOG_MESSAGE(info, "* Found new v-structure %s -> %s <- %s (p-value = %g)",
                            m_data.varName(curr.first), m_data.varName(target), m_data.varName(curr.second), res.second);
          vStructures.push_back(std::make_tuple(res.second, curr.first, target, curr.second));
        }
        LOG_MESSAGE_IF(!res.first, debug,
                       "* Rejected the v-structure %s -> %s <- %s (p-value = %g)",
                       m_data.varName(curr.first), m_data.varName(target), m_data.varName(curr.second), res.second);
      }
    }
  }
  return vStructures;
}

template <typename Data, typename Var, typename Set>
/**
 * @brief Finds all the potential v-structures in the network
 *        by trying to find v-structures centered on each variable.
 */
std::vector<std::tuple<double, Var, Var, Var>>
ConstraintBasedLearning<Data, Var, Set>::findVStructures(
  const bool isParallel
) const
{
  mxx::blk_dist dist(m_allVars.size(), m_comm);
  auto myOffset = isParallel ? dist.eprefix_size() : 0;
  auto mySize = isParallel ? dist.local_size() : dist.global_size();
  auto x = std::next(m_allVars.begin(), myOffset);
  std::vector<std::tuple<double, Var, Var, Var>> myVStructures;
  for (auto v = 0u; v < mySize; ++v, ++x) {
    auto xVStructures = this->findVStructures(*x);
    myVStructures.insert(myVStructures.end(), xVStructures.begin(), xVStructures.end());
  }
  if (isParallel) {
    TIMER_START(m_tMxx);
    auto vStructures = mxx::allgatherv(myVStructures, m_comm);
    TIMER_PAUSE(m_tMxx);
    return vStructures;
  }
  else {
    return myVStructures;
  }
}

template <typename Data, typename Var, typename Set>
/**
 * @brief Top level function for getting the Bayesian network.
 *
 * @param directEdges Specifies if the edges of the network should be directed.
 * @param isParallel Specifies if the skeleton should be constructed in parallel.
 * @param imbalanceThreshold Specifies the amount of imbalance the parallel algorithm should tolerate.
 */
BayesianNetwork<Var>
ConstraintBasedLearning<Data, Var, Set>::getNetwork(
  const bool directEdges,
  const bool isParallel,
  const double imbalanceThreshold
) const
{
  auto bn = isParallel ? this->getSkeleton_parallel(directEdges, imbalanceThreshold) :
                         this->getSkeleton_sequential(directEdges);
  if (directEdges) {
    TIMER_START(m_tDirect);
    // First, orient the v-structures
    auto vStructures = this->findVStructures(isParallel);
    bn.applyVStructures(std::move(vStructures));
    // Then, break any directed cycles in the network
    LOG_MESSAGE_IF(m_comm.is_first() && bn.hasDirectedCycles(),
                   info, "* The initial network contains directed cycles");
    while (bn.hasDirectedCycles()) {
      bn.breakDirectedCycles();
    }
    // Finally, orient edges by applying Meek's rules
    bool changed = true;
    while (changed) {
      changed = bn.applyMeekRules();
    }
    TIMER_PAUSE(m_tDirect);
  }
  return bn;
}

#endif // DETAIL_CONSTRAINTBASEDLEARNING_HPP_
