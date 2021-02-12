/**
 * @file Module.hpp
 * @brief Implementation of the module learning class.
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
#ifndef DETAIL_MODULE_HPP_
#define DETAIL_MODULE_HPP_

#include "TreeNode.hpp"


template <typename Data, typename Var, typename Set>
class Module {
public:
  Module(const Set&&, const mxx::comm&, const Data&);

  const Set&
  variables() const;

  void
  learnTreeStructures(const std::list<std::list<Set>>&&, const bool, const double);

  void
  learnParents(std::mt19937&, const Set&, const OptimalBeta&, const uint32_t);

  const std::unordered_map<Var, double>&
  allParents() const;

  const std::unordered_map<Var, double>&
  randParents() const;

  template <typename Stream>
  void
  toXML(Stream&, const uint32_t) const;

private:
  std::shared_ptr<TreeNode<Data, Var, Set>>
  bestOrderedMerge(const std::list<std::shared_ptr<TreeNode<Data, Var, Set>>>&, const bool) const;

private:
  std::list<std::shared_ptr<TreeNode<Data, Var, Set>>> m_trees;
  std::unordered_map<Var, double> m_allParents;
  std::unordered_map<Var, double> m_randParents;
  const Set m_variables;
  const mxx::comm& m_comm;
  const Data& m_data;
}; // class Module

template <typename Data, typename Var, typename Set>
Module<Data, Var, Set>::Module(
  const Set&& variables,
  const mxx::comm& comm,
  const Data& data
) : m_trees(),
    m_allParents(),
    m_randParents(),
    m_variables(variables),
    m_comm(comm),
    m_data(data)
{
}

template <typename Data, typename Var, typename Set>
const Set&
Module<Data, Var, Set>::variables(
) const
{
  return m_variables;
}

template <typename Data, typename Var, typename Set>
std::shared_ptr<TreeNode<Data, Var, Set>>
Module<Data, Var, Set>::bestOrderedMerge(
  const std::list<std::shared_ptr<TreeNode<Data, Var, Set>>>& treeList,
  const bool scoreBHC
) const
{
  LOG_MESSAGE_IF(treeList.empty(), error, "Empty tree list passed to ordered merge");
  if (treeList.size() > 1) {
    std::shared_ptr<TreeNode<Data, Var, Set>> bestMerged;
    double bestScore = std::numeric_limits<double>::lowest();
    for (auto fit = treeList.begin(), sit = std::next(fit); sit != treeList.end(); ++fit, ++sit) {
      auto mergedTree = std::make_shared<TreeNode<Data, Var, Set>>(*fit, *sit);
      auto mergeScore = mergedTree->mergeScore(scoreBHC);
      if ((bestScore == std::numeric_limits<double>::lowest()) ||
          std::isgreater(mergeScore, bestScore) ||
          (mergeScore == std::numeric_limits<double>::infinity())) {
        bestMerged = mergedTree;
        bestScore = mergeScore;
      }
    }
    LOG_MESSAGE_IF(!std::isfinite(bestScore), warning, "Performed ordered merge using non-finite merge score");
    LOG_MESSAGE(debug, "Best ordered merge score %g", bestScore);
    return bestMerged;
  }
  else {
    return treeList.front();
  }
}

template <typename Data, typename Var, typename Set>
void
Module<Data, Var, Set>::learnTreeStructures(
  const std::list<std::list<Set>>&& sampledClusters,
  const bool scoreBHC,
  const double scoreGain
)
{
  static auto treeCompare = [] (const std::shared_ptr<TreeNode<Data, Var, Set>>& a, const std::shared_ptr<TreeNode<Data, Var, Set>>& b)
                               { return a->mean() < b->mean(); };
  // Learn hierarchical trees for all the sampled observation clusters
  for (const auto& obsClusters : sampledClusters) {
    std::list<std::shared_ptr<TreeNode<Data, Var, Set>>> treeList;
    for (const auto& cluster : obsClusters) {
      auto newTree = std::make_shared<TreeNode<Data, Var, Set>>(this->m_data, m_variables, cluster);
      auto it = std::upper_bound(treeList.begin(), treeList.end(), newTree, treeCompare);
      treeList.insert(it, newTree);
    }
    while (treeList.size() > 2) {
      auto mergedTree = this->bestOrderedMerge(treeList, scoreBHC);
      const auto& mergedChildren = mergedTree->children();
      // Remove both the trees to be merged from the list
      auto firstTree = std::find(treeList.begin(), treeList.end(), mergedChildren.first);
      treeList.erase(firstTree, std::next(firstTree, 2));
      auto it = std::upper_bound(treeList.begin(), treeList.end(), mergedTree, treeCompare);
      // Add the merged merged tree to the list of trees
      treeList.insert(it, mergedTree);
    }
    // Create a tree by merging the remaining two trees and store it
    m_trees.push_back(this->bestOrderedMerge(treeList, scoreBHC));
    m_trees.back()->prune(scoreGain);
  }
}

template <typename Data, typename Var, typename Set>
void
Module<Data, Var, Set>::learnParents(
  std::mt19937& generator,
  const Set& candidateParents,
  const OptimalBeta& ob,
  const uint32_t numSplits
)
{
  for (auto& tree : m_trees) {
    for (auto* node : tree->nodes()) {
      node->learnParentsSplits(generator, candidateParents, ob, numSplits);
      auto nodeObs = node->observations().size();
      auto factor = static_cast<double>(nodeObs) / m_data.numObs();
      for (const auto& w : node->weightSplits()) {
        auto var = std::get<0>(w);
        auto wit = m_allParents.find(var);
        if (wit == m_allParents.end()) {
          wit = m_allParents.insert(std::make_pair(var, 0.0)).first;
        }
        wit->second += exp(std::get<2>(w) / nodeObs) * factor;
      }
      for (const auto& r : node->randomSplits()) {
        auto var = std::get<0>(r);
        auto rit = m_randParents.find(var);
        if (rit == m_randParents.end()) {
          rit = m_randParents.insert(std::make_pair(var, 0.0)).first;
        }
        rit->second += exp(std::get<2>(r) / nodeObs) * factor;
      }
    }
  }
}

template <typename Data, typename Var, typename Set>
const std::unordered_map<Var, double>&
Module<Data, Var, Set>::allParents(
) const
{
  return m_allParents;
}

template <typename Data, typename Var, typename Set>
const std::unordered_map<Var, double>&
Module<Data, Var, Set>::randParents(
) const
{
  return m_randParents;
}

template <typename Data, typename Var, typename Set>
template <typename Stream>
void
Module<Data, Var, Set>::toXML(
  Stream& stream,
  const uint32_t index
) const
{
  stream << "<Module id=\"" << index << "\">" << std::endl;
  stream << "<RegulationTrees>" << std::endl;
  for (const auto& tree : m_trees) {
    stream << "<Root>" << std::endl;
    tree->toXML(stream);
    stream << "</Root>" << std::endl;
  }
  stream << "</RegulationTrees>" << std::endl;
  stream << "</Module>" << std::endl;
}

#endif // DETAIL_MODULE_HPP_
