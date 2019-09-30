/**
 * @file BayesianNetwork.hpp
 * @brief Implementation of the BayesianNetwork functions.
 */
#ifndef DETAIL_BAYESIANNETWORK_HPP_
#define DETAIL_BAYESIANNETWORK_HPP_

#include "utils/Logging.hpp"


template <typename VarType>
/**
 * @brief Constructs empty network with given labels as vertices.
 */
BayesianNetwork<VarType>::BayesianNetwork(
  const std::vector<std::string>& varLabels
) : Graph<BidirectionalAdjacencyList, VertexLabel, VarType>(varLabels),
    m_directed(this->filterAntiParallelEdges())
{
}

template <typename VarType>
/**
 * @brief Function which checks if the network has directed cycles.
 */
bool
BayesianNetwork<VarType>::hasDirectedCycles(
) const
{
  return m_directed.hasCycles();
}

template <typename VarType>
/**
 * @brief Function which orients an edge, if it doesn't create directed cycles.
 *
 * @param e The edge to be removed.
 *
 * @returns true if any changes were made, otherwise returns false.
 */
bool
BayesianNetwork<VarType>::removeEdgeAcyclic(
  Edge&& e
)
{
  auto source = e.source();
  auto target = e.target();
  this->removeEdge(std::move(e));
  if (m_directed.hasCycles(*source)) {
    this->addEdge(source, target);
    return false;
  }
  return true;
}

template <typename VarType>
/**
 * @brief Function which checks if the undirected edge Y - Z can be
 *        oriented as Y -> Z to prevent new unshielded colliders.
 */
bool
BayesianNetwork<VarType>::unshieldedColliderRule(
  const Vertex& y,
  const Vertex& z
) const
{
  if (m_directed.wrap(*y).inDegree() == 0) {
    return false;
  }
  bool potential = false;
  // Examine all the edges incoming into Y for a potential X
  for (const auto inY: m_directed.wrap(*y).inEdges()) {
    auto x = inY.source();
    // Check if an X exists such that no edge exists between X and Z
    if (!(this->edgeExists(*z, *x) || this->edgeExists(*x, *z))) {
      potential = true;
      break;
    }
  }
  return potential;
}

template <typename VarType>
/**
 * @brief Function which checks if the undirected edge X - Z can be
 *        oriented as X -> Z to prevent directed cycles.
 */
bool
BayesianNetwork<VarType>::acyclicityRule(
  const Vertex& x,
  const Vertex& z
) const
{
  // Therefore, try to apply the rule by setting the source as Z and the target as X
  bool orientEdge = false;
  // Iterate over all the outgoing neighbors of X for a potential Y
  for (const auto y: m_directed.wrap(*x).outNeighbors()) {
    if (m_directed.edgeExists(*y, *z)) {
      // Orient this edge as X -> Z
      // ...as long as it does not create an immorality
      orientEdge = true;
      break;
    }
  }
  bool immorality = false;
  if (orientEdge) {
    for (const auto inZ: m_directed.wrap(*z).inEdges()) {
      auto w = inZ.source();
      // Check if there is an edge between the source of another incoming edge into Z and X
      if (!(this->edgeExists(*x, *w) || this->edgeExists(*w, *x))) {
        // If no such edge exists, then the rule can not be applied because an immorality will be created
        immorality = true;
        break;
      }
    }
  }
  return (orientEdge && !immorality);
}

template <typename VarType>
/**
 * @brief Function which checks if the undirected edge X - Z can be
 *        oriented as X -> Z by applying the hybrid rule.
 */
bool
BayesianNetwork<VarType>::hybridRule(
  const Vertex& x,
  const Vertex& z
) const
{
  auto countY = 0u;
  // Iterate over all the incoming neighbors of Z for potential Y
  for (const auto inZ: m_directed.wrap(*z).inEdges()) {
    auto y = inZ.source();
    // Check if an undirected edge exists between X and Y
    if (this->edgeExists(*x, *y) && this->edgeExists(*y, *x)) {
      ++countY;
    }
  }
  // The rule can be applied only if at least two Ys were found
  return (countY >= 2);
}

template <typename VarType>
/**
 * @brief Top level function for orienting edges using Meek's rules.
 *
 * @returns true if any changes were made, otherwise returns false.
 */
bool
BayesianNetwork<VarType>::applyMeekRules(
)
{
  bool changed = false;
  // Iterate over all the undirected edges
  for (auto e: this->antiParallelEdges()) {
    // Check if the anti-parallel edge still exists
    if (!e.hasAntiParallel()) {
      continue;
    }
    // See if this direction of the anti-parallel edge can be removed
    // Therefore, check if the rules can be applied to the reverse direction
    auto source = e.source();
    auto target = e.target();
    if (this->unshieldedColliderRule(target, source)) { // Apply Meek's Rule 1
      if (this->removeEdgeAcyclic(std::move(e))) {
        LOG_MESSAGE(info, "* Directing edge %s -> %s (R1: unshielded colliders)", target.property().label, source.property().label);
        changed = true;
      }
    }
    else if (this->acyclicityRule(target, source)) { // Apply Meek's Rule 2
      if (this->removeEdgeAcyclic(std::move(e))) {
        LOG_MESSAGE(info, "* Directing edge %s -> %s (R2: acyclicity)", target.property().label, source.property().label);
        changed = true;
      }
    }
    else if (this->hybridRule(target, source)) { // Apply Meek's Rule 3
      if (this->removeEdgeAcyclic(std::move(e))) {
        LOG_MESSAGE(info, "* Directing edge %s -> %s (R3: hybrid)", target.property().label, source.property().label);
        changed = true;
      }
    }
  }
  return changed;
}

#endif // DETAIL_BAYESIANNETWORK_HPP_
