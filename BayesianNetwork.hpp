/**
 * @file BayesianNetwork.hpp
 * @brief Declaration of the BayesianNetwork class.
 */
#ifndef BAYESIANNETWORK_HPP_
#define BAYESIANNETWORK_HPP_

#include "graph/Graph.hpp"


/**
 * @brief Class which provides functionality for a Bayesian network.
 *
 * @tparam VarType Type of variable indices (expected to be an integer type).
 */
template <typename VarType>
class BayesianNetwork : public Graph<BidirectionalAdjacencyList, VertexLabel, VarType> {
public:
  using Vertex = typename Graph<BidirectionalAdjacencyList, VertexLabel, VarType>::Vertex;
  using Edge = typename Graph<BidirectionalAdjacencyList, VertexLabel, VarType>::Edge;

public:
  BayesianNetwork(const std::vector<std::string>&);

  bool
  hasDirectedCycles() const;

  bool
  applyMeekRules();

  ~BayesianNetwork() { }

private:
  bool
  removeEdgeAcyclic(Edge&&);

  bool
  unshieldedColliderRule(const Vertex&, const Vertex&) const;

  bool
  acyclicityRule(const Vertex&, const Vertex&) const;

  bool
  hybridRule(const Vertex&, const Vertex&) const;

private:
  using AntiParallelEdgeFilter = typename Graph<BidirectionalAdjacencyList, VertexLabel, VarType>::AntiParallelEdgeFilter;
  using GraphImpl = typename Graph<BidirectionalAdjacencyList, VertexLabel, VarType>::GraphImpl;

private:
  // View of the network with only directed edges
  Graph<GenericBoostGraph, boost::filtered_graph<GraphImpl, AntiParallelEdgeFilter>, VarType> m_directed;
}; // class BayesianNetwork

#include "detail/BayesianNetwork.hpp"

#endif // BAYESIANNETWORK_HPP_
