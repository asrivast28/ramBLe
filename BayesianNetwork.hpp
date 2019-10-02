/**
 * @file BayesianNetwork.hpp
 * @brief Declaration of the BayesianNetwork class.
 */
#ifndef BAYESIANNETWORK_HPP_
#define BAYESIANNETWORK_HPP_

#include "graph/Graph.hpp"

#include <boost/graph/filtered_graph.hpp>


/**
 * @brief Class which provides functionality for a Bayesian network.
 *
 * @tparam VarType Type of variable indices (expected to be an integer type).
 */
template <typename VarType>
class BayesianNetwork : public Graph<BidirectionalAdjacencyList, VertexLabel, VarType> {
private:
  class AntiParallelEdgeFilter;
  class EdgeCycleCounter;

public:
  using Vertex = typename Graph<BidirectionalAdjacencyList, VertexLabel, VarType>::Vertex;
  using Edge = typename Graph<BidirectionalAdjacencyList, VertexLabel, VarType>::Edge;

private:
  using GraphImpl = typename Graph<BidirectionalAdjacencyList, VertexLabel, VarType>::GraphImpl;
  using FilteredGraph = Graph<GenericBoostGraph, boost::filtered_graph<GraphImpl, AntiParallelEdgeFilter>, VarType>;

public:
  BayesianNetwork(const std::vector<std::string>&);

  bool
  hasDirectedCycles() const;

  void
  breakDirectedCycles();

  bool
  applyMeekRules();

  ~BayesianNetwork() { }

private:
  Graph<GenericBoostGraph, boost::filtered_graph<GraphImpl, AntiParallelEdgeFilter>, VarType>
  filterAntiParallelEdges() const;

  std::unordered_map<Edge, size_t, typename Edge::Hash>
  countEdgeCycles() const;

  bool
  removeEdgeAcyclic(Edge&&);

  bool
  unshieldedColliderRule(const Vertex&, const Vertex&) const;

  bool
  acyclicityRule(const Vertex&, const Vertex&) const;

  bool
  hybridRule(const Vertex&, const Vertex&) const;

private:
  // View of the network with only directed edges
  Graph<GenericBoostGraph, boost::filtered_graph<GraphImpl, AntiParallelEdgeFilter>, VarType> m_directed;
}; // class BayesianNetwork

#include "detail/BayesianNetwork.hpp"

#endif // BAYESIANNETWORK_HPP_
