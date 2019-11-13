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
 * @tparam Var Type of variable indices (expected to be an integer type).
 */
template <typename Var>
class BayesianNetwork : public Graph<BidirectionalAdjacencyList, VertexLabel, Var> {
private:
  class AntiParallelEdgeFilter;
  class EdgeCycleCounter;

public:
  using Vertex = typename Graph<BidirectionalAdjacencyList, VertexLabel, Var>::Vertex;
  using Edge = typename Graph<BidirectionalAdjacencyList, VertexLabel, Var>::Edge;

private:
  using GraphImpl = typename Graph<BidirectionalAdjacencyList, VertexLabel, Var>::GraphImpl;
  using FilteredGraph = Graph<GenericBoostGraph, boost::filtered_graph<GraphImpl, AntiParallelEdgeFilter>, Var>;

public:
  BayesianNetwork(const std::vector<std::string>&);

  void
  applyVStructures(const std::multimap<Var, std::pair<Var, Var>>&);

  bool
  hasDirectedCycles() const;

  void
  breakDirectedCycles();

  bool
  applyMeekRules();

  void
  writeGraphviz(const std::string&, const bool) const;

  ~BayesianNetwork() { }

private:
  Graph<GenericBoostGraph, boost::filtered_graph<GraphImpl, AntiParallelEdgeFilter>, Var>
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
  Graph<GenericBoostGraph, boost::filtered_graph<GraphImpl, AntiParallelEdgeFilter>, Var> m_directed;
}; // class BayesianNetwork

#include "detail/BayesianNetwork.hpp"

#endif // BAYESIANNETWORK_HPP_
