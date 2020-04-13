/**
 * @file BayesianNetwork.hpp
 * @brief Declaration of the BayesianNetwork class.
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
  applyVStructures(std::vector<std::tuple<double, Var, Var, Var>>&&);

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
