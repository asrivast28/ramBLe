/**
 * @file BayesianNetwork.hpp
 * @brief Implementation of the BayesianNetwork functions.
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
#ifndef DETAIL_BAYESIANNETWORK_HPP_
#define DETAIL_BAYESIANNETWORK_HPP_

#include "utils/Logging.hpp"

#include <boost/graph/copy.hpp>
#include <boost/graph/tiernan_all_cycles.hpp>


template <typename Var>
/**
  * @brief Helper class that implements the anti-parallel edge filter functionality.
 */
class BayesianNetwork<Var>::AntiParallelEdgeFilter {
public:
  AntiParallelEdgeFilter(
  ) : m_graph(nullptr)
  {
  }

  AntiParallelEdgeFilter(
    const GraphImpl& g
  ) : m_graph(&g)
  {
  }

  template <typename EdgeDescriptor>
  bool
  operator()(
    const EdgeDescriptor& e
  ) const
  {
    auto source = boost::source(e, *m_graph);
    auto target = boost::target(e, *m_graph);
    return !boost::edge(target, source, *m_graph).second;
  }

private:
  const GraphImpl* m_graph;
}; // class AntiParallelEdgeFilter

template <typename Var>
/**
 * @brief Helper class for counting all the simple cycles that
 *        an edge is part of.
 */
class BayesianNetwork<Var>::EdgeCycleCounter {
public:
  EdgeCycleCounter(
    const Graph<BidirectionalAdjacencyList, VertexLabel, Var>& graph,
    std::unordered_map<Edge, size_t, typename Edge::Hash>& counts
  ) : m_graph(graph),
      m_counts(counts)
  {
  }

  template <typename Path, typename DirectedGraph>
  void
  cycle(
    const Path& p,
    const DirectedGraph& dg
  )
  {
    using IndexMap = typename boost::property_map<DirectedGraph, boost::vertex_index_t>::const_type;
    IndexMap indices = boost::get(boost::vertex_index, dg);
    auto u = p.begin();
    auto v = u+1;
    while (v != p.end()) {
      auto e = m_graph.getEdge(boost::get(indices, *u), boost::get(indices, *v));
      if (m_counts.find(e) == m_counts.end()) {
        m_counts[e] = 0;
      }
      m_counts[e] += 1;
      ++u;
      ++v;
    }
    v = p.begin();
    auto e = m_graph.getEdge(boost::get(indices, *u), boost::get(indices, *v));
    if (m_counts.find(e) == m_counts.end()) {
      m_counts[e] = 0;
    }
    m_counts[e] += 1;
  }

private:
  const Graph<BidirectionalAdjacencyList, VertexLabel, Var>& m_graph;
  std::unordered_map<Edge, size_t, typename Edge::Hash>& m_counts;
}; // class EdgeCycleCounter

template <typename Var>
/**
 * @brief Constructs empty network with given labels as vertices.
 */
BayesianNetwork<Var>::BayesianNetwork(
  const std::vector<std::string>& varLabels
) : Graph<BidirectionalAdjacencyList, VertexLabel, Var>(varLabels),
    m_directed(this->filterAntiParallelEdges())
{
}

template <typename Var>
/**
 * @brief Adds a directed or undirected edge between the vertices.
 */
void
BayesianNetwork<Var>::addEdge(
  const Var source,
  const Var target,
  const bool undirected
)
{
  this->addEdge(source, target);
  if (undirected) {
    this->addEdge(target, source);
  }
}

template <typename Var>
/**
 * @brief Returns a filtered view of the current graph with the anti-parallel edges removed.
 */
typename BayesianNetwork<Var>::FilteredGraph
BayesianNetwork<Var>::filterAntiParallelEdges(
) const
{
  AntiParallelEdgeFilter bef(this->m_graph);
  boost::filtered_graph<decltype(this->m_graph), AntiParallelEdgeFilter> fg(this->m_graph, bef);
  return FilteredGraph(std::move(fg), this->m_idVertexMap);
}

template <typename Var>
/**
 * @brief Orients the edges of the network in accordance with v-structures.
 *
 * @param vStructures A tuple with all the discovered v-structures and the corresponding p-values.
 */
void
BayesianNetwork<Var>::applyVStructures(
  std::vector<std::tuple<double, Var, Var, Var>>&& vStructures
)
{
  // First sort the v-structures in the ascending order of the p-values
  std::sort(vStructures.begin(), vStructures.end());
  for (const auto& vs : vStructures) {
    auto y = this->wrap(this->m_idVertexMap.at(std::get<1>(vs)));
    auto x = this->wrap(this->m_idVertexMap.at(std::get<2>(vs)));
    auto z = this->wrap(this->m_idVertexMap.at(std::get<3>(vs)));
    // First check if the reverse edges still exist
    if (!this->edgeExists(y, x) || !this->edgeExists(z, x)) {
      LOG_MESSAGE(warning, "* Could not apply v-structure %s -> %s <- %s (p-value = % g)",
                           y.property().label, x.property().label, z.property().label, std::get<0>(vs));
      LOG_MESSAGE_IF(!this->edgeExists(y, x), debug, "* %s - %s has already been oriented in the opposite direction", y.property().label, x.property().label);
      LOG_MESSAGE_IF(!this->edgeExists(z, x), debug, "* %s - %s has already been oriented in the opposite direction", x.property().label, z.property().label);
      continue;
    }
    LOG_MESSAGE(info, "+ Applying the v-structure %s -> %s <- %s (p-value = %g)",
                      y.property().label, x.property().label, z.property().label, std::get<0>(vs));
    this->removeEdge(x, y);
    this->removeEdge(x, z);
  }
}

template <typename Var>
/**
 * @brief Function which checks if the network has directed cycles.
 */
bool
BayesianNetwork<Var>::hasDirectedCycles(
) const
{
  return m_directed.hasCycles();
}

template <typename Var>
/**
 * @brief Counts the number of simple cycles that each edge is part of.
 */
std::unordered_map<typename BayesianNetwork<Var>::Edge, size_t, typename BayesianNetwork<Var>::Edge::Hash>
BayesianNetwork<Var>::countEdgeCycles(
) const
{
  // Copy the directed view of the graph to a directed graph
  typename DirectedGraph<VertexLabel, Var>::Impl dg;
  boost::copy_graph(*m_directed, dg);
  // Record all the counts
  std::unordered_map<Edge, size_t, typename Edge::Hash> counts;
  EdgeCycleCounter ecc(*this, counts);
  boost::tiernan_all_cycles(dg, ecc);
  return counts;
}

template <typename Var>
/**
 * @brief Function which breaks directed cycles in the network by reversing
 *        the direction of the edge which is part of most cycles.
 */
void
BayesianNetwork<Var>::breakDirectedCycles(
)
{
  Edge e;
  auto maxCount = 0u;
  for (const auto& cc : this->countEdgeCycles()) {
    if (cc.second > maxCount) {
      e = cc.first;
      maxCount = cc.second;
    }
  }
  if (maxCount > 0) {
    auto source = *e.source();
    auto target = *e.target();
    // Reverse the edge
    LOG_MESSAGE(info, "* Reversing the direction of edge %s -> %s", e.source().property().label, e.target().property().label);
    this->removeEdge(source, target);
    this->addEdge(target, source);
  }
}

template <typename Var>
/**
 * @brief Function which orients an edge, if it doesn't create directed cycles.
 *
 * @param e The edge to be removed.
 *
 * @returns true if any changes were made, otherwise returns false.
 */
bool
BayesianNetwork<Var>::removeEdgeAcyclic(
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

template <typename Var>
/**
 * @brief Function which checks if the undirected edge Y - Z can be
 *        oriented as Y -> Z to prevent new unshielded colliders.
 */
bool
BayesianNetwork<Var>::unshieldedColliderRule(
  const Vertex& y,
  const Vertex& z
) const
{
  if (m_directed.wrap(*y).inDegree() == 0) {
    return false;
  }
  bool potential = false;
  // Examine all the edges incoming into Y for a potential X
  for (const auto inY : m_directed.wrap(*y).inEdges()) {
    auto x = inY.source();
    // Check if an X exists such that no edge exists between X and Z
    if (!(this->edgeExists(*z, *x) || this->edgeExists(*x, *z))) {
      potential = true;
      break;
    }
  }
  return potential;
}

template <typename Var>
/**
 * @brief Function which checks if the undirected edge X - Z can be
 *        oriented as X -> Z to prevent directed cycles.
 */
bool
BayesianNetwork<Var>::acyclicityRule(
  const Vertex& x,
  const Vertex& z
) const
{
  // Therefore, try to apply the rule by setting the source as Z and the target as X
  bool orientEdge = false;
  // Iterate over all the outgoing neighbors of X for a potential Y
  for (const auto y : m_directed.wrap(*x).outNeighbors()) {
    if (m_directed.edgeExists(*y, *z)) {
      // Orient this edge as X -> Z
      // ...as long as it does not create an immorality
      orientEdge = true;
      break;
    }
  }
  bool immorality = false;
  if (orientEdge) {
    for (const auto inZ : m_directed.wrap(*z).inEdges()) {
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

template <typename Var>
/**
 * @brief Function which checks if the undirected edge X - Z can be
 *        oriented as X -> Z by applying the hybrid rule.
 */
bool
BayesianNetwork<Var>::hybridRule(
  const Vertex& x,
  const Vertex& z
) const
{
  auto countY = 0u;
  // Iterate over all the incoming neighbors of Z for potential Y
  for (const auto inZ : m_directed.wrap(*z).inEdges()) {
    auto y = inZ.source();
    // Check if an undirected edge exists between X and Y
    if (this->edgeExists(*x, *y) && this->edgeExists(*y, *x)) {
      ++countY;
    }
  }
  // The rule can be applied only if at least two Ys were found
  return (countY >= 2);
}

template <typename Var>
/**
 * @brief Top level function for orienting edges using Meek's rules.
 *
 * @returns true if any changes were made, otherwise returns false.
 */
bool
BayesianNetwork<Var>::applyMeekRules(
)
{
  bool changed = false;
  auto isCollider = [] (const Vertex& v) { return (v.inDegree() > v.outDegree()) &&
                                                  (v.inDegree() - v.outDegree() > 1); };
  // Iterate over all the undirected edges
  for (auto e : this->antiParallelEdges()) {
    // Check if the anti-parallel edge still exists
    if (!e.hasAntiParallel()) {
      continue;
    }
    // See if this direction of the anti-parallel edge can be removed
    // Therefore, check if the rules can be applied to the reverse direction
    auto source = e.source();
    auto target = e.target();
    if (isCollider(source) && isCollider(target)) {
      LOG_MESSAGE_IF(*source < *target, info, "* Fixing edge %s - %s because of conflicting v-structures", source.property().label, target.property().label);
      continue;
    }
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

template <typename Var>
/**
 * @brief Top level function for writing the network in graphviz format.
 *
 * @param fileName Name of the file to which the network should be written.
 */
void
BayesianNetwork<Var>::writeGraphviz(
  const std::string& fileName
) const
{
  std::ofstream out(fileName);
  // Check if there are any directed edges in the graph
  // If not, write the graph as an undirected graph
  auto directed = false;
  for (const auto e : m_directed.edges()) {
    std::ignore = e;
    directed = true;
    break;
  }
  out << (directed ? "digraph" : "graph") << " {" << std::endl;
  for (const auto v : this->vertices()) {
    out << "  ";
    out << boost::escape_dot_string(v.property().label);
    out << " ;" << std::endl;
  }
  auto delimiter = directed ? " -> " : " -- ";
  for (const auto e : this->edges()) {
    bool write = false;
    if (directed && !this->edgeExists(e.target(), e.source())) {
      out << "  edge [dir=forward] ";
      write = true;
    }
    else if (e.source() < e.target()) {
      out << "  edge [dir=none] ";
      write = true;
    }
    if (write) {
      out << boost::escape_dot_string(e.source().property().label);
      out << delimiter;
      out << boost::escape_dot_string(e.target().property().label);
      out << " ;" << std::endl;
    }
  }
  out << "}" << std::endl;
}

#endif // DETAIL_BAYESIANNETWORK_HPP
