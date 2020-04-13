/**
 * @file compare_dot.cpp
 * @brief The script for comparing two graphs
 *        given in the form of dot files.
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
#include <boost/filesystem.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/isomorphism.hpp>

#include <fstream>
#include <iostream>


namespace boost {
  enum edge_dir_t { edge_dir };
  BOOST_INSTALL_PROPERTY(edge, dir);
}

using VertexProperties = boost::property<boost::vertex_name_t, std::string>;
using EdgeProperties = boost::property<boost::edge_dir_t, std::string>;
using DirectedGraph = boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS, VertexProperties, EdgeProperties>;
using UndirectedGraph = boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, VertexProperties, EdgeProperties>;

/**
 * @brief Reads a graph from the given dot file.
 *
 * @tparam Graph Type of the boost graph.
 * @param fileName The name of the dot file.
 *
 * @return An undirected graph corresponding to the dot file.
 */
template <typename Graph>
Graph
readDotFile(
  const std::string& fileName
)
{
  std::ifstream dotFile(fileName);
  Graph g(0);
  boost::dynamic_properties dp;
  dp.property("name", boost::get(boost::vertex_name, g));
  dp.property("dir", boost::get(boost::edge_dir, g));
  if (boost::read_graphviz(dotFile, g, dp, "name")) {
    return g;
  }
  else {
    throw std::runtime_error("Could not read graph from " + std::string(fileName));
  }
}

/**
 * @brief Gets the names of all the vertices in the given boost graph.
 *
 * @tparam Graph Type of the boost graph.
 * @param g The given boost graph.
 *
 * @return A set containing the names of all the vertices in the graph.
 */
template <typename Graph>
std::set<std::string>
getVertexNames(
  const Graph& g
)
{
  std::set<std::string> allNames;
  const auto& vertexNames = boost::get(boost::vertex_name, g);
  typename boost::graph_traits<Graph>::vertex_iterator it, end;
  for (boost::tie(it, end) = boost::vertices(g); it != end; ++it) {
    allNames.insert(vertexNames[*it]);
  }
  return allNames;
}

/**
 * @brief Compares the count and name of vertices in given boost graphs.
 *
 * @tparam Graph Type of the boost graph.
 * @param first The first graph to be compared.
 * @param second The second graph to be compared.
 * @param verbose If the differences should be printed.
 *
 * @return true if the two graphs have the same vertices, otherwise false.
 */
template <typename Graph>
bool
compareVertices(
  const Graph& first,
  const Graph& second,
  const bool verbose
)
{
  bool same = true;
  if (boost::num_vertices(first) != boost::num_vertices(second)) {
    same = false;
    if (verbose) {
      std::cerr << "First graph has " << boost::num_vertices(first) << " vertices while " <<
                   "second graph has " << boost::num_vertices(second) << "vertices" << std::endl;
    }
  }
  auto firstNames = getVertexNames(first);
  auto secondNames = getVertexNames(second);
  decltype(firstNames) missingNames;
  std::set_difference(firstNames.begin(), firstNames.end(),
                      secondNames.begin(), secondNames.end(),
                      std::inserter(missingNames, missingNames.end()));
  if (!missingNames.empty()) {
    same = false;
    if (verbose) {
      std::cerr << "Vertices found only in the first graph: " << std::endl;
      for (const auto& name : missingNames) {
        std::cerr << name << std::endl;
      }
    }
  }
  missingNames.clear();
  std::set_difference(secondNames.begin(), secondNames.end(),
                      firstNames.begin(), firstNames.end(),
                      std::inserter(missingNames, missingNames.end()));
  if (!missingNames.empty()) {
    same = false;
    if (verbose) {
      std::cerr << "Vertices found only in the second graph: " << std::endl;
      for (const auto& name : missingNames) {
        std::cerr << name << std::endl;
      }
    }
  }
  return same;
}

/**
 * @brief Finds the edges which are present in the first graph but
 *        not present in the second graph.
 *
 * @tparam Graph Type of the boost graph.
 * @param first The first graph to be compared.
 * @param second The second graph to be compared.
 *
 * @return A set of pairs with the vertex names corresponding to the
 *         exclusive edges.
 */
template <typename Graph>
std::set<std::pair<std::string, std::string>>
edgeDifference(
  const Graph& first,
  const Graph& second
)
{
  std::set<std::pair<std::string, std::string>> edgeDiff;
  const auto& secondVertices = boost::get(boost::vertex_name, second);
  std::unordered_map<std::string, typename boost::graph_traits<Graph>::vertex_descriptor> secondNameDescriptorMap;
  typename boost::graph_traits<Graph>::vertex_iterator vit, vend;
  for (boost::tie(vit, vend) = boost::vertices(second); vit != vend; ++vit) {
    secondNameDescriptorMap.insert(std::make_pair(secondVertices[*vit], *vit));
  }
  const auto& firstVertices = boost::get(boost::vertex_name, first);
  typename boost::graph_traits<Graph>::edge_iterator eit, eend;
  for (boost::tie(eit, eend) = boost::edges(first); eit != eend; ++eit) {
    auto secondSource = secondNameDescriptorMap.at(firstVertices[boost::source(*eit, first)]);
    auto secondTarget = secondNameDescriptorMap.at(firstVertices[boost::target(*eit, first)]);
    if (!boost::edge(secondSource, secondTarget, second).second) {
      edgeDiff.insert(std::make_pair(firstVertices[boost::source(*eit, first)], firstVertices[boost::target(*eit, first)]));
    }
  }
  return edgeDiff;
}

/**
 * @brief Compares the edges in the given boost graphs.
 *        This function assumes that the two graphs have the same vertex names.
 *
 * @tparam Graph Type of the boost graph.
 * @param first The first graph to be compared.
 * @param second The second graph to be compared.
 * @param verbose If the differences should be printed.
 *
 * @return true if the two graphs have same edges, otherwise false.
 */
template <typename Graph>
bool
compareEdges(
  const Graph& first,
  const Graph& second,
  const bool verbose
)
{
  auto separator = std::is_same<Graph, DirectedGraph>::value ? " -> " : " -- ";
  auto firstOnly = edgeDifference(first, second);
  if (!firstOnly.empty()) {
    if (verbose) {
      std::cerr << "Edges found only in the first graph: " << std::endl;
      for (const auto& e : firstOnly) {
        std::cerr << "(" << e.first << separator << e.second << ")" << std::endl;
      }
    }
  }

  auto secondOnly = edgeDifference(second, first);
  if (!secondOnly.empty()) {
    if (verbose) {
      std::cerr << "Edges found only in the second graph: " << std::endl;
      for (const auto& e : secondOnly) {
        std::cerr << "(" << e.first << separator << e.second << ")" << std::endl;
      }
    }
  }

  if ((boost::num_edges(first) - firstOnly.size()) != (boost::num_edges(second) - secondOnly.size())) {
    throw std::runtime_error("Something went wrong in the edge comparison");
  }

  std::cout << "Edge comparison results: " << std::endl;
  std::cout << "# of edges found only in the first graph: " << firstOnly.size() << std::endl;
  std::cout << "# of edges found only in the second graph: " << secondOnly.size() << std::endl;
  std::cout << "# of edges common to both: " << boost::num_edges(first) - firstOnly.size() << std::endl;

  return firstOnly.empty() && secondOnly.empty();
}

/**
 * @brief Compares two graphs given in dot files.
 *
 * @tparam Graph Type of the boost graph to be used for the comparison.
 * @param firstFile Name of the first dot file.
 * @param secondFile Name of the second dot file.
 * @param verbose If the differences should be printed.
 *
 * @return true if the two graphs have same nodes and edges, otherwise false.
 */
template <typename Graph>
bool
compareGraphs(
  const char* const firstFile,
  const char* const secondFile,
  const bool verbose
)
{
  auto firstGraph = readDotFile<Graph>(boost::filesystem::canonical(firstFile).string());
  auto secondGraph = readDotFile<Graph>(boost::filesystem::canonical(secondFile).string());
  bool same = compareVertices(firstGraph, secondGraph, verbose);
  if (same) {
    same = compareEdges(firstGraph, secondGraph, verbose);
  }
  return same;
}

int
main(
  int argc,
  char** argv
)
{
  if (argc < 3) {
    std::cerr << "Usage: ./compare_dot <first> <second> [-d] [-v]" << std::endl;
    return 1;
  }
  auto directed = false;
  auto verbose = false;
  if (argc >= 5) {
    verbose = true;
    directed = true;
  }
  else if (argc >= 4) {
    if (strcmp(argv[3], "-v") == 0) {
      verbose = true;
    }
    else {
      directed = true;
    }
  }

  auto result = false;
  if (directed) {
    result = compareGraphs<DirectedGraph>(argv[1], argv[2], verbose);
  }
  else {
    result = compareGraphs<UndirectedGraph>(argv[1], argv[2], verbose);
  }
  return static_cast<int>(!result);
}
