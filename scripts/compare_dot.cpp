/**
 * @file compare_dot.cpp
 * @brief The script for comparing two graphs
 *        given in the form of dot files.
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
using Graph = boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, VertexProperties, EdgeProperties>;

/**
 * @brief Reads a graph from the given dot file.
 *
 * @param fileName The name of the dot file.
 *
 * @return An undirected graph corresponding to the dot file.
 */
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
 * @param g The given boost graph.
 *
 * @return A set containing the names of all the vertices in the graph.
 */
std::set<std::string>
getVertexNames(
  const Graph& g
)
{
  std::set<std::string> allNames;
  const auto& vertexNames = boost::get(boost::vertex_name, g);
  boost::graph_traits<Graph>::vertex_iterator it, end;
  for (boost::tie(it, end) = boost::vertices(g); it != end; ++it) {
    allNames.insert(vertexNames[*it]);
  }
  return allNames;
}

/**
 * @brief Compares the count and name of vertices in given boost graphs.
 *
 * @param first The first graph to be compared.
 * @param second The second graph to be compared.
 * @param verbose If the differences should be printed.
 *
 * @return true if the two graphs have the same vertices, otherwise false.
 */
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
 * @param first The first graph to be compared.
 * @param second The second graph to be compared.
 *
 * @return A set of pairs with the vertex names corresponding to the
 *         exclusive edges.
 */
std::set<std::pair<std::string, std::string>>
edgeDifference(
  const Graph& first,
  const Graph& second
)
{
  std::set<std::pair<std::string, std::string>> edgeDiff;
  const auto& secondVertices = boost::get(boost::vertex_name, second);
  std::unordered_map<std::string, typename boost::graph_traits<Graph>::vertex_descriptor> secondNameDescriptorMap;
  boost::graph_traits<Graph>::vertex_iterator vit, vend;
  for (boost::tie(vit, vend) = boost::vertices(second); vit != vend; ++vit) {
    secondNameDescriptorMap.insert(std::make_pair(secondVertices[*vit], *vit));
  }
  const auto& firstVertices = boost::get(boost::vertex_name, first);
  boost::graph_traits<Graph>::edge_iterator eit, eend;
  for (boost::tie(eit, eend) = boost::edges(first); eit != eend; ++eit) {
    auto secondSource = secondNameDescriptorMap.at(firstVertices[boost::source(*eit, first)]);
    auto secondTarget = secondNameDescriptorMap.at(firstVertices[boost::target(*eit, first)]);
    if (!boost::edge(secondSource, secondTarget, second).second) {
      edgeDiff.insert(std::make_pair(firstVertices[boost::source(*eit, first)], firstVertices[boost::source(*eit, first)]));
    }
  }
  return edgeDiff;
}

/**
 * @brief Compares the edges in the given boost graphs.
 *        This function assumes that the two graphs have the same vertex names.
 *
 * @param first The first graph to be compared.
 * @param second The second graph to be compared.
 * @param verbose If the differences should be printed.
 *
 * @return true if the two graphs have same edges, otherwise false.
 */
bool
compareEdges(
  const Graph& first,
  const Graph& second,
  const bool verbose
)
{
  auto firstOnly = edgeDifference(first, second);
  if (!firstOnly.empty()) {
    if (verbose) {
      std::cerr << "Edges found only in the first graph: " << std::endl;
      for (const auto& e : firstOnly) {
        std::cout << "(" << e.first << " -- " << e.second << ")" << std::endl;
      }
    }
  }

  auto secondOnly = edgeDifference(second, first);
  if (!secondOnly.empty()) {
    if (verbose) {
      std::cerr << "Edges found only in the second graph: " << std::endl;
      for (const auto& e : secondOnly) {
        std::cout << "(" << e.first << " -- " << e.second << ")" << std::endl;
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

int
main(
  int argc,
  char** argv
)
{
  if (argc < 3) {
    std::cerr << "Usage: ./compare_dot <first> <second> [-v]" << std::endl;
    return 1;
  }
  auto first = readDotFile(boost::filesystem::canonical(argv[1]).string());
  auto second = readDotFile(boost::filesystem::canonical(argv[2]).string());
  auto verbose = false;
  if (argc >= 4) {
    verbose = true;
  }

  bool same = compareVertices(first, second, verbose);
  if (same) {
    same = compareEdges(first, second, verbose);
  }

  return static_cast<int>(same);
}
