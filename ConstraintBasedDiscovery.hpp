/**
 * @file ConstraintBasedDiscovery.hpp
 * @brief Declaration of the ConstraintBasedDiscovery class.
 */
#ifndef CONSTRAINTBASEDDISCOVERY_HPP_
#define CONSTRAINTBASEDDISCOVERY_HPP_

#include "graph/Graph.hpp"

#include <set>
#include <unordered_map>


/**
 * @brief Abstract base class for causal discovery using constraint-based learning.
 *
 * @tparam DataType Type of the object which is used for querying the data.
 * @tparam VarType Type of variable indices (expected to be an integer type).
 * @tparam SetType Type of set container.
 */
template <typename DataType, typename VarType, typename SetType>
class ConstraintBasedDiscovery {
public:
  ConstraintBasedDiscovery(const DataType&);

  SetType
  getPC(const VarType) const;

  SetType
  getMB(const VarType) const;

  Graph<BidirectionalAdjacencyList, VertexLabel, VarType>
  getNetwork(const bool = false) const;

  virtual
  ~ConstraintBasedDiscovery() { };

protected:
  SetType
  getCandidates(const VarType) const;

  virtual
  SetType
  getCandidatePC(const VarType, SetType) const = 0;

  virtual
  SetType
  getCandidateMB(const VarType, SetType) const = 0;

private:
  SetType
  getCandidatePC_cache(const VarType, SetType) const;

  void
  symmetryCorrectPC(const VarType, SetType&) const;

  SetType
  getCandidateMB_cache(const VarType, SetType) const;

  void
  symmetryCorrectMB(const VarType, SetType&) const;

  bool
  isVStructure(const VarType, const VarType, const VarType) const;

  void
  addVarNeighbors(const VarType, Graph<BidirectionalAdjacencyList, VertexLabel, VarType>&, const bool) const;

protected:
  const DataType m_data;
  SetType m_allVars;

private:
  mutable std::unordered_map<VarType, SetType> m_cachedPC;
  mutable std::unordered_map<VarType, SetType> m_cachedMB;
}; // class ConstraintBasedDiscovery

#include "detail/ConstraintBasedDiscovery.hpp"

#endif // CONSTRAINTBASEDDISCOVERY_HPP_
