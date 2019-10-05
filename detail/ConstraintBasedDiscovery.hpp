/**
 * @file ConstraintBasedDiscovery.hpp
 * @brief Implementation of the ConstraintBasedDiscovery functions.
 */
#ifndef DETAIL_CONSTRAINTBASEDDISCOVERY_HPP_
#define DETAIL_CONSTRAINTBASEDDISCOVERY_HPP_

#include "../SetUtils.hpp"

#include "utils/Logging.hpp"


template <typename DataType, typename VarType, typename SetType>
/**
 * @brief Constructs the object with the given data.
 *
 * @param data Reference to an object of the DataType.
 */
ConstraintBasedDiscovery<DataType, VarType, SetType>::ConstraintBasedDiscovery(
  const DataType& data
) : m_data(data),
    m_allVars(set_init(SetType(), data.numVars()))
{
  for (auto i = 0u; i < data.numVars(); ++i) {
    m_allVars.insert(m_allVars.end(), i);
  }
}

template <typename DataType, typename VarType, typename SetType>
/**
 * @brief Function for getting all the candidates for the given target variable.
 *
 * @param target The index of the target variable.
 *
 * @return The indices of all the variables except the target.
 */
SetType
ConstraintBasedDiscovery<DataType, VarType, SetType>::getCandidates(
  const VarType target
) const
{
  auto candidates = m_allVars;
  candidates.erase(target);
  return candidates;
}

template <typename DataType, typename VarType, typename SetType>
/**
 * @brief Finds the candidate PC for a variable, and caches the result.
 *
 * @param target The index of the target variable.
 * @param candidates The indices of all the candidate variables.
 *
 * @return A pair with the first element being a set containing the indices of
 *         all the variables in the candidate PC of the given target variable,
 *         and the second element specifying if the set has been symmetry corrected.
 */
std::pair<SetType, bool>
ConstraintBasedDiscovery<DataType, VarType, SetType>::getCandidatePC_cache(
  const VarType target,
  SetType candidates
) const
{
  auto cacheIt = m_cachedPC.find(target);
  if (cacheIt == m_cachedPC.end()) {
    auto cpc = std::make_pair(this->getCandidatePC(target, std::move(candidates)), false);
    m_cachedPC.insert(cacheIt, std::make_pair(target, cpc));
    return cpc;
  }
  else {
    LOG_MESSAGE(trace, "* Found candidate PC for %s in the cache", this->m_data.varName(target))
    return cacheIt->second;
  }
}

template <typename DataType, typename VarType, typename SetType>
/**
 * @brief Performs symmetry correction on the given PC set for the given variable.
 *
 * @param target The index of the target variable.
 * @param cpc The set containing the indices of all the variables in
 *            the candidate PC set.
 */
void
ConstraintBasedDiscovery<DataType, VarType, SetType>::symmetryCorrectPC(
  const VarType target,
  SetType& cpc
) const
{
  auto initial = cpc;
  for (const VarType x: initial) {
    auto candidatesX = this->getCandidates(x);
    auto cpcX = this->getCandidatePC_cache(x, std::move(candidatesX)).first;
    if (!set_contains(cpcX, target)) {
      LOG_MESSAGE(info, "- Removing %s from the PC of %s (asymmetry)", this->m_data.varName(x), this->m_data.varName(target));
      cpc.erase(x);
    }
  }
}

template <typename DataType, typename VarType, typename SetType>
/**
 * @brief The top level function for getting the correct PC set
 *        for the given target variable.
 *
 * @param target The index of the target variable.
 *
 * @return A set containing the indices of all the variables
 *         in the correct PC set of the target variable.
 */
SetType
ConstraintBasedDiscovery<DataType, VarType, SetType>::getPC(
  const VarType target
) const
{
  auto candidates = this->getCandidates(target);
  auto cpc = this->getCandidatePC_cache(target, std::move(candidates));
  if (!cpc.second) {
    this->symmetryCorrectPC(target, cpc.first);
    cpc.second = true;
  }
  return cpc.first;
}

template <typename DataType, typename VarType, typename SetType>
/**
 * @brief Finds the candidate MB for a variable, and caches the result.
 *
 * @param target The index of the target variable.
 * @param candidates The indices of all the candidate variables.
 *
 * @return A pair with the first element being a set containing the indices of
 *         all the variables in the candidate MB of the given target variable,
 *         and the second element specifying if the set has been symmetry corrected.
 */
std::pair<SetType, bool>
ConstraintBasedDiscovery<DataType, VarType, SetType>::getCandidateMB_cache(
  const VarType target,
  SetType candidates
) const
{
  auto cacheIt = m_cachedMB.find(target);
  if (cacheIt == m_cachedMB.end()) {
    auto cmb = std::make_pair(this->getCandidateMB(target, std::move(candidates)), false);
    m_cachedMB.insert(cacheIt, std::make_pair(target, cmb));
    return cmb;
  }
  else {
    LOG_MESSAGE(trace, "* Found candidate MB for %s in the cache", this->m_data.varName(target));
    return cacheIt->second;
  }
}

template <typename DataType, typename VarType, typename SetType>
/**
 * @brief Symmetry corrects the candidate MB of the target variable.
 *
 * @param target The index of the target variable.
 * @param cmb The indices of the variables in the candidate MB of the target variable.
 *                The function removes the indices from the candidate set.
 */
void
ConstraintBasedDiscovery<DataType, VarType, SetType>::symmetryCorrectMB(
  const VarType target,
  SetType& cmb
) const
{
  auto initial = cmb;
  for (const VarType x: initial) {
    auto candidatesX = this->getCandidates(x);
    auto cmbX = this->getCandidateMB_cache(x, std::move(candidatesX)).first;
    if (!set_contains(cmbX, target)) {
      LOG_MESSAGE(info, "- Removing %s from the MB of %s (asymmetry)", this->m_data.varName(x), this->m_data.varName(target));
      cmb.erase(x);
    }
  }
}

template <typename DataType, typename VarType, typename SetType>
/**
 * @brief Top level function for getting the MB of the target variable.
 *
 * @param target The index of the target variable.
 */
SetType
ConstraintBasedDiscovery<DataType, VarType, SetType>::getMB(
  const VarType target
) const
{
  auto candidates = this->getCandidates(target);
  auto cmb = this->getCandidateMB_cache(target, std::move(candidates));
  if (!cmb.second) {
    this->symmetryCorrectMB(target, cmb.first);
    cmb.second = true;
  }
  return cmb.first;
}

template <typename DataType, typename VarType, typename SetType>
/**
 * @brief Checks if y-x-z forms a v-structure.
 */
bool
ConstraintBasedDiscovery<DataType, VarType, SetType>::isCollider(
  const VarType y,
  const VarType x,
  const VarType z
) const
{
  static auto smallerSet = [] (const SetType& first, const SetType& second)
                              { return (first.size() <= second.size()) ? first: second; };
  auto setX = set_init(SetType(), this->m_data.numVars());
  setX.insert(x);
  auto mbY = this->getMB(y);
  if (mbY.contains(z)) {
    mbY.erase(z);
  }
  auto mbZ = this->getMB(z);
  if (mbZ.contains(y)) {
    mbZ.erase(y);
  }
  const auto& u = smallerSet(mbY, mbZ);
  return !this->m_data.isIndependentAnySubset(y, z, u, setX);
}

template <typename DataType, typename VarType, typename SetType>
/**
 * @brief Adds the neighbors for the given var to the given network.
 *
 * @tparam GraphType The type of graph used for representing the network.
 * @param x The variable for which edges are added for the graph.
 * @param g The causal network.
 * @param directEdges Specifies if the edges should be directed.
 */
template <typename GraphType>
void
ConstraintBasedDiscovery<DataType, VarType, SetType>::addVarNeighbors(
  const VarType x,
  GraphType& g,
  const bool directEdges
) const
{
  auto pcX = this->getPC(x);
  SetType paX;
  for (const auto y: pcX) {
    if (!directEdges) {
      LOG_MESSAGE_IF(x < y, info, "+ Adding the edge %s <-> %s", this->m_data.varName(x), this->m_data.varName(y));
      g.addEdge(x, y);
    }
    else {
      if (set_contains(paX, y) || g.edgeExists(x, y)) {
        // y must have been determined to be a parent or child of x
        continue;
      }
      // Orient the edge
      auto pcY = this->getPC(y);
      // Candidate parents of x, other than y
      auto cpaX = set_difference(pcX, pcY);
      cpaX.erase(y);
      bool isChild = false;
      for (const auto z: cpaX) {
        if (this->isCollider(y, x, z)) {
          isChild = true;
          LOG_MESSAGE(info, "+ Adding the edges %s -> %s <- %s (collider)", this->m_data.varName(y), this->m_data.varName(x), this->m_data.varName(z));
          g.addEdge(y, x);
          g.addEdge(z, x);
          paX.insert(z);
          paX.insert(y);
          break;
        }
      }
      if (!isChild) {
        LOG_MESSAGE(info, "+ Adding the edge %s -> %s", this->m_data.varName(x), this->m_data.varName(y));
        g.addEdge(x, y);
      }
    }
  }
}

template <typename DataType, typename VarType, typename SetType>
/**
 * @brief Top level function for getting the complete causal network.
 *
 * @param directEdges Specifies if the edges of the network should be directed.
 */
BayesianNetwork<VarType>
ConstraintBasedDiscovery<DataType, VarType, SetType>::getNetwork(
  const bool directEdges
) const
{
  auto varNames = this->m_data.varNames(m_allVars);
  BayesianNetwork<VarType> bn(varNames);
  for (const auto x: m_allVars) {
    this->addVarNeighbors(x, bn, directEdges);
  }
  if (directEdges) {
    // First, break any directed cycles in the network
    LOG_MESSAGE_IF(bn.hasDirectedCycles(), info, "* The initial network contains directed cycles");
    while (bn.hasDirectedCycles()) {
      bn.breakDirectedCycles();
    }
    // Then, orient edges by applying Meek's rules
    bool changed = true;
    while (changed) {
      changed = bn.applyMeekRules();
    }
  }
  return bn;
}

#endif // DETAIL_CONSTRAINTBASEDDISCOVERY_HPP_
