/**
 * @file ConstraintBasedDiscovery.hpp
 * @brief Implementation of the ConstraintBasedDiscovery functions.
 */
#ifndef DETAIL_CONSTRAINTBASEDDISCOVERY_HPP_
#define DETAIL_CONSTRAINTBASEDDISCOVERY_HPP_

#include "../SetUtils.hpp"

#include "utils/Logging.hpp"


template <typename Data, typename Var, typename Set>
/**
 * @brief Constructs the object with the given data.
 *
 * @param data Reference to an object of the Data.
 */
ConstraintBasedDiscovery<Data, Var, Set>::ConstraintBasedDiscovery(
  const mxx::comm& comm,
  const Data& data,
  const Var maxConditioning
) : m_comm(comm),
    m_data(data),
    m_allVars(set_init(Set(), data.numVars())),
    m_maxConditioning(maxConditioning)
{
  for (auto i = 0u; i < data.numVars(); ++i) {
    m_allVars.insert(m_allVars.end(), i);
  }
}

template <typename Data, typename Var, typename Set>
/**
 * @brief Function for getting all the candidates for the given target variable.
 *
 * @param target The index of the target variable.
 *
 * @return The indices of all the variables except the target.
 */
Set
ConstraintBasedDiscovery<Data, Var, Set>::getCandidates(
  const Var target
) const
{
  auto candidates = m_allVars;
  candidates.erase(target);
  return candidates;
}

template <typename Data, typename Var, typename Set>
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
std::pair<Set, bool>
ConstraintBasedDiscovery<Data, Var, Set>::getCandidatePC_cache(
  const Var target,
  Set candidates
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

template <typename Data, typename Var, typename Set>
/**
 * @brief Performs symmetry correction on the given PC set for the given variable.
 *
 * @param target The index of the target variable.
 * @param cpc The set containing the indices of all the variables in
 *            the candidate PC set.
 */
void
ConstraintBasedDiscovery<Data, Var, Set>::symmetryCorrectPC(
  const Var target,
  Set& cpc
) const
{
  auto initial = cpc;
  for (const Var x: initial) {
    auto candidatesX = this->getCandidates(x);
    auto cpcX = this->getCandidatePC_cache(x, std::move(candidatesX)).first;
    if (!set_contains(cpcX, target)) {
      LOG_MESSAGE(info, "- Removing %s from the PC of %s (asymmetry)", this->m_data.varName(x), this->m_data.varName(target));
      cpc.erase(x);
    }
  }
}

template <typename Data, typename Var, typename Set>
/**
 * @brief The top level function for getting the correct PC set
 *        for the given target variable.
 *
 * @param target The index of the target variable.
 *
 * @return A set containing the indices of all the variables
 *         in the correct PC set of the target variable.
 */
Set
ConstraintBasedDiscovery<Data, Var, Set>::getPC(
  const Var target
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

template <typename Data, typename Var, typename Set>
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
std::pair<Set, bool>
ConstraintBasedDiscovery<Data, Var, Set>::getCandidateMB_cache(
  const Var target,
  Set candidates
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

template <typename Data, typename Var, typename Set>
/**
 * @brief Symmetry corrects the candidate MB of the target variable.
 *
 * @param target The index of the target variable.
 * @param cmb The indices of the variables in the candidate MB of the target variable.
 *                The function removes the indices from the candidate set.
 */
void
ConstraintBasedDiscovery<Data, Var, Set>::symmetryCorrectMB(
  const Var target,
  Set& cmb
) const
{
  auto initial = cmb;
  for (const Var x: initial) {
    auto candidatesX = this->getCandidates(x);
    auto cmbX = this->getCandidateMB_cache(x, std::move(candidatesX)).first;
    if (!set_contains(cmbX, target)) {
      LOG_MESSAGE(info, "- Removing %s from the MB of %s (asymmetry)", this->m_data.varName(x), this->m_data.varName(target));
      cmb.erase(x);
    }
  }
}

template <typename Data, typename Var, typename Set>
/**
 * @brief Top level function for getting the MB of the target variable.
 *
 * @param target The index of the target variable.
 */
Set
ConstraintBasedDiscovery<Data, Var, Set>::getMB(
  const Var target
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

template <typename Data, typename Var, typename Set>
/**
 * @brief Function for getting the undirected skeleton network sequentially.
 */
BayesianNetwork<Var>
ConstraintBasedDiscovery<Data, Var, Set>::getSkeleton_sequential(
) const
{
  BayesianNetwork<Var> bn(this->m_data.varNames(m_allVars));
  for (const auto x: m_allVars) {
    auto pcX = this->getPC(x);
    for (const auto y: pcX) {
      if (x < y) {
        LOG_MESSAGE(info, "+ Adding the edge %s <-> %s", this->m_data.varName(x), this->m_data.varName(y));
        bn.addEdge(x, y);
        bn.addEdge(y, x);
      }
    }
  }
  return bn;
}

template <typename Data, typename Var, typename Set>
/**
 * @brief Function for getting the undirected skeleton network in parallel.
 */
BayesianNetwork<Var>
ConstraintBasedDiscovery<Data, Var, Set>::getSkeleton_parallel(
) const
{
  throw std::runtime_error("Getting skeleton in parallel is not implemented for the given algorithm");
}

template <typename Data, typename Var, typename Set>
/**
 * @brief Checks if y-x-z forms a v-structure.
 */
bool
ConstraintBasedDiscovery<Data, Var, Set>::isCollider(
  const Var y,
  const Var x,
  const Var z
) const
{
  static auto smallerSet = [] (const Set& first, const Set& second)
                              { return (first.size() <= second.size()) ? first: second; };
  auto setX = set_init(Set(), this->m_data.numVars());
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
  return !this->m_data.isIndependentAnySubset(y, z, u, setX, this->m_maxConditioning);
}

template <typename Data, typename Var, typename Set>
/**
 * @brief Finds all the potential v-structures in the network.
 */
std::multimap<Var, std::pair<Var, Var>>
ConstraintBasedDiscovery<Data, Var, Set>::findVStructures(
) const
{
  std::multimap<Var, std::pair<Var, Var>> vStructures;
  for (const auto x: m_allVars) {
    auto pcX = this->getPC(x);
    Set paX;
    for (const auto y: pcX) {
      // Candidate parents of x, which are not connected to y
      auto pcY = this->getPC(y);
      auto cpaX = set_difference(pcX, pcY);
      cpaX.erase(y);
      for (const auto z: cpaX) {
        if (!(set_contains(paX, y) && set_contains(paX, z)) && this->isCollider(y, x, z)) {
          LOG_MESSAGE(info, "* Found new v-structure %s -> %s <- %s", this->m_data.varName(y), this->m_data.varName(x), this->m_data.varName(z));
          vStructures.insert(std::make_pair(x, std::make_pair(y, z)));
          paX.insert(y);
          paX.insert(z);
        }
      }
    }
  }
  return vStructures;
}

template <typename Data, typename Var, typename Set>
/**
 * @brief Top level function for getting the complete causal network.
 *
 * @param directEdges Specifies if the edges of the network should be directed.
 * @param isParallel Specifies if the skeleton should be constructed in parallel.
 */
BayesianNetwork<Var>
ConstraintBasedDiscovery<Data, Var, Set>::getNetwork(
  const bool directEdges,
  const bool isParallel
) const
{
  auto bn = isParallel ? this->getSkeleton_parallel() : this->getSkeleton_sequential();
  if (this->m_comm.is_first() && directEdges) {
    // First, orient the v-structures
    auto vStructures = this->findVStructures();
    bn.applyVStructures(vStructures);
    // Then, break any directed cycles in the network
    LOG_MESSAGE_IF(bn.hasDirectedCycles(), info, "* The initial network contains directed cycles");
    while (bn.hasDirectedCycles()) {
      bn.breakDirectedCycles();
    }
    // Finally, orient edges by applying Meek's rules
    bool changed = true;
    while (changed) {
      changed = bn.applyMeekRules();
    }
  }
  return bn;
}

#endif // DETAIL_CONSTRAINTBASEDDISCOVERY_HPP_
