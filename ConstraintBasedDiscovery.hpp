/**
 * @file ConstraintBasedDiscovery.hpp
 * @brief Declaration of the ConstraintBasedDiscovery class.
 */
#ifndef CONSTRAINTBASEDDISCOVERY_HPP_
#define CONSTRAINTBASEDDISCOVERY_HPP_

#include "BayesianNetwork.hpp"

#include "mxx/comm.hpp"

#include <map>
#include <unordered_map>


/**
 * @brief Abstract base class for causal discovery using constraint-based learning.
 *
 * @tparam Data Type of the object which is used for querying the data.
 * @tparam Var Type of variable indices (expected to be an integer type).
 * @tparam Set Type of set container.
 */
template <typename Data, typename Var, typename Set>
class ConstraintBasedDiscovery {
public:
  ConstraintBasedDiscovery(const mxx::comm&, const Data&, const Var);

  Set
  getPC(const Var) const;

  Set
  getMB(const Var) const;

  BayesianNetwork<Var>
  getNetwork(const bool = false) const;

  virtual
  ~ConstraintBasedDiscovery() { }

protected:
  Set
  getCandidates(const Var) const;

  virtual
  Set
  getCandidatePC(const Var, Set) const = 0;

  virtual
  Set
  getCandidateMB(const Var, Set) const = 0;

  virtual
  BayesianNetwork<Var>
  getSkeleton() const;

private:
  std::pair<Set, bool>
  getCandidatePC_cache(const Var, Set) const;

  void
  symmetryCorrectPC(const Var, Set&) const;

  std::pair<Set, bool>
  getCandidateMB_cache(const Var, Set) const;

  void
  symmetryCorrectMB(const Var, Set&) const;

  bool
  isCollider(const Var, const Var, const Var) const;

  std::multimap<Var, std::pair<Var, Var>>
  findVStructures() const;

protected:
  const mxx::comm& m_comm;
  const Data m_data;
  Set m_allVars;
  const Var m_maxConditioning;

private:
  mutable std::unordered_map<Var, std::pair<Set, bool>> m_cachedPC;
  mutable std::unordered_map<Var, std::pair<Set, bool>> m_cachedMB;
}; // class ConstraintBasedDiscovery

#include "detail/ConstraintBasedDiscovery.hpp"

#endif // CONSTRAINTBASEDDISCOVERY_HPP_
