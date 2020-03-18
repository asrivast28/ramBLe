/**
 * @file DirectDiscovery.hpp
 * @brief Declaration of the DirectDiscovery class and all the
 *        classes that are derived from it.
 */
#ifndef DIRECTDISCOVERY_HPP_
#define DIRECTDISCOVERY_HPP_

#include "ConstraintBasedDiscovery.hpp"


/**
 * @brief Abstract base class for causal discovery by first learning MB sets.
 *
 * @tparam Data Type of the object which is used for querying the data.
 * @tparam Var Type of variable indices (expected to be an integer type).
 * @tparam Set Type of set container.
 */
template <typename Data, typename Var, typename Set>
class DirectDiscovery : public ConstraintBasedDiscovery<Data, Var, Set> {
public:
  DirectDiscovery(const mxx::comm&, const Data&, const Var);

  virtual
  ~DirectDiscovery() { }

protected:
  Set
  shrinkMB(const Var, Set&) const;

  virtual
  std::pair<Var, double>
  pickBestCandidate(const Var, const Set&, const Set&) const;

  virtual
  void
  updatePValues(std::vector<std::tuple<Var, Var, double>>&, const std::unordered_map<Var, Set>&) const;

  std::set<std::tuple<Var, Var, double>>
  growAll(const std::vector<std::tuple<Var, Var, double>>&, std::unordered_map<Var, Set>&) const;

  std::set<std::pair<Var, Var>>
  shrinkAll(std::unordered_map<Var, Set>&) const;

  void
  syncBlankets(std::unordered_map<Var, Set>&) const;

  bool
  fixImbalance(std::vector<std::tuple<Var, Var, double>>&, const double) const;

  void
  syncMissingBlankets(const std::vector<std::tuple<Var, Var, double>>&, std::unordered_map<Var, Set>&) const;

  virtual
  void
  growShrink(std::vector<std::tuple<Var, Var, double>>&&, std::unordered_map<Var, Set>&, std::set<std::pair<Var, Var>>&, const double) const;

private:
  Set
  getCandidatePC(const Var, Set&&) const override;

  bool
  evaluateCandidatePC(const Var, const Var, const Set&, const Set&) const;

  std::vector<std::pair<Var, Var>>
  symmetryCorrect(const std::unordered_map<Var, Set>&&, const std::set<std::pair<Var, Var>>&&) const;

  BayesianNetwork<Var>
  getSkeleton_parallel(const double) const override;
}; // class DirectDiscovery

/**
 * @brief Class that implements Grow-Shrink strategy for MB discovery,
 *        as described by Margaritis & Thrun.
 *
 * @tparam Data Type of the object which is used for querying the data.
 * @tparam Var Type of variable indices (expected to be an integer type).
 * @tparam Set Type of set container.
 */
template <typename Data, typename Var, typename Set>
class GSMB: public DirectDiscovery<Data, Var, Set> {
public:
  GSMB(const mxx::comm&, const Data&, const Var = std::numeric_limits<Var>::max());

private:
  Set
  getCandidateMB(const Var, Set&&) const override;

  void
  updatePValues(std::vector<std::tuple<Var, Var, double>>&, const std::unordered_map<Var, Set>&) const override;

  std::pair<Var, double>
  pickBestCandidate(const Var, const Set&, const Set&) const override;
}; // class GSMB

/**
 * @brief Class that implements Incremental Association strategy for MB discovery,
 *        as described by Tsamardinos et al.
 *
 * @tparam Data Type of the object which is used for querying the data.
 * @tparam Var Type of variable indices (expected to be an integer type).
 * @tparam Set Type of set container.
 */
template <typename Data, typename Var, typename Set>
class IAMB: public DirectDiscovery<Data, Var, Set> {
public:
  IAMB(const mxx::comm&, const Data&, const Var = std::numeric_limits<Var>::max());

private:
  Set
  getCandidateMB(const Var, Set&&) const override;
}; // class IAMB

/**
 * @brief Class that implements Interleaved Incremental Association strategy
 *        for MB discovery, as described by Tsamardinos et al.
 *
 * @tparam Data Type of the object which is used for querying the data.
 * @tparam Var Type of variable indices (expected to be an integer type).
 * @tparam Set Type of set container.
 */
template <typename Data, typename Var, typename Set>
class InterIAMB: public DirectDiscovery<Data, Var, Set> {
public:
  InterIAMB(const mxx::comm&, const Data&, const Var = std::numeric_limits<Var>::max());

private:
  Set
  getCandidateMB(const Var, Set&&) const override;

  void
  growShrink(std::vector<std::tuple<Var, Var, double>>&&, std::unordered_map<Var, Set>&, std::set<std::pair<Var, Var>>&, const double) const override;
}; // class InterIAMB

#include "detail/DirectDiscovery.hpp"

#endif // DIRECTDISCOVERY_HPP_
