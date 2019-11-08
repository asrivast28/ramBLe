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
  DirectDiscovery(const Data&);

  Set
  getCandidatePC(const Var, Set) const override;

  Set
  shrinkMB(const Var, Set&) const;

  virtual
  ~DirectDiscovery() { }

private:
  bool
  evaluateCandidatePC(const Var, const Var, const Set&, const Set&) const;
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
  GSMB(const Data&);

  Set
  getCandidateMB(const Var, Set) const override;
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
  IAMB(const Data&);

  Set
  getCandidateMB(const Var, Set) const override;
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
  InterIAMB(const Data&);

  Set
  getCandidateMB(const Var, Set) const override;
}; // class InterIAMB

#include "detail/DirectDiscovery.hpp"

#endif // DIRECTDISCOVERY_HPP_
