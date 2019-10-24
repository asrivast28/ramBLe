/**
 * @file TopologicalDiscovery.hpp
 * @brief Declaration of the TopologicalDiscovery class and all the
 *        classes that are derived from it.
 */
#ifndef TOPOLOGICALDISCOVERY_HPP_
#define TOPOLOGICALDISCOVERY_HPP_

#include "ConstraintBasedDiscovery.hpp"


/**
 * @brief Abstract base class for causal discovery by first learning PC sets.
 *
 * @tparam Data Type of the object which is used for querying the data.
 * @tparam Var Type of variable indices (expected to be an integer type).
 * @tparam Set Type of set container.
 */
template <typename Data, typename Var, typename Set>
class TopologicalDiscovery : public ConstraintBasedDiscovery<Data, Var, Set> {
public:
  TopologicalDiscovery(const Data&);

  Set
  removeFalsePC(const Var, Set&) const;

  Set
  getCandidateMB(const Var, Set) const override;

  virtual
  ~TopologicalDiscovery() { }
}; // class TopologicalDiscovery


/**
 * @brief Class that implements Max-Min algorithm for PC discovery,
 *        as described by Tsamardinos et al.
 *
 * @tparam Data Type of the object which is used for querying the data.
 * @tparam Var Type of variable indices (expected to be an integer type).
 * @tparam Set Type of set container.
 */
template <typename Data, typename Var, typename Set>
class MMPC: public TopologicalDiscovery<Data, Var, Set> {
public:
  MMPC(const Data&);

  Set
  getCandidatePC(const Var, Set) const override;
}; // class MMPC

/**
 * @brief Class that implements HITON algorithm for PC discovery,
 *        as described by Aliferis et al.
 *
 * @tparam Data Type of the object which is used for querying the data.
 * @tparam Var Type of variable indices (expected to be an integer type).
 * @tparam Set Type of set container.
 */
template <typename Data, typename Var, typename Set>
class HITON: public TopologicalDiscovery<Data, Var, Set> {
public:
  HITON(const Data&);

  Set
  getCandidatePC(const Var, Set) const override;
}; // class HITON

/**
 * @brief Class that implements Semi-interleaved HITON algorithm for PC discovery,
 *        as described by Aliferis et al.
 *
 * @tparam Data Type of the object which is used for querying the data.
 * @tparam Var Type of variable indices (expected to be an integer type).
 * @tparam Set Type of set container.
 */
template <typename Data, typename Var, typename Set>
class SemiInterleavedHITON: public TopologicalDiscovery<Data, Var, Set> {
public:
  SemiInterleavedHITON(const Data&);

  Set
  getCandidatePC(const Var, Set) const override;
}; // class SemiInterleavedHITON

/**
 * @brief Class that implements GetPC algorithm for PC discovery,
 *        as described by Pena et al.
 *
 * @tparam Data Type of the object which is used for querying the data.
 * @tparam Var Type of variable indices (expected to be an integer type).
 * @tparam Set Type of set container.
 */
template <typename Data, typename Var, typename Set>
class GetPC: public TopologicalDiscovery<Data, Var, Set> {
public:
  GetPC(const Data&);

  Set
  getCandidatePC(const Var, Set) const override;
}; // class GetPC

#include "detail/TopologicalDiscovery.hpp"

#endif // TOPOLOGICALDISCOVERY_HPP_
