/**
 * @file DirectDiscovery.hpp
 * @brief Declaration of the DirectDiscovery class and all the
 *        classes that are derived from it.
 */
#ifndef DIRECTDISCOVERY_HPP_
#define DIRECTDISCOVERY_HPP_

#include "ConstraintBasedDiscovery.hpp"

#include <set>


/**
 * @brief Abstract base class for causal discovery by first learning MB sets.
 *
 * @tparam DataType Type of the object which is used for querying the data.
 * @tparam VarType Type of variable indices (expected to be an integer type).
 * @tparam SetType Type of set container.
 */
template <typename DataType, typename VarType, typename SetType>
class DirectDiscovery : public ConstraintBasedDiscovery<DataType, VarType, SetType> {
public:
  DirectDiscovery(const DataType&);

  SetType
  getCandidatePC(const VarType, SetType) const override;

  SetType
  shrinkMB(const VarType, SetType&) const;

  virtual
  ~DirectDiscovery() { }
}; // class DirectDiscovery

/**
 * @brief Class that implements Grow-Shrink strategy for MB discovery,
 *        as described by Margaritis & Thrun.
 *
 * @tparam DataType Type of the object which is used for querying the data.
 * @tparam VarType Type of variable indices (expected to be an integer type).
 * @tparam SetType Type of set container.
 */
template <typename DataType, typename VarType, typename SetType>
class GSMB: public DirectDiscovery<DataType, VarType, SetType> {
public:
  GSMB(const DataType&);

  SetType
  getCandidateMB(const VarType, SetType) const override;
}; // class GSMB

/**
 * @brief Class that implements Incremental Association strategy for MB discovery,
 *        as described by Tsamardinos et al.
 *
 * @tparam DataType Type of the object which is used for querying the data.
 * @tparam VarType Type of variable indices (expected to be an integer type).
 * @tparam SetType Type of set container.
 */
template <typename DataType, typename VarType, typename SetType>
class IAMB: public DirectDiscovery<DataType, VarType, SetType> {
public:
  IAMB(const DataType&);

  SetType
  getCandidateMB(const VarType, SetType) const override;
}; // class IAMB

/**
 * @brief Class that implements Interleaved Incremental Association strategy
 *        for MB discovery, as described by Tsamardinos et al.
 *
 * @tparam DataType Type of the object which is used for querying the data.
 * @tparam VarType Type of variable indices (expected to be an integer type).
 * @tparam SetType Type of set container.
 */
template <typename DataType, typename VarType, typename SetType>
class InterIAMB: public DirectDiscovery<DataType, VarType, SetType> {
public:
  InterIAMB(const DataType&);

  SetType
  getCandidateMB(const VarType, SetType) const override;
}; // class InterIAMB

#include "detail/DirectDiscovery.hpp"

#endif // DIRECTDISCOVERY_HPP_
