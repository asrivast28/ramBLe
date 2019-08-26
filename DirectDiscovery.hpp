/**
 * @file DirectDiscovery.hpp
 * @brief Declaration of the DirectDiscovery class and all the
 *        classes that are derived from it.
 */
#ifndef DIRECTDISCOVERY_HPP_
#define DIRECTDISCOVERY_HPP_

#include "MBDiscovery.hpp"

#include <set>


/**
 * @brief Abstract base class for direct discovery of MBs.
 *
 * @tparam DataType Type of the object which is used for querying the data.
 * @tparam VarType Type of variable indices (expected to be an integer type).
 */
template <typename DataType, typename VarType>
class DirectDiscovery : public MBDiscovery<DataType, VarType> {
public:
  DirectDiscovery(const DataType&);

  std::set<VarType>
  shrinkMB(const VarType, std::set<VarType>&) const;

  virtual
  ~DirectDiscovery() { }
}; // class DirectDiscovery

/**
 * @brief Class that implements Grow-Shrink strategy for MB discovery,
 *         as described by Margaritis & Thrun.
 *
 * @tparam DataType Type of the object which is used for querying the data.
 * @tparam VarType Type of variable indices (expected to be an integer type).
 */
template <typename DataType, typename VarType>
class GSMB: public DirectDiscovery<DataType, VarType> {
public:
  GSMB(const DataType&);

  std::set<VarType>
  getCandidateMB(const VarType, std::set<VarType>) const override;
}; // class GSMB

/**
 * @brief Class that implements Incremental Association strategy for MB discovery,
 *         as described by Tsamardinos et al.
 *
 * @tparam DataType Type of the object which is used for querying the data.
 * @tparam VarType Type of variable indices (expected to be an integer type).
 */
template <typename DataType, typename VarType>
class IAMB: public DirectDiscovery<DataType, VarType> {
public:
  IAMB(const DataType&);

  std::set<VarType>
  getCandidateMB(const VarType, std::set<VarType>) const override;
}; // class IAMB

/**
 * @brief Class that implements Interleaved Incremental Association strategy
 *         for MB discovery, as described by Tsamardinos et al.
 *
 * @tparam DataType Type of the object which is used for querying the data.
 * @tparam VarType Type of variable indices (expected to be an integer type).
 */
template <typename DataType, typename VarType>
class InterIAMB: public DirectDiscovery<DataType, VarType> {
public:
  InterIAMB(const DataType&);

  std::set<VarType>
  getCandidateMB(const VarType, std::set<VarType>) const override;
}; // class InterIAMB

#include "detail/DirectDiscovery.hpp"

#endif // DIRECTDISCOVERY_HPP_
