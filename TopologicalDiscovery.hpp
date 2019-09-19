/**
 * @file TopologicalDiscovery.hpp
 * @brief Declaration of the TopologicalDiscovery class and all the
 *        classes that are derived from it.
 */
#ifndef TOPOLOGICALDISCOVERY_HPP_
#define TOPOLOGICALDISCOVERY_HPP_

#include "ConstraintBasedDiscovery.hpp"

#include <set>
#include <unordered_map>


/**
 * @brief Abstract base class for causal discovery by first learning PC sets.
 *
 * @tparam DataType Type of the object which is used for querying the data.
 * @tparam VarType Type of variable indices (expected to be an integer type).
 * @tparam SetType Type of set container.
 */
template <typename DataType, typename VarType, typename SetType>
class TopologicalDiscovery : public ConstraintBasedDiscovery<DataType, VarType, SetType> {
public:
  TopologicalDiscovery(const DataType&);

  SetType
  removeFalsePC(const VarType, SetType&) const;

  SetType
  getCandidateMB(const VarType, SetType) const override;

  virtual
  ~TopologicalDiscovery() { }
}; // class TopologicalDiscovery


/**
 * @brief Class that implements Max-Min algorithm for PC discovery,
 *        as described by Tsamardinos et al.
 *
 * @tparam DataType Type of the object which is used for querying the data.
 * @tparam VarType Type of variable indices (expected to be an integer type).
 * @tparam SetType Type of set container.
 */
template <typename DataType, typename VarType, typename SetType>
class MMPC: public TopologicalDiscovery<DataType, VarType, SetType> {
public:
  MMPC(const DataType&);

  SetType
  getCandidatePC(const VarType, SetType) const override;
}; // class MMPC

/**
 * @brief Class that implements HITON algorithm for PC discovery,
 *        as described by Aliferis et al.
 *
 * @tparam DataType Type of the object which is used for querying the data.
 * @tparam VarType Type of variable indices (expected to be an integer type).
 * @tparam SetType Type of set container.
 */
template <typename DataType, typename VarType, typename SetType>
class HITON: public TopologicalDiscovery<DataType, VarType, SetType> {
public:
  HITON(const DataType&);

  SetType
  getCandidatePC(const VarType, SetType) const override;
}; // class HITON

/**
 * @brief Class that implements Semi-interleaved HITON algorithm for PC discovery,
 *        as described by Aliferis et al.
 *
 * @tparam DataType Type of the object which is used for querying the data.
 * @tparam VarType Type of variable indices (expected to be an integer type).
 * @tparam SetType Type of set container.
 */
template <typename DataType, typename VarType, typename SetType>
class SemiInterleavedHITON: public TopologicalDiscovery<DataType, VarType, SetType> {
public:
  SemiInterleavedHITON(const DataType&);

  SetType
  getCandidatePC(const VarType, SetType) const override;
}; // class SemiInterleavedHITON

/**
 * @brief Class that implements GetPC algorithm for PC discovery,
 *        as described by Pena et al.
 *
 * @tparam DataType Type of the object which is used for querying the data.
 * @tparam VarType Type of variable indices (expected to be an integer type).
 * @tparam SetType Type of set container.
 */
template <typename DataType, typename VarType, typename SetType>
class GetPC: public TopologicalDiscovery<DataType, VarType, SetType> {
public:
  GetPC(const DataType&);

  SetType
  getCandidatePC(const VarType, SetType) const override;
}; // class GetPC

#include "detail/TopologicalDiscovery.hpp"

#endif // TOPOLOGICALDISCOVERY_HPP_
