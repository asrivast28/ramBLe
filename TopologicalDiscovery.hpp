/**
 * @file TopologicalDiscovery.hpp
 * @brief Declaration of the TopologicalDiscovery class and all the
 *        classes that are derived from it.
 */
#ifndef TOPOLOGICALDISCOVERY_HPP_
#define TOPOLOGICALDISCOVERY_HPP_

#include "MBDiscovery.hpp"

#include <set>
#include <unordered_map>


/**
 * @brief Abstract base class for topological discovery of MBs.
 *
 * @tparam DataType Type of the object which is used for querying the data.
 * @tparam VarType Type of variable indices (expected to be an integer type).
 */
template <typename DataType, typename VarType>
class TopologicalDiscovery : public MBDiscovery<DataType, VarType> {
public:
  TopologicalDiscovery(const DataType&);

  std::set<VarType>
  removeFalsePC(const VarType, std::set<VarType>&) const;

  std::set<VarType>
  getCandidateMB(const VarType, std::set<VarType>) const override;

  virtual
  ~TopologicalDiscovery() { }
}; // class TopologicalDiscovery


/**
 * @brief Class that implements Max-Min algorithm for PC discovery,
 *        as described by Tsamardinos et al.
 *
 * @tparam DataType Type of the object which is used for querying the data.
 * @tparam VarType Type of variable indices (expected to be an integer type).
 */
template <typename DataType, typename VarType>
class MMPC: public TopologicalDiscovery<DataType, VarType> {
public:
  MMPC(const DataType&);

  std::set<VarType>
  getCandidatePC(const VarType, std::set<VarType>) const override;
}; // class MMPC

/**
 * @brief Class that implements HITON algorithm for PC discovery,
 *        as described by Aliferis et al.
 *
 * @tparam DataType Type of the object which is used for querying the data.
 * @tparam VarType Type of variable indices (expected to be an integer type).
 */
template <typename DataType, typename VarType>
class HITON: public TopologicalDiscovery<DataType, VarType> {
public:
  HITON(const DataType&);

  std::set<VarType>
  getCandidatePC(const VarType, std::set<VarType>) const override;
}; // class HITON

/**
 * @brief Class that implements Semi-interleaved HITON algorithm for PC discovery,
 *        as described by Aliferis et al.
 *
 * @tparam DataType Type of the object which is used for querying the data.
 * @tparam VarType Type of variable indices (expected to be an integer type).
 */
template <typename DataType, typename VarType>
class SemiInterleavedHITON: public TopologicalDiscovery<DataType, VarType> {
public:
  SemiInterleavedHITON(const DataType&);

  std::set<VarType>
  getCandidatePC(const VarType, std::set<VarType>) const override;
}; // class SemiInterleavedHITON

/**
 * @brief Class that implements GetPC algorithm for PC discovery,
 *        as described by Pena et al.
 *
 * @tparam DataType Type of the object which is used for querying the data.
 * @tparam VarType Type of variable indices (expected to be an integer type).
 */
template <typename DataType, typename VarType>
class GetPC: public TopologicalDiscovery<DataType, VarType> {
public:
  GetPC(const DataType&);

  std::set<VarType>
  getCandidatePC(const VarType, std::set<VarType>) const override;
}; // class GetPC

#include "detail/TopologicalDiscovery.hpp"

#endif // TOPOLOGICALDISCOVERY_HPP_
