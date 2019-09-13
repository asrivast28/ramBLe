/**
 * @file MBDiscovery.hpp
 * @brief Declaration of the MBDiscovery class.
 */
#ifndef MBDISCOVERY_HPP_
#define MBDISCOVERY_HPP_

#include <set>
#include <unordered_map>


/**
 * @brief Abstract base class for discovering Markov Blankets (MB).
 *
 * @tparam DataType Type of the object which is used for querying the data.
 * @tparam VarType Type of variable indices (expected to be an integer type).
 * @tparam SetType Type of set container.
 */
template <typename DataType, typename VarType, typename SetType>
class MBDiscovery {
public:
  MBDiscovery(const DataType&);

  SetType
  getPC(const VarType) const;

  SetType
  getMB(const VarType) const;

  virtual
  ~MBDiscovery() { };

protected:
  SetType
  getCandidates(const VarType) const;

  virtual
  SetType
  getCandidatePC(const VarType, SetType) const = 0;

  virtual
  SetType
  getCandidateMB(const VarType, SetType) const = 0;

private:
  SetType
  getCandidatePC_cache(const VarType, SetType) const;

  void
  symmetryCorrectPC(const VarType, SetType&) const;

  SetType
  getCandidateMB_cache(const VarType, SetType) const;

  void
  symmetryCorrectMB(const VarType, SetType&) const;

protected:
  const DataType m_data;
  SetType m_vars;

private:
  mutable std::unordered_map<VarType, SetType> m_cachedPC;
  mutable std::unordered_map<VarType, SetType> m_cachedMB;
}; // class MBDiscovery

#include "detail/MBDiscovery.hpp"

#endif // MBDISCOVERY_HPP_
