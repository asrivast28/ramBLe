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
 */
template <typename DataType, typename VarType>
class MBDiscovery {
public:
  MBDiscovery(const DataType&);

  std::set<VarType>
  getPC(const VarType) const;

  std::set<VarType>
  getMB(const VarType) const;

  virtual
  ~MBDiscovery() { };

protected:
  std::set<VarType>
  getCandidates(const VarType) const;

  virtual
  std::set<VarType>
  getCandidatePC(const VarType, std::set<VarType>) const = 0;

  virtual
  std::set<VarType>
  getCandidateMB(const VarType, std::set<VarType>) const = 0;

private:
  std::set<VarType>
  getCandidatePC_cache(const VarType, std::set<VarType>) const;

  void
  symmetryCorrectPC(const VarType, std::set<VarType>&) const;

  std::set<VarType>
  getCandidateMB_cache(const VarType, std::set<VarType>) const;

  void
  symmetryCorrectMB(const VarType, std::set<VarType>&) const;

protected:
  const DataType m_data;
  std::set<VarType> m_vars;

private:
  mutable std::unordered_map<VarType, std::set<VarType>> m_cachedPC;
  mutable std::unordered_map<VarType, std::set<VarType>> m_cachedMB;
}; // class MBDiscovery

#include "detail/MBDiscovery.hpp"

#endif // MBDISCOVERY_HPP_
