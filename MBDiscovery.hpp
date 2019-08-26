/**
 * @file MBDiscovery.hpp
 * @brief Declaration of the MBDiscovery class.
 */
#ifndef MBDISCOVERY_HPP_
#define MBDISCOVERY_HPP_

#include <set>


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
  getCandidates(const VarType) const;

  virtual
  std::set<VarType>
  getCandidateMB(const VarType, std::set<VarType>) const = 0;

  void
  symmetryCorrectMB(const VarType, std::set<VarType>&) const;

  std::set<VarType>
  getMB(const VarType) const;

  virtual
  ~MBDiscovery() { };

protected:
  const DataType m_data;
  std::set<VarType> m_vars;
}; // class MBDiscovery

#include "detail/MBDiscovery.hpp"

#endif // MBDISCOVERY_HPP_
