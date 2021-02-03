/**
 * @file ModuleNetworkLearning.hpp
 * @brief Declaration of the classes for learning module networks.
 * @author Ankit Srivastava <asrivast@gatech.edu>
 *
 * Copyright 2020 Georgia Institute of Technology
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
#ifndef MODULENETWORKLEARNING_HPP_
#define MODULENETWORKLEARNING_HPP_

#include "ModuleNetwork.hpp"

#include "mxx/comm.hpp"

#include <boost/property_tree/ptree.hpp>

#include <stdexcept>

namespace pt = boost::property_tree;

/**
 * @brief Exception to be used for not implemented functions.
 */
class NotImplementedError : public std::logic_error {
public:
  NotImplementedError(const std::string&);
}; // class NotImplementedError

/**
 * @brief Abstract base class for learning module networks.
 *        All the module network learning algorithm implementations should be a
 *        descendant of this class.
 *
 * @tparam Data Type of the object which is used for querying the data.
 * @tparam Var Type of variable indices (expected to be an integer type).
 * @tparam Set Type of set container.
 */
template <typename Data, typename Var, typename Set>
class ModuleNetworkLearning {
public:
  ModuleNetworkLearning(const mxx::comm&, const Data&);

  virtual
  void
  learnNetwork(const bool, const pt::ptree&, const std::string&) const;

  virtual
  ~ModuleNetworkLearning() { }

protected:
  virtual
  void
  learnNetwork_sequential(const pt::ptree&, const std::string&) const = 0;

  virtual
  void
  learnNetwork_parallel(const pt::ptree&, const std::string&) const = 0;

protected:
  const mxx::comm& m_comm;
  const Data m_data;
  Set m_allVars;
}; // class ModuleNetworkLearning

#include "detail/ModuleNetworkLearning.hpp"

#endif // MODULENETWORKLEARNING_HPP_
