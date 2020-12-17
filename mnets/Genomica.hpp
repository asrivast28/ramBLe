/**
 * @file Genomica.hpp
 * @brief Declaration of the class for learning module networks
 *        using the approach of Genomica.
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
#ifndef GENOMICA_HPP_
#define GENOMICA_HPP_

#include "ModuleNetworkLearning.hpp"


/**
 * @brief Class that implements learning of module networks
 *        using the approach of Genomica.
 *
 * @tparam Data Type of the object which is used for querying the data.
 * @tparam Var Type of variable indices (expected to be an integer type).
 * @tparam Set Type of set container.
 */
template <typename Data, typename Var, typename Set>
class Genomica : public ModuleNetworkLearning<Data, Var, Set> {
public:
  Genomica(const mxx::comm&, const Data&);

  ~Genomica();

protected:
  ModuleNetwork<Var>
  getNetwork_sequential(const pt::ptree&) const;

  ModuleNetwork<Var>
  getNetwork_parallel(const pt::ptree&) const;
}; // class Genomica

#include "detail/Genomica.hpp"

#endif // GENOMICA_HPP_
