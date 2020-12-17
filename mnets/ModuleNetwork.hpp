/**
 * @file ModuleNetwork.hpp
 * @brief Declaration of the module network class.
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
#ifndef MODULENETWORK_HPP_
#define MODULENETWORK_HPP_

#include "BayesianNetwork.hpp"


template <typename Var>
class ModuleNetwork : public BayesianNetwork<Var> {
public:
  ModuleNetwork(const std::vector<std::string>&);

  ~ModuleNetwork() { }
}; // class ModuleNetwork

#include "detail/ModuleNetwork.hpp"

#endif // MODULENETWORK_HPP_
