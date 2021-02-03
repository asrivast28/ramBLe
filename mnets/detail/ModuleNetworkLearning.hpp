/**
 * @file ModuleNetworkLearning.hpp
 * @brief Implementation of the classes for constraint-based learning.
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
#ifndef DETAIL_MODULENETWORKLEARNING_HPP_
#define DETAIL_MODULENETWORKLEARNING_HPP_

#include "SetUtils.hpp"


NotImplementedError::NotImplementedError(
  const std::string& message
) : std::logic_error(message)
{
}

template <typename Data, typename Var, typename Set>
/**
 * @brief Constructs the object with the given data.
 *
 * @param data Reference to an object of the Data.
 */
ModuleNetworkLearning<Data, Var, Set>::ModuleNetworkLearning(
  const mxx::comm& comm,
  const Data& data
) : m_comm(comm),
    m_data(data),
    m_allVars(set_init(Set(), data.numVars()))
{
  for (auto i = 0u; i < data.numVars(); ++i) {
    m_allVars.insert(m_allVars.end(), i);
  }
}

template <typename Data, typename Var, typename Set>
/**
 * @brief Top level function for getting the module network.
 *
 * @param isParallel Specifies if the network should be learned in parallel.
 */
void
ModuleNetworkLearning<Data, Var, Set>::learnNetwork(
  const bool isParallel,
  const pt::ptree& algoConfigs,
  const std::string& outputDir
) const
{
  isParallel ? this->learnNetwork_parallel(algoConfigs, outputDir) :
               this->learnNetwork_sequential(algoConfigs, outputDir);
}

#endif // DETAIL_MODULENETWORKLEARNING_HPP_
