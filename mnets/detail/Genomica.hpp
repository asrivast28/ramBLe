/**
 * @file Genomica.hpp
 * @brief Implementation of the class for learning module networks
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
#ifndef DETAIL_GENOMICA_HPP_
#define DETAIL_GENOMICA_HPP_

#include "utils/Logging.hpp"


template <typename Data, typename Var, typename Set>
Genomica<Data, Var, Set>::Genomica(
  const mxx::comm& comm,
  const Data& data
) : ModuleNetworkLearning<Data, Var, Set>(comm, data)
{
}

template <typename Data, typename Var, typename Set>
/**
 * @brief Default destructor.
 */
Genomica<Data, Var, Set>::~Genomica(
)
{
}

template <typename Data, typename Var, typename Set>
void
Genomica<Data, Var, Set>::learnNetwork_sequential(
  const pt::ptree&,
  const std::string&
) const
{
  throw NotImplementedError("Genomica: Sequential algorithm is not implemented yet");
}

template <typename Data, typename Var, typename Set>
void
Genomica<Data, Var, Set>::learnNetwork_parallel(
  const pt::ptree&,
  const std::string&
) const
{
  throw NotImplementedError("Genomica: Parallel algorithm is not implemented yet");
}

#endif // DETAIL_GENOMICA_HPP_
