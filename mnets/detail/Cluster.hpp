/**
 * @file Cluster.hpp
 * @brief Implementation of functionality for storing clusters.
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
#ifndef DETAIL_CLUSTER_HPP_
#define DETAIL_CLUSTER_HPP_

#include "SetUtils.hpp"

#include "utils/Logging.hpp"

#include <random>


/**
 * @brief Base class that provides functionality for cluster storage.
 *
 * @tparam Data Type of the data provider.
 * @tparam Var Type of variables stored in the cluster.
 * @tparam Set Type of container used to store the clusters.
 */
template <typename Data, typename Var, typename Set>
class Cluster {
public:
  Cluster(
    const Data& data,
    std::mt19937* const generator,
    const Var max
  ) : m_elements(max),
      m_data(data),
      m_generator(generator)
  {
  }

  Cluster(
    const Data& data,
    std::mt19937* const generator,
    const Set& elements
  ) : m_elements(elements),
      m_data(data),
      m_generator(generator)
  {
  }

  Cluster(const Cluster& other)
    : m_elements(other.m_elements),
      m_data(other.m_data),
      m_generator(other.m_generator)
  {
  }

  Cluster(
    const Cluster& first,
    const Cluster& second
  ) : m_elements(set_union(first.m_elements, second.m_elements)),
      m_data(first.m_data),
      m_generator(first.m_generator)
  {
  }

  void
  merge(const Cluster& other)
  {
    // Merge the elements from the other cluster into this cluster
    m_elements = set_union(m_elements, other.m_elements);
  }

  void
  insert(const Var e)
  {
    m_elements.insert(e);
  }

  void
  erase(const Var e)
  {
    m_elements.erase(e);
  }

  void
  clear()
  {
    m_elements.clear();
  }

  bool
  empty() const
  {
    return m_elements.empty();
  }

  Var
  max() const
  {
    return m_elements.max();
  }

  Var
  size() const
  {
    return m_elements.size();
  }

  const Set&
  elements() const
  {
    return m_elements;
  }

  ~Cluster()
  {
  }

public:
  template <typename D, typename V, typename S>
  friend
  std::ostream&
  operator<<(std::ostream&, const Cluster<D, V, S>&);

protected:
  Set m_elements;
  const Data& m_data;
  std::mt19937* const m_generator;
}; // class Cluster

/**
 * @brief Function for getting the output represention of a cluster.
 */
template <typename Data, typename Var, typename Set>
std::ostream&
operator<<(
  std::ostream& stream,
  const Cluster<Data, Var, Set>& cluster
)
{
  stream << cluster.m_elements;
  return stream;
}

#endif // DETAIL_CLUSTER_HPP_
