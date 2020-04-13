/**
 * @file CTCounter.hpp
 * @brief Declaration of the functions used for counting using contingency table.
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
#ifndef CT_COUNTER_HPP
#define CT_COUNTER_HPP

#include <algorithm>
#include <vector>


/**
 * @brief Class that provides functionality for counting using contingency tables.
 *
 * @tparam Data The type of data to be stored.
 */
template <typename Data = uint8_t>
class CTCounter {
public:
  using data_type = Data;

public:
  // declaration of the static method to create the counter
  template <typename Iter>
  static
  CTCounter<Data> create(const uint32_t, const uint32_t, Iter);

public:
  uint32_t n() const { return m_nvars; }

  uint32_t m() const { return m_nobs; }

  uint32_t r(const uint32_t xi) const { return m_arity[xi]; }

  const std::vector<data_type>& data() const { return m_data; }

private:
  std::vector<data_type> m_data;
  std::vector<data_type> m_arity; // we do not expect more than 255 states
  uint32_t m_nvars;
  uint32_t m_nobs;
}; // class CTCounter


template <typename Data>
template <typename Iter>
CTCounter<Data>
CTCounter<Data>::create(const uint32_t n, const uint32_t m, Iter it) {
  CTCounter<data_type> ct;
  ct.m_nvars = n;
  ct.m_nobs = m;
  ct.m_data.resize(n * m);
  ct.m_arity.resize(n, 0);
  for (uint32_t i = 0; i < n * m; ++i, ++it) {
    ct.m_data[i] = *it;
  }

  for (uint32_t xi = 0; xi < n; ++xi) {
    auto result = std::minmax_element(ct.m_data.begin() + xi * m, ct.m_data.begin() + (xi + 1) * m);
    auto min = *(result.first);
    ct.m_arity[xi] = (*result.second + 1) - min;
    std::for_each(ct.m_data.begin() + xi * m, ct.m_data.begin() + (xi + 1) * m, [min](data_type& x) { return x -= min; } );
  }

  return ct;
} // CTCounter<Data>::create

#endif // CT_COUNTER_HPP
