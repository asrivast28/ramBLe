/**
 * @file StdSetUtils.hpp
 * @brief Implementation of functions for operations on std::set.
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
#ifndef DETAIL_STDSETUTILS_HPP_
#define DETAIL_STDSETUTILS_HPP_

#include <algorithm>
#include <set>


// Definition of all the operations on std::set
#define DEFINE_STD_SET_OPERATIONS(type) \
template <> \
std::set<type> \
set_init<std::set<type>, type>( \
  std::set<type>&& set, \
  const type \
) \
{ \
  return set; \
} \
 \
template <> \
bool \
set_contains<std::set<type>, type>( \
  const std::set<type>& set, \
  const type value \
) \
{ \
  return (set.find(value) != set.end()); \
} \
 \
template <> \
std::set<type> \
set_union<std::set<type>>( \
  const std::set<type>& first, \
  const std::set<type>& second \
) \
{ \
  std::set<type> result; \
  std::set_union(first.begin(), first.end(), second.begin(), second.end(), std::inserter(result, result.begin())); \
  return result; \
} \
 \
template <> \
std::set<type> \
set_difference<std::set<type>>( \
  const std::set<type>& first, \
  const std::set<type>& second \
) \
{ \
  std::set<type> result; \
  std::set_difference(first.begin(), first.end(), second.begin(), second.end(), std::inserter(result, result.begin())); \
  return result; \
}

DEFINE_STD_SET_OPERATIONS(uint8_t)
DEFINE_STD_SET_OPERATIONS(uint16_t)

/**
 * @brief Function for getting the output represention of a set.
 */
template <typename Element>
std::ostream&
operator<<(
  std::ostream& stream,
  const std::set<Element>& set
)
{
  stream << "{";
  for (auto elem : set) {
    stream << static_cast<uint32_t>(elem) << ",";
  }
  if (set.size() > 0) {
    stream << '\b';
  }
  stream << "}";
  return stream;
}

#endif // DETAIL_STDSETUTILS_HPP_
