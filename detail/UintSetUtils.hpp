/**
 * @file UintSetUtils.hpp
 * @brief Implementation of functions for operations on UintSet.
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
#ifndef DETAIL_UINTSETUTILS_HPP_
#define DETAIL_UINTSETUTILS_HPP_

#include "../UintSet.hpp"

#include "mxx/reduction.hpp"

#include <algorithm>
#include <type_traits>

// MVAPICH2 (pre release-2.3.3)'s MPI_Allreduce segfaults in the function MPIR_Reduce_scatter_ring_2lvl.
// Therefore, we use bounded-size MPI_Allreduce so that the function is never called.
#if defined MVAPICH2_NUMVERSION && MVAPICH2_NUMVERSION < 20303300
#define ALLREDUCE_MAXSIZE 65535u
#else
#define ALLREDUCE_MAXSIZE std::numeric_limits<uint32_t>::max()
#endif


/**
 * @brief Class that iterates over all the subsets of the given size of a given UintSet.
 *
 * @tparam Element Type of the variable (expected to be an integral type).
 */
template <typename Element, typename... Args>
class Subsets<UintSet, Element, Args...> {
public:
  class Iterator : public std::iterator<std::forward_iterator_tag, UintSet<Element, Args...>> {
  public:
    Iterator(
      const UintSet<Element, Args...>& given,
      const std::vector<bool>& candidates,
      const bool valid = true
    ) : m_given(given),
        m_candidates(candidates),
        m_curr(0),
        m_valid(valid)
    {
    }

    Iterator&
    operator++()
    {
      if (m_valid) {
        m_valid = std::prev_permutation(m_candidates.begin(), m_candidates.end());
        ++m_curr;
      }
      return *this;
    }

    bool
    operator==(const Iterator& other) const
    {
      if (m_valid && other.m_valid) {
        // If both the iterators are valid then
        // check the permutation count
        // XXX: This is done in order to avoid relatively
        // expensive vector comparisons
        return m_curr == other.m_curr;
      }
      else {
        // Otherwise, the iterators are equal if
        // both of the are invalidated
        return !(m_valid || other.m_valid);
      }
    }

    bool
    operator!=(const Iterator& other) const
    {
      return !(*this == other);
    }

    UintSet<Element, Args...>
    operator*() const
    {
      return m_given.subset(m_candidates);
    }

  private:
    const UintSet<Element, Args...>& m_given;
    std::vector<bool> m_candidates;
    size_t m_curr;
    bool m_valid;
  };

public:
  // Required typedefs
  using value_type = Element;
  using iterator = Iterator;

public:
  Subsets(const UintSet<Element, Args...>& given, const uint32_t k)
    : m_given(given),
      m_candidates(given.size(), false)
  {
    auto it = m_candidates.begin();
    for (auto i = 0u; i < k; ++i, ++it) {
      *it = true;
    }
  }

  Iterator
  begin() const
  {
    return Iterator(m_given, m_candidates);
  }

  Iterator
  end() const
  {
    return Iterator(m_given, m_candidates, false);
  }

private:
  const UintSet<Element, Args...>& m_given;
  std::vector<bool> m_candidates;
}; // class Subsets

template <typename Element, typename USet, template <typename> class Functor, typename ReduceType>
void
uintset_allreduce(
  USet& set,
  Functor<ReduceType> func,
  const mxx::comm& comm
)
{
  auto size = (set.max() + 63) / 64;
  auto b = new ReduceType[size];
  mxx::allreduce((*set).b, size, b, func, comm);
  memcpy((*set).b, b, size * sizeof(ReduceType));
  set.m_size = static_cast<Element>(set_size(*set));
  delete[] b;
}

template <typename Element, typename USet, template <typename> class Functor, typename ReduceType>
void
uintset_allreduce_indexed(
  std::unordered_map<Element, USet>& indexedSets,
  const USet& allIndices,
  const Element max,
  Functor<ReduceType> func,
  const mxx::comm& comm
)
{
  auto blockSize = static_cast<uint32_t>(max + 63) / 64;
  auto totalSize = blockSize * allIndices.max();
  auto bitSets = new ReduceType[totalSize]();
  auto b = bitSets;
  for (const auto x : allIndices) {
    auto it = indexedSets.find(x);
    if (it != indexedSets.end()) {
      memcpy(b, (*(it->second)).b, blockSize * sizeof(ReduceType));
    }
    b += blockSize;
  }
  b = bitSets;
  while (totalSize > 0) {
    auto reduceSize = std::min(totalSize, ALLREDUCE_MAXSIZE);
    mxx::allreduce(static_cast<ReduceType*>(MPI_IN_PLACE), reduceSize, b, func, comm);
    b += reduceSize;
    totalSize -= reduceSize;
  }
  b = bitSets;
  for (const auto x : allIndices) {
    auto it = indexedSets.find(x);
    if (it != indexedSets.end()) {
      memcpy((*(it->second)).b, b, blockSize * sizeof(ReduceType));
      (it->second).m_size = static_cast<Element>(set_size(*(it->second)));
    }
    b += blockSize;
  }
  delete[] bitSets;
}

// Definition of all the operations on UintSet
#define DEFINE_UINT_SET_OPERATIONS(Element, N) \
template <> \
UintSet<Element, std::integral_constant<int, N>> \
set_init<UintSet<Element, std::integral_constant<int, N>>, Element>( \
  UintSet<Element, std::integral_constant<int, N>>&& set, \
  const Element max \
) \
{ \
  set.m_max = max; \
  return set; \
} \
 \
template <> \
bool \
set_contains<UintSet<Element, std::integral_constant<int, N>>, Element>( \
  const UintSet<Element, std::integral_constant<int, N>>& set, \
  const Element value \
) \
{ \
  return set.contains(value); \
} \
 \
template <> \
UintSet<Element, std::integral_constant<int, N>> \
set_union<UintSet<Element, std::integral_constant<int, N>>>( \
  const UintSet<Element, std::integral_constant<int, N>>& first, \
  const UintSet<Element, std::integral_constant<int, N>>& second \
) \
{ \
  auto result = *first | *second; \
  auto max = std::max(first.max(), second.max()); \
  return UintSet<Element, std::integral_constant<int, N>>(result, max); \
} \
 \
template <> \
UintSet<Element, std::integral_constant<int, N>> \
set_intersection<UintSet<Element, std::integral_constant<int, N>>>( \
  const UintSet<Element, std::integral_constant<int, N>>& first, \
  const UintSet<Element, std::integral_constant<int, N>>& second \
) \
{ \
  auto result = *first & *second; \
  auto max = std::max(first.max(), second.max()); \
  return UintSet<Element, std::integral_constant<int, N>>(result, max); \
} \
 \
template <> \
UintSet<Element, std::integral_constant<int, N>> \
set_difference<UintSet<Element, std::integral_constant<int, N>>>( \
  const UintSet<Element, std::integral_constant<int, N>>& first, \
  const UintSet<Element, std::integral_constant<int, N>>& second \
) \
{ \
  auto result = set_diff(*first, *second); \
  auto max = std::max(first.max(), second.max()); \
  return UintSet<Element, std::integral_constant<int, N>>(result, max); \
} \
 \
template <> \
void \
set_bcast<UintSet<Element, std::integral_constant<int, N>>>( \
  UintSet<Element, std::integral_constant<int, N>>& set, \
  const int root, \
  const mxx::comm& comm \
) \
{ \
  auto size = (set.max() + 63) / 64; \
  mxx::bcast((*set).b, size, root, comm); \
  if (comm.rank() != root) { \
    set.m_size = static_cast<Element>(set_size(*set)); \
  } \
} \
 \
template <> \
void \
set_allunion<UintSet<Element, std::integral_constant<int, N>>>( \
  UintSet<Element, std::integral_constant<int, N>>& set, \
  const mxx::comm& comm \
) \
{ \
  using ReduceType = typename UintSet<Element, std::integral_constant<int, N>>::ReduceType; \
  uintset_allreduce<Element>(set, std::bit_or<ReduceType>(), comm); \
} \
 \
template <> \
void \
set_allunion_indexed<UintSet<Element, std::integral_constant<int, N>>>( \
  std::unordered_map<Element, UintSet<Element, std::integral_constant<int, N>>>& indexedSets, \
  const UintSet<Element, std::integral_constant<int, N>>& allIndices, \
  const Element max, \
  const mxx::comm& comm \
) \
{ \
  using ReduceType = typename UintSet<Element, std::integral_constant<int, N>>::ReduceType; \
  uintset_allreduce_indexed(indexedSets, allIndices, max, std::bit_or<ReduceType>(), comm); \
} \
 \
template <> \
void \
set_allintersect<UintSet<Element, std::integral_constant<int, N>>>( \
  UintSet<Element, std::integral_constant<int, N>>& set, \
  const mxx::comm& comm \
) \
{ \
  using ReduceType = typename UintSet<Element, std::integral_constant<int, N>>::ReduceType; \
  uintset_allreduce<Element>(set, std::bit_and<ReduceType>(), comm); \
} \
 \
template <> \
void \
set_allintersect_indexed<UintSet<Element, std::integral_constant<int, N>>>( \
  std::unordered_map<Element, UintSet<Element, std::integral_constant<int, N>>>& indexedSets, \
  const UintSet<Element, std::integral_constant<int, N>>& allIndices, \
  const Element max, \
  const mxx::comm& comm \
) \
{ \
  using ReduceType = typename UintSet<Element, std::integral_constant<int, N>>::ReduceType; \
  uintset_allreduce_indexed(indexedSets, allIndices, max, std::bit_and<ReduceType>(), comm); \
} \
 \
template <> \
std::vector<UintSet<Element, std::integral_constant<int, N>>> \
set_allgatherv<UintSet<Element, std::integral_constant<int, N>>>( \
  const std::vector<UintSet<Element, std::integral_constant<int, N>>>& mySets, \
  const std::vector<size_t>& allSizes, \
  const Element max, \
  const mxx::comm& comm \
) \
{ \
  using SetType = UintSet<Element, std::integral_constant<int, N>>; \
  using ReduceType = typename SetType::ReduceType; \
  auto blockSize = static_cast<uint32_t>(max + 63) / 64; \
  std::vector<size_t> allOffsets(allSizes.size()); \
  std::partial_sum(allSizes.begin(), allSizes.end(), allOffsets.begin()); \
  auto totalSize = static_cast<uint32_t>(blockSize * allOffsets.back()); \
  auto bitSets = new ReduceType[totalSize](); \
  auto myOffset = allOffsets[comm.rank()] - allSizes[comm.rank()]; \
  auto b = bitSets + (myOffset * blockSize); \
  for (const auto& set : mySets) { \
    memcpy(b, (*set).b, blockSize * sizeof(ReduceType)); \
    b += blockSize; \
  } \
  b = bitSets; \
  while (totalSize > 0) { \
    auto reduceSize = std::min(totalSize, ALLREDUCE_MAXSIZE); \
    mxx::allreduce(static_cast<ReduceType*>(MPI_IN_PLACE), reduceSize, b, std::bit_or<ReduceType>(), comm); \
    b += reduceSize; \
    totalSize -= reduceSize; \
  } \
  std::vector<SetType> allSets(allOffsets.back(), set_init(SetType(), max)); \
  b = bitSets; \
  for (auto& set : allSets) { \
    memcpy((*set).b, b, blockSize * sizeof(ReduceType)); \
    set.m_size = static_cast<Element>(set_size(*set)); \
    b += blockSize; \
  } \
  delete[] bitSets; \
  return allSets; \
}

DEFINE_UINT_SET_OPERATIONS(uint8_t, (maxSize<uint8_t>() >> 2))
DEFINE_UINT_SET_OPERATIONS(uint8_t, (maxSize<uint8_t>() >> 1))
DEFINE_UINT_SET_OPERATIONS(uint8_t, maxSize<uint8_t>())
DEFINE_UINT_SET_OPERATIONS(uint16_t, (maxSize<uint16_t>() >> 7))
DEFINE_UINT_SET_OPERATIONS(uint16_t, (maxSize<uint16_t>() >> 6))
DEFINE_UINT_SET_OPERATIONS(uint16_t, (maxSize<uint16_t>() >> 5))
DEFINE_UINT_SET_OPERATIONS(uint16_t, (maxSize<uint16_t>() >> 4))
DEFINE_UINT_SET_OPERATIONS(uint16_t, (maxSize<uint16_t>() >> 3))
DEFINE_UINT_SET_OPERATIONS(uint16_t, (maxSize<uint16_t>() >> 2))
DEFINE_UINT_SET_OPERATIONS(uint16_t, (maxSize<uint16_t>() >> 1))
DEFINE_UINT_SET_OPERATIONS(uint16_t, maxSize<uint16_t>())

/**
 * @brief Function for getting the output represention of a set.
 */
template <typename Element, typename Size>
std::ostream&
operator<<(
  std::ostream& stream,
  const UintSet<Element, Size>& set
)
{
  auto first = true;
  for (const auto elem : set) {
    stream << (first ? "" : ";") << static_cast<uint32_t>(elem);
    first = false;
  }
  return stream;
}

#endif // DETAIL_UINTSETUTILS_HPP_
