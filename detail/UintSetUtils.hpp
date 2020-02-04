/**
 * @file UintSetUtils.hpp
 * @brief Implementation of functions for operations on UintSet.
 */
#ifndef DETAIL_UINTSETUTILS_HPP_
#define DETAIL_UINTSETUTILS_HPP_

#include "../UintSet.hpp"

#include "mxx/reduction.hpp"

#include <algorithm>
#include <type_traits>


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
    auto i = 0u;
    for (auto it = m_candidates.begin(); i < k; ++it, ++i) {
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
  UintSet<Element, std::integral_constant<int, N>> result; \
  *result = *first | *second; \
  return result; \
} \
 \
template <> \
UintSet<Element, std::integral_constant<int, N>> \
set_difference<UintSet<Element, std::integral_constant<int, N>>>( \
  const UintSet<Element, std::integral_constant<int, N>>& first, \
  const UintSet<Element, std::integral_constant<int, N>>& second \
) \
{ \
  UintSet<Element, std::integral_constant<int, N>> result; \
  *result = set_diff(*first, *second); \
  return result; \
} \
 \
template <> \
void \
set_bcast<UintSet<Element, std::integral_constant<int, N>>>( \
  UintSet<Element, std::integral_constant<int, N>>& set, \
  const int source, \
  const mxx::comm& comm \
) \
{ \
  auto size = (set.max() + 63) / 64; \
  mxx::bcast(&set.m_set.b[0], size, source, comm); \
} \
 \
template <> \
void \
set_allunion<UintSet<Element, std::integral_constant<int, N>>>( \
  UintSet<Element, std::integral_constant<int, N>>& set, \
  const mxx::comm& comm \
) \
{ \
  auto size = (set.max() + 63) / 64; \
  using ReduceType = std::remove_reference<decltype(set.m_set.b[0])>::type; \
  auto b = new ReduceType[size]; \
  mxx::allreduce(set.m_set.b, size, b, std::bit_or<ReduceType>(), comm); \
  memcpy(set.m_set.b, b, size * sizeof(ReduceType)); \
  delete[] b; \
} \
 \
template <> \
void \
set_allintersect<UintSet<Element, std::integral_constant<int, N>>>( \
  UintSet<Element, std::integral_constant<int, N>>& set, \
  const mxx::comm& comm \
) \
{ \
  auto size = (set.max() + 63) / 64; \
  using ReduceType = std::remove_reference<decltype(set.m_set.b[0])>::type; \
  auto b = new ReduceType[size]; \
  mxx::allreduce(set.m_set.b, size, b, std::bit_and<ReduceType>(), comm); \
  memcpy(set.m_set.b, b, size * sizeof(ReduceType)); \
  delete[] b; \
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
  stream << "{";
  for (auto elem: set) {
    stream << static_cast<uint32_t>(elem) << ",";
  }
  if (set.size() > 0) {
    stream << '\b';
  }
  stream << "}";
  return stream;
}

#endif // DETAIL_UINTSETUTILS_HPP_
