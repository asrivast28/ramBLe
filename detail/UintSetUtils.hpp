/**
 * @file UintSetUtils.hpp
 * @brief Implementation of functions for operations on UintSet.
 */
#ifndef DETAIL_UINTSETUTILS_HPP_
#define DETAIL_UINTSETUTILS_HPP_

#include "../UintSet.hpp"

#include <algorithm>


/**
 * @brief Class that iterates over all the subsets of the given size of a given UintSet.
 *
 * @tparam SetType Type of the set container.
 * @tparam ValueType Type of the variable (expected to be an integral type).
 */
template <typename ValueType>
class SubsetIterator<UintSet<ValueType>, ValueType> {
public:
  static constexpr int N = UintTypeTrait<ValueType>::N;

public:
  SubsetIterator(const UintSet<ValueType>& given, const uint32_t k)
    : m_given(given),
      m_subset(),
      m_candidates(given.max(), false),
      m_valid(!((given.size() == 0) || (given.size() == k)))
  {
    auto i = 0u;
    for (auto it = m_candidates.begin(); i < k; ++it, ++i) {
      *it = true;
    }
    m_subset = UintSet<ValueType>(m_candidates, m_subset.max());
    if (!is_subset(m_subset, m_given)) {
      next();
    }
  }

  void
  next()
  {
    if (!m_valid) {
      return;
    }
    m_valid = false;
    while (std::prev_permutation(m_candidates.begin(), m_candidates.end())) {
      m_subset = UintSet<ValueType>(m_candidates, m_subset.max());
      if (is_subset(m_subset, m_given)) {
        m_valid = true;
        break;
      }
    }
  }

  bool
  valid()
  {
    return m_valid;
  }

  const UintSet<ValueType>&
  get() const
  {
    return (!m_valid) ? m_given: m_subset;
  }

private:
  const UintSet<ValueType>& m_given;
  UintSet<ValueType> m_subset;
  std::vector<bool> m_candidates;
  bool m_valid;
}; // class SubsetIterator

template <>
UintSet<uint8_t>
set_init<UintSet<uint8_t>, uint8_t>(
  UintSet<uint8_t>&& set,
  const uint8_t max
)
{
  set.m_max = max;
  return set;
}

template <>
UintSet<uint16_t>
set_init<UintSet<uint16_t>, uint16_t>(
  UintSet<uint16_t>&& set,
  const uint16_t max
)
{
  set.m_max = max;
  return set;
}

template <>
bool
set_contains<UintSet<uint8_t>, uint8_t>(
  const UintSet<uint8_t>& set,
  const uint8_t value
)
{
  return set.contains(value);
}

template <>
bool
set_contains<UintSet<uint16_t>, uint16_t>(
  const UintSet<uint16_t>& set,
  const uint16_t value
)
{
  return set.contains(value);
}

template <>
UintSet<uint8_t>
set_union<UintSet<uint8_t>>(
  const UintSet<uint8_t>& first,
  const UintSet<uint8_t>& second
)
{
  UintSet<uint8_t> result;
  *result = *first | *second;
  return result;
}

template <>
UintSet<uint16_t>
set_union<UintSet<uint16_t>>(
  const UintSet<uint16_t>& first,
  const UintSet<uint16_t>& second
)
{
  UintSet<uint16_t> result;
  *result = *first | *second;
  return result;
}

template <>
UintSet<uint8_t>
set_difference<UintSet<uint8_t>>(
  const UintSet<uint8_t>& first,
  const UintSet<uint8_t>& second
)
{
  UintSet<uint8_t> result;
  *result = set_diff(*first, *second);
  return result;
}

template <>
UintSet<uint16_t>
set_difference<UintSet<uint16_t>>(
  const UintSet<uint16_t>& first,
  const UintSet<uint16_t>& second
)
{
  UintSet<uint16_t> result;
  *result = set_diff(*first, *second);
  return result;
}

template <>
bool
is_subset<UintSet<uint8_t>>(
  const UintSet<uint8_t>& first,
  const UintSet<uint8_t>& second
)
{
  return is_superset(*second, *first);
}

template <>
bool
is_subset<UintSet<uint16_t>>(
  const UintSet<uint16_t>& first,
  const UintSet<uint16_t>& second
)
{
  return is_superset(*second, *first);
}

/**
 * @brief Function for getting the output represention of a set.
 */
template <typename ValueType>
std::ostream&
operator<<(
  std::ostream& stream,
  const UintSet<ValueType>& set
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
