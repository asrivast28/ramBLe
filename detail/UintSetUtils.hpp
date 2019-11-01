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
 * @tparam Element Type of the variable (expected to be an integral type).
 */
template <typename Element>
class SubsetIterator<UintSet<Element>, Element> {
public:
  SubsetIterator(const UintSet<Element>& given, const uint32_t k)
    : m_given(given),
      m_subset(),
      m_candidates(given.max(), false),
      m_valid(!((given.size() == 0) || (given.size() == k)))
  {
    auto i = 0u;
    for (auto it = m_candidates.begin(); i < k; ++it, ++i) {
      *it = true;
    }
    m_subset = UintSet<Element>(m_candidates, m_subset.max());
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
      m_subset = UintSet<Element>(m_candidates, m_subset.max());
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

  const UintSet<Element>&
  get() const
  {
    return (!m_valid) ? m_given: m_subset;
  }

private:
  const UintSet<Element>& m_given;
  UintSet<Element> m_subset;
  std::vector<bool> m_candidates;
  bool m_valid;
}; // class SubsetIterator


// Definition of all the operations on UintSet
#define DEFINE_UINT_SET_OPERATIONS(type) \
template <> \
UintSet<type> \
set_init<UintSet<type>, type>( \
  UintSet<type>&& set, \
  const type max \
) \
{ \
  set.m_max = max; \
  return set; \
} \
 \
template <> \
bool \
set_contains<UintSet<type>, type>( \
  const UintSet<type>& set, \
  const type value \
) \
{ \
  return set.contains(value); \
} \
 \
template <> \
UintSet<type> \
set_union<UintSet<type>>( \
  const UintSet<type>& first, \
  const UintSet<type>& second \
) \
{ \
  UintSet<type> result; \
  *result = *first | *second; \
  return result; \
} \
 \
template <> \
UintSet<type> \
set_difference<UintSet<type>>( \
  const UintSet<type>& first, \
  const UintSet<type>& second \
) \
{ \
  UintSet<type> result; \
  *result = set_diff(*first, *second); \
  return result; \
} \
 \
template <> \
bool \
is_subset<UintSet<type>>( \
  const UintSet<type>& first, \
  const UintSet<type>& second \
) \
{ \
  return is_superset(*second, *first); \
}

DEFINE_UINT_SET_OPERATIONS(uint8_t)
DEFINE_UINT_SET_OPERATIONS(uint16_t)

/**
 * @brief Function for getting the output represention of a set.
 */
template <typename Element>
std::ostream&
operator<<(
  std::ostream& stream,
  const UintSet<Element>& set
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
