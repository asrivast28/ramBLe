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
class Subsets<UintSet<Element>, Element> {
public:
  class Iterator : public std::iterator<std::forward_iterator_tag, UintSet<Element>> {
  public:
    Iterator(
      const UintSet<Element>& given,
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

    UintSet<Element>
    operator*() const
    {
      return m_given.subset(m_candidates);
    }

  private:
    const UintSet<Element>& m_given;
    std::vector<bool> m_candidates;
    size_t m_curr;
    bool m_valid;
  };

public:
  // Required typedefs
  using value_type = Element;
  using iterator = Iterator;

public:
  Subsets(const UintSet<Element>& given, const uint32_t k)
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
  const UintSet<Element>& m_given;
  std::vector<bool> m_candidates;
}; // class Subsets


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
