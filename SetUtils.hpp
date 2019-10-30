/**
 * @file SetUtils.hpp
 * @brief Declaration of functions for common set operations.
 */
#ifndef SETUTILS_HPP_
#define SETUTILS_HPP_

#include "mxx/comm.hpp"

#include <algorithm>
#include <iterator>
#include <vector>
#include <ostream>

/**
 * @brief Class that provides a lightweight subset wrapper over an STL container.
 *
 * @tparam Iterator The type of the iterator.
 *
 * This class allows iterating over a contiguous subset of a container as defined
 * by the first iterator, the last iterator, and the size of the slice.
 */
template <typename Iterator>
class SubsetWrapper {
public:
  SubsetWrapper(const Iterator first, const Iterator last, const uint32_t size)
    : m_begin(first),
      m_end(last),
      m_size(size)
  {
  }

  const Iterator&
  begin() const
  {
    return m_begin;
  }

  const Iterator&
  end() const
  {
    return m_end;
  }

  uint32_t
  size() const
  {
    return m_size;
  }

private:
  const Iterator m_begin;
  const Iterator m_end;
  const uint32_t m_size;
}; // class SubsetWrapper

/**
 * @brief Class that iterates over all the subsets of the given size of a given set.
 *
 * @tparam Set Type of the set container.
 * @tparam Element Type of the variable (expected to be an integral type).
 */
template <typename Set, typename Element>
class SubsetIterator {
public:
  SubsetIterator(const Set& given, const uint32_t k)
    : m_given(given.begin(), given.end()),
      m_first(m_given.begin()),
      m_current(m_given.begin()+k),
      m_last(m_given.end()),
      m_k(k),
      m_valid(!((given.size() < 2) || (k == 0) || (k == given.size())))
  {
  }

  void
  next()
  {
    if (!m_valid) {
      return;
    }
    m_valid = false;
    auto it2 = (m_last-1);
    for (auto it1 = (m_current-1); (it1+1) != m_first; --it1) {
      if (*it1 < *it2) {
        auto j = m_current;
        while (*j <= *it1) {
          ++j;
        }
        std::iter_swap(it1, j);
        ++it1;
        ++j;
        std::rotate(it1, j, m_last);
        it2 = m_current;
        while (j != m_last) {
          ++j;
          ++it2;
        }
        std::rotate(m_current, it2, m_last);
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

  SubsetWrapper<typename std::vector<Element>::iterator>
  get() const
  {
    return SubsetWrapper<typename std::vector<Element>::iterator>(m_first, m_current, m_k);
  }

private:
  std::vector<Element> m_given;
  const typename std::vector<Element>::iterator m_first;
  const typename std::vector<Element>::iterator m_current;
  const typename std::vector<Element>::iterator m_last;
  const uint32_t m_k;
  bool m_valid;
}; // class SubsetIterator


/**
 * @brief Function for initializing a given set.
 */
template <typename Set, typename Element>
Set
set_init(Set&&, const Element);

/**
 * @brief Function for checking if a given set contains a value.
 */
template <typename Set, typename Element>
bool
set_contains(const Set&, const Element);

/**
 * @brief Function for getting the union of two given sets.
 */
template <typename Set>
Set
set_union(const Set&, const Set&);

/**
 * @brief Function for getting the difference of the second set from the first.
 */
template <typename Set>
Set
set_difference(const Set&, const Set&);

/**
 * @brief Function for checking if the first set is a subset of the second.
 */
template <typename Set>
bool
is_subset(
  const Set& first,
  const Set& second
)
{
  auto diff = set_difference(first, second);
  return diff.empty();
}

#include "detail/StdSetUtils.hpp"
#include "detail/UintSetUtils.hpp"

#endif // SETUTILS_HPP_
