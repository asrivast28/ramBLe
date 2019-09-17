/**
 * @file SetUtils.hpp
 * @brief Declaration of functions for common set operations.
 */
#ifndef SETUTILS_HPP_
#define SETUTILS_HPP_

#include <algorithm>
#include <iterator>
#include <vector>

/**
 * @brief Class that provides a lightweight subset wrapper over an STL container.
 *
 * @tparam IteratorType The type of the iterator.
 *
 * This class allows iterating over a contiguous subset of a container as defined
 * by the first iterator, the last iterator, and the size of the slice.
 */
template <typename IteratorType>
class SubsetWrapper {
public:
  SubsetWrapper(const IteratorType first, const IteratorType last, const uint32_t size)
    : m_begin(first),
      m_end(last),
      m_size(size)
  {
  }

  const IteratorType&
  begin() const
  {
    return m_begin;
  }

  const IteratorType&
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
  const IteratorType m_begin;
  const IteratorType m_end;
  const uint32_t m_size;
}; // class SubsetWrapper

/**
 * @brief Class that iterates over all the subsets of the given size of a given set.
 *
 * @tparam SetType Type of the set container.
 * @tparam ValueType Type of the variable (expected to be an integral type).
 */
template <typename SetType, typename ValueType>
class SubsetIterator {
public:
  SubsetIterator(const SetType& given, const uint32_t k)
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

  SubsetWrapper<typename std::vector<ValueType>::iterator>
  subset() const
  {
    return SubsetWrapper<typename std::vector<ValueType>::iterator>(m_first, m_current, m_k);
  }

private:
  std::vector<ValueType> m_given;
  const typename std::vector<ValueType>::iterator m_first;
  const typename std::vector<ValueType>::iterator m_current;
  const typename std::vector<ValueType>::iterator m_last;
  const uint32_t m_k;
  bool m_valid;
}; // class SubsetIterator


/**
 * @brief Function for initializing a given set.
 */
template <typename SetType, typename ValueType>
SetType
set_init(SetType&&, const ValueType);

/**
 * @brief Function for checking if a given set contains a value.
 */
template <typename SetType, typename ValueType>
bool
set_contains(const SetType&, const ValueType);

/**
 * @brief Function for getting the union of two given sets.
 */
template <typename SetType>
SetType
set_union(const SetType&, const SetType&);

#include "detail/StdSetUtils.hpp"
#include "detail/UintSetUtils.hpp"

#endif // SETUTILS_HPP_
