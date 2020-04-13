/**
 * @file UintSet.hpp
 * @brief Implementation of the UintSet and related classes.
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
#ifndef DETAIL_UINTSET_HPP_
#define DETAIL_UINTSET_HPP_

#include "bit_util.hpp"

#include "utils/Logging.hpp"

#include <iterator>


template <typename Size>
class UintTypeTrait {
public:
  using Set = uint_type<Size::value>;

public:
  static
  constexpr
  uint32_t
  max()
  {
    // The maximum variable index that this set can hold
    return static_cast<uint32_t>(set_max_size<Set>() - 1);
  }
};

template <typename Element>
constexpr
int
maxSize() {
  return static_cast<int>(std::numeric_limits<Element>::max() / 64) + 1;
}

template <typename Element, typename Size>
/**
 * @brief Enumerates the elements of the UintSet.
 */
class UintSet<Element, Size>::Enumerator : public std::iterator<std::forward_iterator_tag, Element> {
public:
  using Impl = typename UintTypeTrait<Size>::Set;

public:
  Enumerator(const Impl& set, const Element max)
    : m_state(set),
      m_maxN(max >> 6),
      m_currN(0),
      m_max(max),
      m_curr(0)
  {
    next();
  }

  Enumerator(const Impl& set, const Element max, const Element curr)
    : m_state((curr < max) ? set : set_empty<Impl>()),
      m_maxN(max >> 6),
      m_currN(curr >> 6),
      m_max(max),
      m_curr((curr < max) ? 0 : max)
  {
    if (curr < max) {
      // Empty out all the bitsets smaller than the index of
      // the current element
      for (auto i = 0; i < m_currN; ++i) {
        m_state.b[i] = 0;
      }
      // Now, unset the bits smaller than or equal to the one
      // corresponding to the current element
      m_curr = (m_currN << 6);
      while (m_curr <= curr) {
        auto x = lsb(m_state.b[m_currN]);
        m_curr = static_cast<Element>((m_currN << 6) + x);
        m_state.b[m_currN] = set_remove(m_state.b[m_currN], x);
      }
    }
  }

  Enumerator&
  operator++()
  {
    next();
    return *this;
  }

  bool
  operator==(const Enumerator& other) const
  {
    return m_curr == *other;
  }

  bool
  operator!=(const Enumerator& other) const
  {
    return m_curr != *other;
  }

  Element
  operator*() const
  {
    return m_curr;
  }

private:
  void
  next()
  {
    if (m_curr < m_max) {
      // Find the first index with a non-empty bitset
      while ((m_currN <= m_maxN) && is_emptyset(m_state.b[m_currN])) {
        ++m_currN;
      }
      if (m_currN <= m_maxN) {
        // There must be at least one set bit in the current index
        // Get the index of the lowest such bit
        auto x = lsb(m_state.b[m_currN]);
        // Multiply the current bitset index by 64 and then add the
        // index of the set bit to get the current element
        m_curr = static_cast<Element>((m_currN << 6) + x);
        // Set the bit to zero so that the next bit can be found
        m_state.b[m_currN] = set_remove(m_state.b[m_currN], x);
      }
      else {
        // No new elements were found
        m_curr = m_max;
      }
    }
  }

private:
  Impl m_state;
  int m_maxN;
  int m_currN;
  Element m_max;
  Element m_curr;
};


template <typename Element, typename Size>
/**
 * @brief The capacity of a set of type UintSet<Element, Size>.
 */
constexpr
Element
UintSet<Element, Size>::capacity(
)
{
  return static_cast<Element>(UintTypeTrait<Size>::max());
}

template <typename Element, typename Size>
UintSet<Element, Size>::UintSet(
  const Element max
) : m_set(set_empty<Impl>()),
    m_max(max),
    m_size()
{
}

template <typename Element, typename Size>
UintSet<Element, Size>::UintSet(
  const std::initializer_list<Element>& s,
  const Element max
) : m_set(as_set<Impl>(s.begin(), s.end())),
    m_max(max),
    m_size()
{
  m_size = static_cast<Element>(set_size(m_set));
}

template <typename Element, typename Size>
UintSet<Element, Size>::UintSet(
  const typename UintSet<Element, Size>::iterator& first,
  const typename UintSet<Element, Size>::iterator& last,
  const Element max
) : m_set(as_set<Impl>(first, last)),
    m_max(max),
    m_size()
{
  m_size = static_cast<Element>(set_size(m_set));
}

template <typename Element, typename Size>
UintSet<Element, Size>::UintSet(
  const typename std::vector<Element>::iterator& first,
  const typename std::vector<Element>::iterator& last,
  const Element max
) : m_set(as_set<Impl>(first, last)),
    m_max(max),
    m_size()
{
  m_size = static_cast<Element>(set_size(m_set));
}

template <typename Element, typename Size>
const typename UintSet<Element, Size>::Impl&
UintSet<Element, Size>::operator*(
) const
{
  return m_set;
}

template <typename Element, typename Size>
typename UintSet<Element, Size>::Impl&
UintSet<Element, Size>::operator*(
)
{
  return m_set;
}

template <typename Element, typename Size>
bool
UintSet<Element, Size>::operator==(
  const UintSet<Element, Size>& other
) const
{
  return (m_set == other.m_set);
}

template <typename Element, typename Size>
bool
UintSet<Element, Size>::operator!=(
  const UintSet<Element, Size>& other
) const
{
  return (m_set != other.m_set);
}

template <typename Element, typename Size>
typename UintSet<Element, Size>::iterator
UintSet<Element, Size>::insert(
  const Element x
)
{
  LOG_MESSAGE_IF(x > m_max, error, "Inserting a value (%d) which is greater than the max (%d)", static_cast<int>(x), static_cast<int>(m_max));
  if (!contains(x)) {
    ++m_size;
  }
  m_set = set_add(std::move(m_set), static_cast<int>(x));
  return end();
}

template <typename Element, typename Size>
typename UintSet<Element, Size>::iterator
UintSet<Element, Size>::insert(
  const iterator&,
  const Element x
)
{
  LOG_MESSAGE_IF(x > m_max, error, "Inserting a value (%d) which is greater than the max (%d)", static_cast<int>(x), static_cast<int>(m_max));
  if (!contains(x)) {
    ++m_size;
  }
  m_set = set_add(std::move(m_set), static_cast<int>(x));
  return end();
}

template <typename Element, typename Size>
void
UintSet<Element, Size>::erase(
  const Element x
)
{
  LOG_MESSAGE_IF(!contains(x), warning, "Removing a value (%d) which does not exist in the set", static_cast<int>(x));
  if ((m_size != 0) && contains(x)) {
    --m_size;
  }
  m_set = set_remove(std::move(m_set), static_cast<int>(x));
}

template <typename Element, typename Size>
Element
UintSet<Element, Size>::max(
) const
{
  return m_max;
}

template <typename Element, typename Size>
Element
UintSet<Element, Size>::size(
) const
{
  if (m_size == 0) {
    m_size = static_cast<Element>(set_size(m_set));
  }
  return m_size;
}

template <typename Element, typename Size>
bool
UintSet<Element, Size>::empty(
) const
{
  return is_emptyset(m_set);
}

template <typename Element, typename Size>
bool
UintSet<Element, Size>::contains(
  const Element x
) const
{
  return in_set(m_set, static_cast<int>(x));
}

template <typename Element, typename Size>
typename UintSet<Element, Size>::iterator
UintSet<Element, Size>::begin(
) const
{
  return Enumerator(m_set, m_max);
}

template <typename Element, typename Size>
typename UintSet<Element, Size>::iterator
UintSet<Element, Size>::end(
) const
{
  return Enumerator(m_set, m_max, m_max);
}

template <typename Element, typename Size>
typename UintSet<Element, Size>::iterator
UintSet<Element, Size>::find(
  const Element x
) const
{
  return contains(x) ? Enumerator(m_set, m_max, x): end();
}

template <typename Element, typename Size>
UintSet<Element, Size>
UintSet<Element, Size>::subset(
  const std::vector<bool>& bitset
) const
{
  UintSet subset(m_max);
  auto it = this->begin();
  for (auto i = 0u; i < bitset.size(); ++i, ++it) {
    if (bitset[i]) {
      subset.insert(*it);
    }
  }
  return subset;
}

#endif // DETAIL_UINTSET_HPP_
