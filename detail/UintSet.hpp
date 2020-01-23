/**
 * @file UintSet.hpp
 * @brief Implementation of the UintSet and related classes.
 */
#ifndef DETAIL_UINTSET_HPP_
#define DETAIL_UINTSET_HPP_

#include "bit_util.hpp"

#include "utils/Logging.hpp"

#include <iterator>


template <int N>
class UintTypeTrait {
public:
  using Set = uint_type<N>;

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
maxN() {
  return static_cast<int>(std::numeric_limits<Element>::max() / 64) + 1;
}

template <typename Element, int N>
/**
 * @brief Enumerates the elements of the UintSet.
 */
class UintSet<Element, N>::Enumerator : public std::iterator<std::forward_iterator_tag, Element> {
public:
  using Set = typename UintTypeTrait<N>::Set;

public:
  Enumerator(const Set& set, const Element max)
    : m_state(set),
      m_maxN(max >> 6),
      m_currN(0),
      m_max(max),
      m_curr(0)
  {
    next();
  }

  Enumerator(const Set& set, const Element max, const Element curr)
    : m_state((curr < max) ? set : set_empty<Set>()),
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
  Set m_state;
  int m_maxN;
  int m_currN;
  Element m_max;
  Element m_curr;
};


template <typename Element, int N>
/**
 * @brief The capacity of a set of type UintSet<Element, N>.
 */
constexpr
Element
UintSet<Element, N>::capacity(
)
{
  return static_cast<Element>(UintTypeTrait<N>::max());
}

template <typename Element, int N>
UintSet<Element, N>::UintSet(
  const Element max
) : m_set(set_empty<Set>()),
    m_max(max)
{
}

template <typename Element, int N>
UintSet<Element, N>::UintSet(
  const std::initializer_list<Element>& s,
  const Element max
) : m_set(as_set<Set>(s.begin(), s.end())),
    m_max(max)
{
}

template <typename Element, int N>
UintSet<Element, N>::UintSet(
  const typename UintSet<Element, N>::iterator& first,
  const typename UintSet<Element, N>::iterator& last,
  const Element max
) : m_set(as_set<Set>(first, last)),
    m_max(max)
{
}

template <typename Element, int N>
UintSet<Element, N>::UintSet(
  const typename std::vector<Element>::iterator& first,
  const typename std::vector<Element>::iterator& last,
  const Element max
) : m_set(as_set<Set>(first, last)),
    m_max(max)
{
}

template <typename Element, int N>
const typename UintSet<Element, N>::Set&
UintSet<Element, N>::operator*(
) const
{
  return m_set;
}

template <typename Element, int N>
typename UintSet<Element, N>::Set&
UintSet<Element, N>::operator*(
)
{
  return m_set;
}

template <typename Element, int N>
bool
UintSet<Element, N>::operator==(
  const UintSet<Element, N>& other
) const
{
  return (m_set == other.m_set);
}

template <typename Element, int N>
bool
UintSet<Element, N>::operator!=(
  const UintSet<Element, N>& other
) const
{
  return (m_set != other.m_set);
}

template <typename Element, int N>
typename UintSet<Element, N>::iterator
UintSet<Element, N>::insert(
  const Element x
)
{
  LOG_MESSAGE_IF(x > m_max, error, "Inserting a value (%d) which is greater than the max (%d)", static_cast<int>(x), static_cast<int>(m_max));
  m_set = set_add(std::move(m_set), static_cast<int>(x));
  return end();
}

template <typename Element, int N>
typename UintSet<Element, N>::iterator
UintSet<Element, N>::insert(
  const iterator&,
  const Element x
)
{
  LOG_MESSAGE_IF(x > m_max, error, "Inserting a value (%d) which is greater than the max (%d)", static_cast<int>(x), static_cast<int>(m_max));
  m_set = set_add(std::move(m_set), static_cast<int>(x));
  return end();
}

template <typename Element, int N>
void
UintSet<Element, N>::erase(
  const Element x
)
{
  LOG_MESSAGE_IF(!contains(x), warning, "Removing a value (%d) which does not exist in the set", static_cast<int>(x));
  m_set = set_remove(std::move(m_set), static_cast<int>(x));
}

template <typename Element, int N>
Element
UintSet<Element, N>::max(
) const
{
  return m_max;
}

template <typename Element, int N>
uint32_t
UintSet<Element, N>::size(
) const
{
  return static_cast<Element>(set_size(m_set));
}

template <typename Element, int N>
bool
UintSet<Element, N>::empty(
) const
{
  return is_emptyset(m_set);
}

template <typename Element, int N>
bool
UintSet<Element, N>::contains(
  const Element x
) const
{
  return in_set(m_set, static_cast<int>(x));
}

template <typename Element, int N>
typename UintSet<Element, N>::iterator
UintSet<Element, N>::begin(
) const
{
  return Enumerator(m_set, m_max);
}

template <typename Element, int N>
typename UintSet<Element, N>::iterator
UintSet<Element, N>::end(
) const
{
  return Enumerator(m_set, m_max, m_max);
}

template <typename Element, int N>
typename UintSet<Element, N>::iterator
UintSet<Element, N>::find(
  const Element x
) const
{
  return contains(x) ? Enumerator(m_set, m_max, x): end();
}

template <typename Element, int N>
UintSet<Element, N>
UintSet<Element, N>::subset(
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
