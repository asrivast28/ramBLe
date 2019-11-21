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
 * @brief Iterator for UintSet.
 */
class UintSet<Element, N>::Iterator : public std::iterator<std::forward_iterator_tag, Element> {
public:
  using Set = typename UintTypeTrait<N>::Set;

public:
  Iterator(const Set* set, const Element max, const Element curr = UintSet<Element, N>::capacity())
    : m_set(set),
      m_max(max),
      m_curr((curr > max)? max: curr)
  {
    next();
  }

  Iterator&
  operator++()
  {
    ++m_curr;
    next();
    return *this;
  }

  bool
  operator==(const Iterator& other) const
  {
    return m_curr == *other;
  }

  bool
  operator!=(const Iterator& other) const
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
    while ((m_curr < m_max) && !in_set(*m_set, static_cast<int>(m_curr))) {
      ++m_curr;
    }
  }

private:
  const Set* m_set;
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
  const typename UintSet<Element, N>::Iterator& first,
  const typename UintSet<Element, N>::Iterator& last,
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
typename UintSet<Element, N>::Iterator
UintSet<Element, N>::insert(
  const Element x
)
{
  LOG_MESSAGE_IF(x > m_max, error, "Inserting a value (%d) which is greater than the max (%d)", static_cast<int>(x), static_cast<int>(m_max));
  m_set = set_add(std::move(m_set), static_cast<int>(x));
  return end();
}

template <typename Element, int N>
typename UintSet<Element, N>::Iterator
UintSet<Element, N>::insert(
  const Iterator&,
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
typename UintSet<Element, N>::Iterator
UintSet<Element, N>::begin(
) const
{
  return Iterator(&m_set, m_max, 0);
}

template <typename Element, int N>
typename UintSet<Element, N>::Iterator
UintSet<Element, N>::end(
) const
{
  return Iterator(&m_set, m_max);
}

template <typename Element, int N>
typename UintSet<Element, N>::Iterator
UintSet<Element, N>::find(
  const Element x
) const
{
  return contains(x) ? Iterator(&m_set, m_max, x): end();
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
