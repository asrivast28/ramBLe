/**
 * @file UintSet.hpp
 * @brief Implementation of the UintSet and related classes.
 */
#ifndef DETAIL_UINTSET_HPP_
#define DETAIL_UINTSET_HPP_

#include "bit_util.hpp"

#include "utils/Logging.hpp"

#include <iterator>


/**
 * @brief Specialization of the type trait class.
 */
template <>
class UintTypeTrait<uint8_t> {
public:
  static constexpr int N = 1;

public:
  using Set = uint_type<N>;

public:
  static
  constexpr
  uint8_t
  max()
  {
    return static_cast<uint32_t>(set_max_size<Set>() - 1);
  }
};

/**
 * @brief Specialization of the type trait class.
 */
template <>
class UintTypeTrait<uint16_t> {
public:
  static constexpr int N = 2;

public:
  using Set = uint_type<N>;

public:
  static
  constexpr
  uint16_t
  max()
  {
    return static_cast<uint16_t>(set_max_size<Set>() - 1);
  }
};

/**
 * @brief Specialization of the type trait class.
 */
template <>
class UintTypeTrait<uint32_t> {
public:
  static constexpr int N = 4;

public:
  using Set = uint_type<N>;

public:
  static
  constexpr
  uint32_t
  max()
  {
    return static_cast<uint32_t>(set_max_size<Set>());
  }
};


template <typename Element>
/**
 * @brief Iterator for UintSet.
 */
class UintSet<Element>::Iterator : public std::iterator<std::forward_iterator_tag, Element> {
public:
  using Set = typename UintTypeTrait<Element>::Set;

public:
  Iterator(const Set* set, const Element max, const Element curr = UintSet<Element>::capacity())

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


template <typename Element>
/**
 * @brief The capacity of a set of type UintSet<Element>.
 */
constexpr
Element
UintSet<Element>::capacity(
)
{
  return UintTypeTrait<Element>::max();
}

template <typename Element>
UintSet<Element>::UintSet(
  const Element max
) : m_set(set_empty<Set>()),
    m_max(max)
{
}

template <typename Element>
UintSet<Element>::UintSet(
  const std::initializer_list<Element>& s,
  const Element max
) : m_set(as_set<Set>(s.begin(), s.end())),
    m_max(max)
{
}

template <typename Element>
UintSet<Element>::UintSet(
  const typename UintSet<Element>::Iterator& first,
  const typename UintSet<Element>::Iterator& last,
  const Element max
) : m_set(as_set<Set>(first, last)),
    m_max(max)
{
}

template <typename Element>
UintSet<Element>::UintSet(
  const typename std::vector<Element>::iterator& first,
  const typename std::vector<Element>::iterator& last,
  const Element max
) : m_set(as_set<Set>(first, last)),
    m_max(max)
{
}

template <typename Element>
UintSet<Element>::UintSet(
  const std::vector<bool>& bitset,
  const Element max
) : m_set(set_empty<Set>()),
    m_max(max)
{
  auto i = 0u;
  for (auto elem: bitset) {
    if (elem) {
      insert(i);
    }
    ++i;
  }
}

template <typename Element>
const typename UintSet<Element>::Set&
UintSet<Element>::operator*(
) const
{
  return m_set;
}

template <typename Element>
typename UintSet<Element>::Set&
UintSet<Element>::operator*(
)
{
  return m_set;
}

template <typename Element>
bool
UintSet<Element>::operator==(
  const UintSet<Element>& other
) const
{
  return (m_set == other.m_set);
}

template <typename Element>
bool
UintSet<Element>::operator!=(
  const UintSet<Element>& other
) const
{
  return (m_set != other.m_set);
}

template <typename Element>
typename UintSet<Element>::Iterator
UintSet<Element>::insert(
  const Element x
)
{
  LOG_MESSAGE_IF(x > m_max, error, "Inserting a value (%d) which is greater than the max (%d)", static_cast<int>(x), static_cast<int>(m_max));
  m_set = set_add(std::move(m_set), static_cast<int>(x));
  return end();
}

template <typename Element>
typename UintSet<Element>::Iterator
UintSet<Element>::insert(
  const Iterator&,
  const Element x
)
{
  LOG_MESSAGE_IF(x > m_max, error, "Inserting a value (%d) which is greater than the max (%d)", static_cast<int>(x), static_cast<int>(m_max));
  m_set = set_add(std::move(m_set), static_cast<int>(x));
  return end();
}

template <typename Element>
void
UintSet<Element>::erase(
  const Element x
)
{
  LOG_MESSAGE_IF(!contains(x), warning, "Removing a value (%d) which does not exist in the set", static_cast<int>(x));
  m_set = set_remove(std::move(m_set), static_cast<int>(x));
}

template <typename Element>
Element
UintSet<Element>::max(
) const
{
  return m_max;
}

template <typename Element>
uint32_t
UintSet<Element>::size(
) const
{
  return static_cast<Element>(set_size(m_set));
}

template <typename Element>
bool
UintSet<Element>::empty(
) const
{
  return is_emptyset(m_set);
}

template <typename Element>
bool
UintSet<Element>::contains(
  const Element x
) const
{
  return in_set(m_set, static_cast<int>(x));
}

template <typename Element>
typename UintSet<Element>::Iterator
UintSet<Element>::begin(
) const
{
  return Iterator(&m_set, m_max, 0);
}

template <typename Element>
typename UintSet<Element>::Iterator
UintSet<Element>::end(
) const
{
  return Iterator(&m_set, m_max);
}

template <typename Element>
typename UintSet<Element>::Iterator
UintSet<Element>::find(
  const Element x
) const
{
  return contains(x) ? Iterator(&m_set, m_max, x): end();
}

#endif // DETAIL_UINTSET_HPP_
