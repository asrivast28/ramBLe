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
  using SetType = uint_type<N>;

public:
  static
  constexpr
  uint32_t
  max()
  {
    return static_cast<uint32_t>(set_max_size<SetType>());
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
  using SetType = uint_type<N>;

public:
  static
  constexpr
  uint32_t
  max()
  {
    return static_cast<uint32_t>(set_max_size<SetType>());
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
  using SetType = uint_type<N>;

public:
  static
  constexpr
  uint32_t
  max()
  {
    return static_cast<uint32_t>(set_max_size<SetType>());
  }
};


template <typename ValueType>
/**
 * @brief Iterator for UintSet.
 */
class UintSet<ValueType>::Iterator : public std::iterator<std::forward_iterator_tag, ValueType> {
public:
  using SetType = typename UintTypeTrait<ValueType>::SetType;

public:
  Iterator(const SetType* set, const ValueType max, const ValueType curr = UintTypeTrait<ValueType>::max())
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

  ValueType
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
  const SetType* m_set;
  ValueType m_max;
  ValueType m_curr;
};


template <typename ValueType>
/**
 * @brief The capacity of a set of type UintSet<ValueType>.
 */
constexpr
uint32_t
UintSet<ValueType>::capacity(
)
{
  return UintTypeTrait<ValueType>::max();
}

template <typename ValueType>
UintSet<ValueType>::UintSet(
  const ValueType max
) : m_set(set_empty<SetType>()),
    m_max(max)
{
}

template <typename ValueType>
UintSet<ValueType>::UintSet(
  const std::initializer_list<ValueType>& s,
  const ValueType max
) : m_set(as_set<SetType>(s.begin(), s.end())),
    m_max(max)
{
}

template <typename ValueType>
UintSet<ValueType>::UintSet(
  const typename UintSet<ValueType>::Iterator& first,
  const typename UintSet<ValueType>::Iterator& last,
  const ValueType max
) : m_set(as_set<SetType>(first, last)),
    m_max(max)
{
}

template <typename ValueType>
UintSet<ValueType>::UintSet(
  const typename std::vector<ValueType>::iterator& first,
  const typename std::vector<ValueType>::iterator& last,
  const ValueType max
) : m_set(as_set<SetType>(first, last)),
    m_max(max)
{
}

template <typename ValueType>
UintSet<ValueType>::UintSet(
  const std::vector<bool>& bitset,
  const ValueType max
) : m_set(set_empty<SetType>()),
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

template <typename ValueType>
const typename UintSet<ValueType>::SetType&
UintSet<ValueType>::operator*(
) const
{
  return m_set;
}

template <typename ValueType>
typename UintSet<ValueType>::SetType&
UintSet<ValueType>::operator*(
)
{
  return m_set;
}

template <typename ValueType>
bool
UintSet<ValueType>::operator==(
  const UintSet<ValueType>& other
) const
{
  return (m_set == other.m_set);
}

template <typename ValueType>
bool
UintSet<ValueType>::operator!=(
  const UintSet<ValueType>& other
) const
{
  return (m_set != other.m_set);
}

template <typename ValueType>
typename UintSet<ValueType>::Iterator
UintSet<ValueType>::insert(
  const ValueType x
)
{
  LOG_MESSAGE_IF(x >= m_max, error, "Inserting a value (%d) which is greater than the max (%d)", static_cast<int>(x), static_cast<int>(m_max));
  m_set = set_add(std::move(m_set), static_cast<int>(x));
  return end();
}

template <typename ValueType>
typename UintSet<ValueType>::Iterator
UintSet<ValueType>::insert(
  const Iterator&,
  const ValueType x
)
{
  LOG_MESSAGE_IF(x >= m_max, error, "Inserting a value (%d) which is greater than the max (%d)", static_cast<int>(x), static_cast<int>(m_max));
  m_set = set_add(std::move(m_set), static_cast<int>(x));
  return end();
}

template <typename ValueType>
void
UintSet<ValueType>::erase(
  const ValueType x
)
{
  LOG_MESSAGE_IF(!contains(x), warning, "Removing a value (%d) which does not exist in the set", static_cast<int>(x));
  m_set = set_remove(std::move(m_set), static_cast<int>(x));
}

template <typename ValueType>
ValueType
UintSet<ValueType>::max(
) const
{
  return m_max;
}

template <typename ValueType>
uint32_t
UintSet<ValueType>::size(
) const
{
  return static_cast<ValueType>(set_size(m_set));
}

template <typename ValueType>
bool
UintSet<ValueType>::empty(
) const
{
  return is_emptyset(m_set);
}

template <typename ValueType>
bool
UintSet<ValueType>::contains(
  const ValueType x
) const
{
  return in_set(m_set, static_cast<int>(x));
}

template <typename ValueType>
typename UintSet<ValueType>::Iterator
UintSet<ValueType>::begin(
) const
{
  return Iterator(&m_set, m_max, 0);
}

template <typename ValueType>
typename UintSet<ValueType>::Iterator
UintSet<ValueType>::end(
) const
{
  return Iterator(&m_set, m_max);
}

template <typename ValueType>
typename UintSet<ValueType>::Iterator
UintSet<ValueType>::find(
  const ValueType x
) const
{
  return contains(x) ? Iterator(&m_set, m_max, x): end();
}

#endif // DETAIL_UINTSET_HPP_
