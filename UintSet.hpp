/**
 * @file UintSet.hpp
 * @brief Declaration of the UintSet and related classes.
 */
#ifndef UINTSET_HPP_
#define UINTSET_HPP_

#include <cstdint>
#include <initializer_list>
#include <type_traits>
#include <vector>


/**
 * @brief Type trait class for uint_type sets.
 */
template <typename Size>
class UintTypeTrait;

/**
 * @brief Computes the number of 64-bit elements that are required
 *        for safely storing the given datatype.
 */
template <typename Element>
constexpr
int
maxSize();

/**
 * @brief STL style interface for the uint_type sets from the SABNAtk library.
 *
 * @tparam Element Datatype of the value contained by the set.
 * @tparam Size Number of 64-bit elements used for storage by the set.
 */
template <typename Element, typename Size = std::integral_constant<int, maxSize<Element>()>>
class UintSet {
  static_assert(Size::value <= maxSize<Element>(), "Provided size is bigger than the maximum required");

public:
  class Enumerator;

public:
  using Set = typename UintTypeTrait<Size>::Set;
  // Required typedefs
  using value_type = Element;
  using iterator = Enumerator;

public:
  static
  constexpr
  Element
  capacity();

public:
  UintSet(const Element = capacity());

  UintSet(const std::initializer_list<Element>&, const Element = capacity());

  UintSet(const iterator&, const iterator&, const Element = capacity());

  UintSet(const typename std::vector<Element>::iterator&, const typename std::vector<Element>::iterator&, const Element = capacity());

  const Set&
  operator*() const;

  Set&
  operator*();

  bool
  operator==(const UintSet&) const;

  bool
  operator!=(const UintSet&) const;

  iterator
  insert(const Element);

  iterator
  insert(const iterator&, const Element);

  void
  erase(const Element);

  Element
  max() const;

  uint32_t
  size() const;

  bool
  empty() const;

  bool
  contains(const Element) const;

  iterator
  begin() const;

  iterator
  end() const;

  iterator
  find(const Element) const;

  UintSet
  subset(const std::vector<bool>&) const;

public:
  template <typename Set, typename Var>
  friend
  Set
  set_init(Set&&, const Var);

private:
  Set m_set;
  Element m_max;
  mutable Element m_size;
}; // class UintSet

#include "detail/UintSet.hpp"

#endif // UINTSET_HPP_
