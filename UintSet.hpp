/**
 * @file UintSet.hpp
 * @brief Declaration of the UintSet and related classes.
 */
#ifndef UINTSET_HPP_
#define UINTSET_HPP_

#include "bit_util.hpp"

#include <cstdint>
#include <initializer_list>
#include <vector>


/**
 * @brief Type trait class for uint_type sets.
 *        Specialized for all the supported values of N.
 *
 * @tparam N Number of 64-bit elements used for storage by the set.
 */
template <int N>
class UintTypeTrait;

/**
 * @brief Computes the N that is required for safely
 *        storing the given datatype.
 */
template <typename Element>
constexpr
int
maxN();

/**
 * @brief STL style interface for the uint_type sets from the SABNAtk library.
 *
 * @tparam Element Datatype of the value contained by the set.
 * @tparam N Number of 64-bit elements used for storage by the set.
 */
template <typename Element, int N = maxN<Element>()>
class UintSet {
  static_assert(N <= maxN<Element>(), "Provided N is bigger than the maximum required");

public:
  class Iterator;

public:
  using Set = typename UintTypeTrait<N>::Set;
  // Required typedefs
  using value_type = Element;
  using iterator = Iterator;

public:
  static
  constexpr
  Element
  capacity();

public:
  UintSet(const Element = capacity());

  UintSet(const std::initializer_list<Element>&, const Element = capacity());

  UintSet(const Iterator&, const Iterator&, const Element = capacity());

  UintSet(const typename std::vector<Element>::iterator&, const typename std::vector<Element>::iterator&, const Element = capacity());

  const Set&
  operator*() const;

  Set&
  operator*();

  bool
  operator==(const UintSet&) const;

  bool
  operator!=(const UintSet&) const;

  Iterator
  insert(const Element);

  Iterator
  insert(const Iterator&, const Element);

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

  Iterator
  begin() const;

  Iterator
  end() const;

  Iterator
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
}; // class UintSet

#include "detail/UintSet.hpp"

#endif // UINTSET_HPP_
