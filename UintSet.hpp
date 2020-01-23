/**
 * @file UintSet.hpp
 * @brief Declaration of the UintSet and related classes.
 */
#ifndef UINTSET_HPP_
#define UINTSET_HPP_

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
  class Enumerator;

public:
  using Set = typename UintTypeTrait<N>::Set;
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
}; // class UintSet

#include "detail/UintSet.hpp"

#endif // UINTSET_HPP_
