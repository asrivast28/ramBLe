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
 * @tparam ValueType Datatype of the value contained by the set.
 */
template <typename ValueType>
class UintTypeTrait;

/**
 * @brief STL style interface for the uint_type sets from the SABNAtk library.
 *
 * @tparam ValueType Datatype of the value contained by the set.
 */
template <typename ValueType>
class UintSet {
public:
  class Iterator;

public:
  using SetType = typename UintTypeTrait<ValueType>::SetType;
  // Required typedefs
  using value_type = ValueType;
  using iterator = Iterator;

public:
  static
  constexpr
  uint32_t
  capacity();

public:
  UintSet(const ValueType = capacity());

  UintSet(const std::initializer_list<ValueType>&, const ValueType = capacity());

  UintSet(const typename std::vector<ValueType>::iterator&, const typename std::vector<ValueType>::iterator&, const ValueType = capacity());

  bool
  operator==(const UintSet&) const;

  bool
  operator!=(const UintSet&) const;

  Iterator
  insert(const ValueType);

  Iterator
  insert(const Iterator&, const ValueType);

  void
  erase(const ValueType);

  uint32_t
  size() const;

  bool
  empty() const;

  bool
  contains(const ValueType) const;

  Iterator
  begin() const;

  Iterator
  end() const;

  Iterator
  find(const ValueType) const;

public:
  template <typename SetType, typename VarType>
  friend
  SetType
  set_init(SetType&&, const VarType);

  template <typename SetType>
  friend
  SetType
  set_union(const SetType&, const SetType&);

private:
  SetType m_set;
  ValueType m_max;
}; // class UintSet

#include "detail/UintSet.hpp"

#endif // UINTSET_HPP_
