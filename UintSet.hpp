/**
 * @file UintSet.hpp
 * @brief Declaration of the UintSet and related classes.
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
#ifndef UINTSET_HPP_
#define UINTSET_HPP_

#include "mxx/comm.hpp"

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
  using Impl = typename UintTypeTrait<Size>::Set;
  using ReduceType = typename UintTypeTrait<Size>::ArrayType;
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

  UintSet(const Impl&, const Element = capacity());

  const Impl&
  operator*() const;

  Impl&
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

  void
  clear();

  Element
  max() const;

  Element
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

  template <typename Set>
  friend
  void
  set_bcast(Set&, const int, const mxx::comm&);

  template <typename Var, typename Set, template <typename> class Functor, typename RType>
  friend
  void
  uintset_allreduce(Set&, Functor<RType>, const mxx::comm&);

  template <typename Var, typename Set, template <typename> class Functor, typename RType>
  friend
  void
  uintset_allreduce_indexed(std::unordered_map<Var, Set>&, const Set&, const Var, Functor<RType>, const mxx::comm&);

  template <typename Set, typename Var>
  friend
  std::vector<Set>
  set_allgatherv(const std::vector<Set>&, const std::vector<size_t>&, const Var, const mxx::comm&);

private:
  Impl m_set;
  Element m_max;
  mutable Element m_size;
}; // class UintSet

#include "detail/UintSet.hpp"

#endif // UINTSET_HPP_
