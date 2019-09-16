/**
 * @file UintSetUtils.hpp
 * @brief Implementation of functions for operations on UintSet.
 */
#ifndef DETAIL_UINTSETUTILS_HPP_
#define DETAIL_UINTSETUTILS_HPP_

#include "../UintSet.hpp"


template <>
bool
set_contains<UintSet<uint8_t>, uint8_t>(
  const UintSet<uint8_t>& set,
  const uint8_t value
)
{
  return set.contains(value);
}

template <>
bool
set_contains<UintSet<uint16_t>, uint16_t>(
  const UintSet<uint16_t>& set,
  const uint16_t value
)
{
  return set.contains(value);
}

template <>
UintSet<uint8_t>
set_union<UintSet<uint8_t>>(
  const UintSet<uint8_t>& first,
  const UintSet<uint8_t>& second
)
{
  UintSet<uint8_t> result;
  result.m_set = first.m_set | second.m_set;
  return result;
}

template <>
UintSet<uint16_t>
set_union<UintSet<uint16_t>>(
  const UintSet<uint16_t>& first,
  const UintSet<uint16_t>& second
)
{
  UintSet<uint16_t> result;
  result.m_set = first.m_set | second.m_set;
  return result;
}

#endif // DETAIL_UINTSETUTILS_HPP_
