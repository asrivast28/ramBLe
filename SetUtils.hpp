/**
 * @file SetUtils.hpp
 * @brief Declaration of functions for common set operations.
 */
#ifndef SETUTILS_HPP_
#define SETUTILS_HPP_


/**
 * @brief Function for initializing a given set.
 */
template <typename SetType, typename ValueType>
SetType
set_init(SetType&&, const ValueType);

/**
 * @brief Function for checking if a given set contains a value.
 */
template <typename SetType, typename ValueType>
bool
set_contains(const SetType&, const ValueType);

/**
 * @brief Function for getting the union of two given sets.
 */
template <typename SetType>
SetType
set_union(const SetType&, const SetType&);

#include "detail/StdSetUtils.hpp"
#include "detail/UintSetUtils.hpp"

#endif // SETUTILS_HPP_
