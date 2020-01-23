/**
 * @file SetUtilsTests.hpp
 * @brief Unit tests for set util functions.
 */
#ifndef TEST_SETUTILS_HPP_
#define TEST_SETUTILS_HPP_

#include "SetUtils.hpp"
#include "UintSet.hpp"

#include <gtest/gtest.h>

#include <cstdint>
#include <numeric>
#include <set>


template <typename Set>
class SetUtilsTests: public testing::Test {
protected:
  void
  SetUp() override {
    max = UintSet<typename Set::value_type>::capacity();
  }

  typename Set::value_type max;
};

typedef testing::Types<std::set<uint8_t>, std::set<uint16_t>,
                       UintSet<uint8_t>, UintSet<uint16_t>> SetTypes;

TYPED_TEST_CASE(SetUtilsTests, SetTypes);

TYPED_TEST(SetUtilsTests, SetUnion) {
  std::vector<typename TypeParam::value_type> v(3);
  std::iota(v.begin(), v.end(), this->max - 5);
  TypeParam first(v.begin(), v.end());

  std::iota(v.begin(), v.end(), this->max - 4);
  TypeParam second(v.begin(), v.end());

  v.resize(4);
  std::iota(v.begin(), v.end(), this->max - 5);
  TypeParam expectedUnion(v.begin(), v.end());

  auto computedUnion = set_union(first, second);
  EXPECT_EQ(expectedUnion, computedUnion);
}

TYPED_TEST(SetUtilsTests, SetDifference) {
  std::vector<typename TypeParam::value_type> v(3);
  std::iota(v.begin(), v.end(), this->max - 5);
  TypeParam first(v.begin(), v.end());

  std::iota(v.begin(), v.end(), this->max - 4);
  TypeParam second(v.begin(), v.end());

  TypeParam expectedDifference;
  expectedDifference.insert(this->max - 5);

  auto computedDifference = set_difference(first, second);
  EXPECT_EQ(expectedDifference, computedDifference);
}

#endif // TEST_SETUTILS_HPP_
