/**
 * @file SetUtilsTests.hpp
 * @brief Unit tests for set util functions.
 */
#ifndef TEST_SETUTILS_HPP_
#define TEST_SETUTILS_HPP_

#include "SetUtils.hpp"
#include "UintSet.hpp"

#include <cstdint>
#include <set>


TEST(SetUtilsTests, StdSetContains) {
  std::set<uint8_t> set_8{1, 2, 3, 4, 5, 6};
  EXPECT_TRUE(set_contains(set_8, static_cast<uint8_t>(3)));
  EXPECT_FALSE(set_contains(set_8, static_cast<uint8_t>(7)));

  std::set<uint16_t> set_16{64, 65, 66, 67, 68, 69};
  EXPECT_TRUE(set_contains(set_16, static_cast<uint16_t>(65)));
  EXPECT_FALSE(set_contains(set_16, static_cast<uint16_t>(70)));
}

TEST(SetUtilsTests, UintSetContains) {
  UintSet<uint8_t> set_8{1, 2, 3, 4, 5, 6};
  EXPECT_TRUE(set_contains(set_8, static_cast<uint8_t>(3)));
  EXPECT_FALSE(set_contains(set_8, static_cast<uint8_t>(7)));

  UintSet<uint16_t> set_16{64, 65, 66, 67, 68, 69};
  EXPECT_TRUE(set_contains(set_16, static_cast<uint16_t>(65)));
  EXPECT_FALSE(set_contains(set_16, static_cast<uint16_t>(70)));
}

TEST(SetUtilsTests, StdSetUnion) {
  std::set<uint8_t> first_8{1, 2, 3};
  std::set<uint8_t> second_8{2, 3, 4};
  std::set<uint8_t> expectedUnion_8{1, 2, 3, 4};

  auto computedUnion_8 = set_union(first_8, second_8);
  EXPECT_EQ(expectedUnion_8, computedUnion_8);

  std::set<uint16_t> first_16{64, 65, 66};
  std::set<uint16_t> second_16{65, 66, 67};
  std::set<uint16_t> expectedUnion_16{64, 65, 66, 67};

  auto computedUnion_16 = set_union(first_16, second_16);
  EXPECT_EQ(expectedUnion_16, computedUnion_16);
}

TEST(SetUtilsTests, UintSetUnion) {
  UintSet<uint8_t> first_8{1, 2, 3};
  UintSet<uint8_t> second_8{2, 3, 4};
  UintSet<uint8_t> expectedUnion_8{1, 2, 3, 4};

  auto computedUnion_8 = set_union(first_8, second_8);
  EXPECT_EQ(expectedUnion_8, computedUnion_8);


  UintSet<uint16_t> first_16{64, 65, 66};
  UintSet<uint16_t> second_16{65, 66, 67};
  UintSet<uint16_t> expectedUnion_16{64, 65, 66, 67};

  auto computedUnion_16 = set_union(first_16, second_16);
  EXPECT_EQ(expectedUnion_16, computedUnion_16);
}

#endif // TEST_SETUTILS_HPP_
