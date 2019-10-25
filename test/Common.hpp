/**
 * @file Common.hpp
 * @brief Implementation of the helper classes for the unit tests.
 */
#ifndef TEST_COMMON_HPP_
#define TEST_COMMON_HPP_

#include "DataQuery.hpp"
#include "DataReader.hpp"

#include "BVCounter.hpp"
#include "CTCounter.hpp"
#include "RadCounter.hpp"

#include <gtest/gtest.h>

#include <set>
#include <string>
#include <vector>


// All the different counter implementations
typedef testing::Types<BVCounter<1>, CTCounter<1>, RadCounter<1>> AllCounters;
// Default counter implementation
typedef testing::Types<BVCounter<1>> DefaultCounter;

/**
 * @brief Helper for running unit tests on data from the examples
 *        in the Neapolitan text book.
 */
template <typename Counter>
class NeapolitanTest: public testing::Test {
protected:
  void
  SetUp() override {
    uint32_t n = 2;
    uint32_t m = 8;
    data.resize(3);
    for (uint8_t i = 0; i < 3; ++i) {
      auto fileName = "neapolitan_" + std::to_string(i+1) + ".txt";
      RowObservationReader<uint8_t> dataFile(fileName, n, m, ',', false, false, true);
      auto counter = Counter::create(n, m, std::begin(dataFile.data()));
      data[i] = DataQuery<Counter, uint8_t>(counter, dataFile.varNames());
    }
  }

  std::vector<DataQuery<Counter, uint8_t>> data;
}; // class NeapolitanTest

TYPED_TEST_CASE(NeapolitanTest, AllCounters);

/**
 * @brief Helper for running unit tests on the lizards dataset.
 */
template <typename Counter>
class LizardsTest: public testing::Test {
protected:
  void
  SetUp() override {
    uint32_t n = 3;
    uint32_t m = 409;
    RowObservationReader<uint8_t> dataFile("lizards.csv", n, m, ',', true, false, true);
    auto counter = Counter::create(n, m, std::begin(dataFile.data()));
    data = DataQuery<Counter, uint8_t>(counter, dataFile.varNames());
  }

  DataQuery<Counter, uint8_t> data;
}; // class LizardsTest

TYPED_TEST_CASE(LizardsTest, AllCounters);

/**
 * @brief Helper for running unit tests on the coronary dataset.
 */
template <typename Counter>
class CoronaryTest: public testing::Test {
protected:
  void
  SetUp() override {
    uint32_t n = 6;
    uint32_t m = 1841;
    RowObservationReader<uint8_t> dataFile("coronary.csv", n, m, ',', true, false, true);
    auto counter = Counter::create(n, m, std::begin(dataFile.data()));
    data = DataQuery<Counter, uint8_t>(counter, dataFile.varNames());
  }

  DataQuery<Counter, uint8_t> data;
}; // class CoronaryTest

TYPED_TEST_CASE(CoronaryTest, DefaultCounter);

/**
 * @brief Helper for running unit tests on the coronary dataset.
 */
template <typename Counter>
class AsiaTest: public testing::Test {
protected:
  void
  SetUp() override {
    uint32_t n = 8;
    uint32_t m = 5000;
    RowObservationReader<uint8_t> dataFile("asia.csv", n, m, ',', true, false, true);
    auto counter = Counter::create(n, m, std::begin(dataFile.data()));
    data = DataQuery<Counter, uint8_t>(counter, dataFile.varNames());
  }

  DataQuery<Counter, uint8_t> data;
}; // class AsiaTest

TYPED_TEST_CASE(AsiaTest, DefaultCounter);

#endif // TEST_COMMON_HPP_
