/**
 * @file Common.hpp
 * @brief Implementation of the helper classes for the unit tests.
 */
#ifndef TEST_COMMON_HPP_
#define TEST_COMMON_HPP_

#include "Data.hpp"
#include "DataFile.hpp"

#include "BVCounter.hpp"

#include <gtest/gtest.h>

#include <set>
#include <string>
#include <vector>


/**
 * @brief Helper for running unit tests on data from the examples
 *        in the Neapolitan text book.
 */
class NeapolitanTest: public testing::Test {
protected:
  void
  SetUp() override {
    uint32_t n = 2;
    uint32_t m = 8;
    data.resize(3);
    for (uint8_t i = 0; i < 3; ++i) {
      auto fileName = "neapolitan_" + std::to_string(i+1) + ".txt";
      SeparatedFile<uint8_t> dataFile(fileName, n, m, ',', false, true);
      auto bvc = create_BVCounter<1>(n, m, std::begin(dataFile.data()));
      data[i] = Data<BVCounter<1>, uint8_t>(bvc, dataFile.varNames());
    }
  }

  std::vector<Data<BVCounter<1>, uint8_t>> data;
}; // class NeapolitanTest

/**
 * @brief Helper for running unit tests on the lizards dataset.
 */
class LizardsTest: public testing::Test {
protected:
  void
  SetUp() override {
    uint32_t n = 3;
    uint32_t m = 409;
    SeparatedFile<uint8_t> dataFile("lizards.csv", n, m, ',', true, true);
    auto bvc = create_BVCounter<1>(n, m, std::begin(dataFile.data()));
    data = Data<BVCounter<1>, uint8_t>(bvc, dataFile.varNames());
  }

  Data<BVCounter<1>, uint8_t> data;
}; // class LizardsTest

/**
 * @brief Helper for running unit tests on the coronary dataset.
 */
class CoronaryTest: public testing::Test {
protected:
  void
  SetUp() override {
    uint32_t n = 6;
    uint32_t m = 1841;
    SeparatedFile<uint8_t> dataFile("coronary.csv", n, m, ',', true, true);
    auto bvc = create_BVCounter<1>(n, m, std::begin(dataFile.data()));
    data = Data<BVCounter<1>, uint8_t>(bvc, dataFile.varNames());
  }

  Data<BVCounter<1>, uint8_t> data;
}; // class CoronaryTest

/**
 * @brief Helper for running unit tests on the coronary dataset.
 */
class AsiaTest: public testing::Test {
protected:
  void
  SetUp() override {
    uint32_t n = 8;
    uint32_t m = 5000;
    SeparatedFile<uint8_t> dataFile("asia.csv", n, m, ',', true, true);
    auto bvc = create_BVCounter<1>(n, m, std::begin(dataFile.data()));
    data = Data<BVCounter<1>, uint8_t>(bvc, dataFile.varNames());
  }

  Data<BVCounter<1>, uint8_t> data;
}; // class AsiaTest

#endif // TEST_COMMON_HPP_
