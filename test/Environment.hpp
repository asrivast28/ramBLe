/**
 * @file Environment.hpp
 * @brief Implementation of the environment classes for the unit tests.
 */
#ifndef TEST_ENVIRONMENT_HPP_
#define TEST_ENVIRONMENT_HPP_

#include "DataReader.hpp"

#include <gtest/gtest.h>

#include <string>
#include <vector>


/**
 * @brief Environment for running unit tests on data from the examples
 *        in the Neapolitan text book.
 */
class NeapolitanEnvironment: public testing::Environment {
public:
  void
  SetUp() override {
    n = 2;
    m = 8;
    reader.resize(3);
    for (uint8_t i = 0; i < 3; ++i) {
      auto fileName = "neapolitan_" + std::to_string(i+1) + ".txt";
      reader[i] = new RowObservationReader<uint8_t>(fileName, n, m, ',', false, false, true);
    }
  }

  void
  TearDown() override {
    for (auto i = 0u; i < 3u; ++i) {
      delete reader[i];
    }
  }

  static std::vector<RowObservationReader<uint8_t>*> reader;
  static uint32_t n;
  static uint32_t m;
}; // class NeapolitanEnvironment

// Static member initializations
std::vector<RowObservationReader<uint8_t>*> NeapolitanEnvironment::reader;
uint32_t NeapolitanEnvironment::n = 0;
uint32_t NeapolitanEnvironment::m = 0;


/**
 * @brief Environment for running unit tests on the lizards dataset.
 */
class LizardsEnvironment: public testing::Environment {
public:
  void
  SetUp() override {
    n = 3;
    m = 409;
    reader = new RowObservationReader<uint8_t>("lizards.csv", n, m, ',', true, false, true);
  }

  void
  TearDown() override {
    delete reader;
  }

  static RowObservationReader<uint8_t>* reader;
  static uint32_t n;
  static uint32_t m;
}; // class LizardsEnvironment

// Static member initializations
RowObservationReader<uint8_t>* LizardsEnvironment::reader;
uint32_t LizardsEnvironment::n = 0;
uint32_t LizardsEnvironment::m = 0;


/**
 * @brief Environment for running unit tests on the coronary dataset.
 */
class CoronaryEnvironment: public testing::Environment {
public:
  void
  SetUp() override {
    n = 6;
    m = 1841;
    reader = new RowObservationReader<uint8_t>("coronary.csv", n, m, ',', true, false, true);
  }

  void
  TearDown() override {
    delete reader;
  }

  static RowObservationReader<uint8_t>* reader;
  static uint32_t n;
  static uint32_t m;
}; // class CoronaryEnvironment

// Static member initializations
RowObservationReader<uint8_t>* CoronaryEnvironment::reader;
uint32_t CoronaryEnvironment::n = 0;
uint32_t CoronaryEnvironment::m = 0;


/**
 * @brief Environment for running unit tests on the coronary dataset.
 */
class AsiaEnvironment: public testing::Environment {
public:
  void
  SetUp() override {
    n = 8;
    m = 5000;
    reader = new RowObservationReader<uint8_t>("asia.csv", n, m, ',', true, false, true);
  }

  void
  TearDown() override {
    delete reader;
  }

  static RowObservationReader<uint8_t>* reader;
  static uint32_t n;
  static uint32_t m;
}; // class AsiaEnvironment

// Static member initializations
RowObservationReader<uint8_t>* AsiaEnvironment::reader;
uint32_t AsiaEnvironment::n = 0;
uint32_t AsiaEnvironment::m = 0;

#endif // TEST_ENVIRONMENT_HPP_
