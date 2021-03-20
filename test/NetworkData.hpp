/**
 * @file NetworkData.hpp
 * @brief Data classes for testing network learning.
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
#ifndef TEST_NETWORKDATA_HPP_
#define TEST_NETWORKDATA_HPP_

#include "DiscreteData.hpp"

#include "common/CTCounter.hpp"


using Counter = CTCounter<>;

template <typename Algorithm>
class ChildData : public testing::Test {
protected:
  void
  SetUp() override {
    constexpr uint32_t n = 20;
    constexpr uint32_t m = 10000;
    ColumnObservationReader<uint8_t> reader("child.csv", n, m, ' ', false, false, true);
    auto counter = Counter::create(n, m, std::begin(reader.data()));
    data = new DiscreteData<Counter, uint8_t>(counter, reader.varNames());
    algo = new Algorithm(comm, *data);
    expected = this->expectedBN();
  }

  BayesianNetwork<uint8_t>*
  expectedBN();

  void
  TearDown() override {
    delete expected;
    delete algo;
    delete data;
  }

  mxx::comm comm;
  DiscreteData<Counter, uint8_t>* data;
  Algorithm* algo;
  BayesianNetwork<uint8_t>* expected;
};


template <typename Algorithm>
class InsuranceData : public testing::Test {
protected:
  void
  SetUp() override {
    constexpr uint32_t n = 27;
    constexpr uint32_t m = 10000;
    ColumnObservationReader<uint8_t> reader("insurance.csv", n, m, ' ', false, false, true);
    auto counter = Counter::create(n, m, std::begin(reader.data()));
    data = new DiscreteData<Counter, uint8_t>(counter, reader.varNames());
    algo = new Algorithm(comm, *data);
    expected = this->expectedBN();
  }

  BayesianNetwork<uint8_t>*
  expectedBN();

  void
  TearDown() override {
    delete expected;
    delete algo;
    delete data;
  }

  mxx::comm comm;
  DiscreteData<Counter, uint8_t>* data;
  Algorithm* algo;
  BayesianNetwork<uint8_t>* expected;
};


template <typename Algorithm>
class MildewData : public testing::Test {
protected:
  void
  SetUp() override {
    constexpr uint32_t n = 35;
    constexpr uint32_t m = 10000;
    ColumnObservationReader<uint8_t> reader("mildew.csv", n, m, ' ', false, false, true);
    auto counter = Counter::create(n, m, std::begin(reader.data()));
    data = new DiscreteData<Counter, uint8_t>(counter, reader.varNames());
    algo = new Algorithm(comm, *data);
    expected = this->expectedBN();
  }

  BayesianNetwork<uint8_t>*
  expectedBN();

  void
  TearDown() override {
    delete expected;
    delete algo;
    delete data;
  }

  mxx::comm comm;
  DiscreteData<Counter, uint8_t>* data;
  Algorithm* algo;
  BayesianNetwork<uint8_t>* expected;
};


template <typename Algorithm>
class AlarmData : public testing::Test {
protected:
  void
  SetUp() override {
    constexpr uint32_t n = 37;
    constexpr uint32_t m = 10000;
    ColumnObservationReader<uint8_t> reader("alarm.csv", n, m, ' ', false, false, true);
    auto counter = Counter::create(n, m, std::begin(reader.data()));
    data = new DiscreteData<Counter, uint8_t>(counter, reader.varNames());
    algo = new Algorithm(comm, *data);
    expected = this->expectedBN();
  }

  BayesianNetwork<uint8_t>*
  expectedBN();

  void
  TearDown() override {
    delete expected;
    delete algo;
    delete data;
  }

  mxx::comm comm;
  DiscreteData<Counter, uint8_t>* data;
  Algorithm* algo;
  BayesianNetwork<uint8_t>* expected;
};


#include "GSNetworks.hpp"
#include "IAMBNetworks.hpp"
#include "InterIAMBNetworks.hpp"
#include "MMPCNetworks.hpp"
#include "HITONNetworks.hpp"
#include "SemiInterleavedHITONNetworks.hpp"
#include "GetPCNetworks.hpp"
#include "PCStableNetworks.hpp"

#endif // TEST_NETWORKDATA_HPP_
