/**
 * @file NetworkTests.hpp
 * @brief Unit tests for testing learning of directed networks.
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
#ifndef TEST_NETWORKTESTS_HPP_
#define TEST_NETWORKTESTS_HPP_

#include "BlanketLearning.hpp"
#include "DirectLearning.hpp"
#include "NetworkData.hpp"


// Declaration of type-parameterized tests
TYPED_TEST_CASE_P(ChildData);
TYPED_TEST_CASE_P(InsuranceData);
TYPED_TEST_CASE_P(MildewData);
TYPED_TEST_CASE_P(AlarmData);

// Type-parameterized sequential tests
TYPED_TEST_P(ChildData, Sequential) {
  auto sequentialBN = this->algo->getNetwork(true, false);
  EXPECT_EQ(*(this->expected), sequentialBN);
}

TYPED_TEST_P(InsuranceData, Sequential) {
  auto sequentialBN = this->algo->getNetwork(true, false);
  EXPECT_EQ(*(this->expected), sequentialBN);
}

TYPED_TEST_P(MildewData, Sequential) {
  auto sequentialBN = this->algo->getNetwork(true, false);
  EXPECT_EQ(*(this->expected), sequentialBN);
}

TYPED_TEST_P(AlarmData, Sequential) {
  auto sequentialBN = this->algo->getNetwork(true, false);
  EXPECT_EQ(*(this->expected), sequentialBN);
}

// Type-parameterized parallel tests
TYPED_TEST_P(ChildData, Parallel) {
  try {
    auto parallelBN = this->algo->getNetwork(true, true);
    EXPECT_EQ(*(this->expected), parallelBN);
  }
  catch (const NotImplementedError& e) {
    if (this->comm.is_first()) {
#ifdef GTEST_SKIP
      GTEST_SKIP() << e.what();
#else
      std::cerr << "[  SKIPPED ] " << e.what() << std::endl;
#endif
    }
  }
}

TYPED_TEST_P(InsuranceData, Parallel) {
  try {
    auto parallelBN = this->algo->getNetwork(true, true);
    EXPECT_EQ(*(this->expected), parallelBN);
  }
  catch (const NotImplementedError& e) {
    if (this->comm.is_first()) {
#ifdef GTEST_SKIP
      GTEST_SKIP() << e.what();
#else
      std::cerr << "[  SKIPPED ] " << e.what() << std::endl;
#endif
    }
  }
}

TYPED_TEST_P(MildewData, Parallel) {
  try {
    auto parallelBN = this->algo->getNetwork(true, true);
    EXPECT_EQ(*(this->expected), parallelBN);
  }
  catch (const NotImplementedError& e) {
    if (this->comm.is_first()) {
#ifdef GTEST_SKIP
      GTEST_SKIP() << e.what();
#else
      std::cerr << "[  SKIPPED ] " << e.what() << std::endl;
#endif
    }
  }
}

TYPED_TEST_P(AlarmData, Parallel) {
  try {
    auto parallelBN = this->algo->getNetwork(true, true);
    EXPECT_EQ(*(this->expected), parallelBN);
  }
  catch (const NotImplementedError& e) {
    if (this->comm.is_first()) {
#ifdef GTEST_SKIP
      GTEST_SKIP() << e.what();
#else
      std::cerr << "[  SKIPPED ] " << e.what() << std::endl;
#endif
    }
  }
}

// Register sequential as well as parallel tests
REGISTER_TYPED_TEST_CASE_P(ChildData, Sequential, Parallel);
REGISTER_TYPED_TEST_CASE_P(InsuranceData, Sequential, Parallel);
REGISTER_TYPED_TEST_CASE_P(MildewData, Sequential, Parallel);
REGISTER_TYPED_TEST_CASE_P(AlarmData, Sequential, Parallel);

// All the different learning algorithms
using LearningAlgorithms = testing::Types<GS<DiscreteData<Counter, uint8_t>, uint8_t, UintSet<uint8_t>>,
                                          IAMB<DiscreteData<Counter, uint8_t>, uint8_t, UintSet<uint8_t>>,
                                          InterIAMB<DiscreteData<Counter, uint8_t>, uint8_t, UintSet<uint8_t>>,
                                          MMPC<DiscreteData<Counter, uint8_t>, uint8_t, UintSet<uint8_t>>,
                                          HITON<DiscreteData<Counter, uint8_t>, uint8_t, UintSet<uint8_t>>,
                                          SemiInterleavedHITON<DiscreteData<Counter, uint8_t>, uint8_t, UintSet<uint8_t>>,
                                          GetPC<DiscreteData<Counter, uint8_t>, uint8_t, UintSet<uint8_t>>>;

// Instantiate all the type-parameterized tests for the learning algorithms
INSTANTIATE_TYPED_TEST_CASE_P(Network, ChildData, LearningAlgorithms);
INSTANTIATE_TYPED_TEST_CASE_P(Network, InsuranceData, LearningAlgorithms);
INSTANTIATE_TYPED_TEST_CASE_P(Network, MildewData, LearningAlgorithms);
INSTANTIATE_TYPED_TEST_CASE_P(Network, AlarmData, LearningAlgorithms);

#endif // TEST_NETWORKTESTS_HPP_
