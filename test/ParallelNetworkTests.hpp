/**
 * @file ParallelNetworkTests.hpp
 * @brief Unit tests for testing directed network in parallel.
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
#ifndef TEST_PARALLELNETWORK_HPP_
#define TEST_PARALLELNETWORK_HPP_

#include "BlanketLearning.hpp"
#include "DirectLearning.hpp"
#include "NetworkData.hpp"


// All the different parallel learning algorithms
using ParallelLearningAlgorithms = testing::Types<GS<DiscreteData<Counter, uint8_t>, uint8_t, UintSet<uint8_t>>,
                                                  IAMB<DiscreteData<Counter, uint8_t>, uint8_t, UintSet<uint8_t>>,
                                                  InterIAMB<DiscreteData<Counter, uint8_t>, uint8_t, UintSet<uint8_t>>,
                                                  MMPC<DiscreteData<Counter, uint8_t>, uint8_t, UintSet<uint8_t>>,
                                                  SemiInterleavedHITON<DiscreteData<Counter, uint8_t>, uint8_t, UintSet<uint8_t>>>;


TYPED_TEST_CASE(ChildData, ParallelLearningAlgorithms);

TYPED_TEST(ChildData, ParallelNetwork) {
  auto sequentialBN = this->algo->getNetwork(true, false);
  auto parallelBN = this->algo->getNetwork(true, true);
  if (this->comm.is_first()) {
    EXPECT_EQ(sequentialBN, parallelBN);
  }
}


TYPED_TEST_CASE(InsuranceData, ParallelLearningAlgorithms);

TYPED_TEST(InsuranceData, ParallelNetwork) {
  auto sequentialBN = this->algo->getNetwork(true, false);
  auto parallelBN = this->algo->getNetwork(true, true);
  if (this->comm.is_first()) {
    EXPECT_EQ(sequentialBN, parallelBN);
  }
}


TYPED_TEST_CASE(MildewData, ParallelLearningAlgorithms);

TYPED_TEST(MildewData, ParallelNetwork) {
  auto sequentialBN = this->algo->getNetwork(true, false);
  auto parallelBN = this->algo->getNetwork(true, true);
  if (this->comm.is_first()) {
    EXPECT_EQ(sequentialBN, parallelBN);
  }
}


TYPED_TEST_CASE(AlarmData, ParallelLearningAlgorithms);

TYPED_TEST(AlarmData, ParallelNetwork) {
  auto sequentialBN = this->algo->getNetwork(true, false);
  auto parallelBN = this->algo->getNetwork(true, true);
  if (this->comm.is_first()) {
    EXPECT_EQ(sequentialBN, parallelBN);
  }
}

#endif // TEST_PARALLELNETWORK_HPP_
