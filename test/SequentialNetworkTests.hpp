/**
 * @file SequentialNetworkTests.hpp
 * @brief Unit tests for learning network sequentially.
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
#ifndef TEST_SEQUENTIALNETWORK_HPP_
#define TEST_SEQUENTIALNETWORK_HPP_

#include "BlanketLearning.hpp"
#include "DirectLearning.hpp"
#include "NetworkData.hpp"


class ChildData_GS : public ChildData<GS<DiscreteData<Counter, uint8_t>, uint8_t, UintSet<uint8_t>>>,
                     public testing::WithParamInterface<bool> {
};

TEST_P(ChildData_GS, DirectedNetwork) {
  auto computedBN = this->algo->getNetwork(true, GetParam());

  auto expectedBN = BayesianNetwork<uint8_t>(this->data->varNames());
  expectedBN.addEdge(static_cast<uint8_t>(1), static_cast<uint8_t>(2));
  expectedBN.addEdge(static_cast<uint8_t>(1), static_cast<uint8_t>(7));
  expectedBN.addEdge(static_cast<uint8_t>(2), static_cast<uint8_t>(7));
  expectedBN.addEdge(static_cast<uint8_t>(2), static_cast<uint8_t>(8));
  expectedBN.addEdge(static_cast<uint8_t>(3), static_cast<uint8_t>(4));
  expectedBN.addEdge(static_cast<uint8_t>(3), static_cast<uint8_t>(9));
  expectedBN.addEdge(static_cast<uint8_t>(4), static_cast<uint8_t>(2));
  expectedBN.addEdge(static_cast<uint8_t>(4), static_cast<uint8_t>(3));
  expectedBN.addEdge(static_cast<uint8_t>(4), static_cast<uint8_t>(10));
  expectedBN.addEdge(static_cast<uint8_t>(5), static_cast<uint8_t>(12));
  expectedBN.addEdge(static_cast<uint8_t>(5), static_cast<uint8_t>(17));
  expectedBN.addEdge(static_cast<uint8_t>(6), static_cast<uint8_t>(14));
  expectedBN.addEdge(static_cast<uint8_t>(9), static_cast<uint8_t>(3));
  expectedBN.addEdge(static_cast<uint8_t>(10), static_cast<uint8_t>(4));
  expectedBN.addEdge(static_cast<uint8_t>(11), static_cast<uint8_t>(15));
  expectedBN.addEdge(static_cast<uint8_t>(12), static_cast<uint8_t>(5));
  expectedBN.addEdge(static_cast<uint8_t>(14), static_cast<uint8_t>(6));
  expectedBN.addEdge(static_cast<uint8_t>(15), static_cast<uint8_t>(11));
  // Check false positives
  EXPECT_NE(expectedBN, computedBN);

  expectedBN.addEdge(static_cast<uint8_t>(17), static_cast<uint8_t>(5));
  EXPECT_EQ(expectedBN, computedBN);
}

INSTANTIATE_TEST_CASE_P(Sequential, ChildData_GS, testing::Values(false));
INSTANTIATE_TEST_CASE_P(Parallel, ChildData_GS, testing::Values(true));

class ChildData_MMPC : public ChildData<MMPC<DiscreteData<Counter, uint8_t>, uint8_t, UintSet<uint8_t>>>,
                       public testing::WithParamInterface<bool> {
};

TEST_P(ChildData_MMPC, DirectedNetwork) {
  auto computedBN = this->algo->getNetwork(true, GetParam());

  auto expectedBN = BayesianNetwork<uint8_t>(this->data->varNames());
  expectedBN.addEdge(static_cast<uint8_t>(1), static_cast<uint8_t>(7));
  expectedBN.addEdge(static_cast<uint8_t>(15), static_cast<uint8_t>(1));
  expectedBN.addEdge(static_cast<uint8_t>(16), static_cast<uint8_t>(1));
  expectedBN.addEdge(static_cast<uint8_t>(2), static_cast<uint8_t>(7));
  expectedBN.addEdge(static_cast<uint8_t>(2), static_cast<uint8_t>(8));
  expectedBN.addEdge(static_cast<uint8_t>(2), static_cast<uint8_t>(16));
  expectedBN.addEdge(static_cast<uint8_t>(3), static_cast<uint8_t>(9));
  expectedBN.addEdge(static_cast<uint8_t>(3), static_cast<uint8_t>(17));
  expectedBN.addEdge(static_cast<uint8_t>(4), static_cast<uint8_t>(10));
  expectedBN.addEdge(static_cast<uint8_t>(5), static_cast<uint8_t>(12));
  expectedBN.addEdge(static_cast<uint8_t>(6), static_cast<uint8_t>(14));
  expectedBN.addEdge(static_cast<uint8_t>(8), static_cast<uint8_t>(2));
  expectedBN.addEdge(static_cast<uint8_t>(9), static_cast<uint8_t>(3));
  expectedBN.addEdge(static_cast<uint8_t>(11), static_cast<uint8_t>(15));
  expectedBN.addEdge(static_cast<uint8_t>(11), static_cast<uint8_t>(16));
  expectedBN.addEdge(static_cast<uint8_t>(14), static_cast<uint8_t>(6));
  expectedBN.addEdge(static_cast<uint8_t>(15), static_cast<uint8_t>(11));
  expectedBN.addEdge(static_cast<uint8_t>(17), static_cast<uint8_t>(3));
  expectedBN.addEdge(static_cast<uint8_t>(17), static_cast<uint8_t>(4));
  expectedBN.addEdge(static_cast<uint8_t>(18), static_cast<uint8_t>(4));
  expectedBN.addEdge(static_cast<uint8_t>(17), static_cast<uint8_t>(5));
  expectedBN.addEdge(static_cast<uint8_t>(19), static_cast<uint8_t>(5));
  expectedBN.addEdge(static_cast<uint8_t>(13), static_cast<uint8_t>(19));
  // Check false positives
  EXPECT_NE(expectedBN, computedBN);

  expectedBN.addEdge(static_cast<uint8_t>(19), static_cast<uint8_t>(13));
  EXPECT_EQ(expectedBN, computedBN);
}

INSTANTIATE_TEST_CASE_P(Sequential, ChildData_MMPC, testing::Values(false));
INSTANTIATE_TEST_CASE_P(Parallel, ChildData_MMPC, testing::Values(true));


class InsuranceData_GS : public InsuranceData<GS<DiscreteData<Counter, uint8_t>, uint8_t, UintSet<uint8_t>>>,
                         public testing::WithParamInterface<bool> {
};

TEST_P(InsuranceData_GS, DirectedNetwork) {
  auto computedBN = this->algo->getNetwork(true, GetParam());

  auto expectedBN = BayesianNetwork<uint8_t>(this->data->varNames());
  expectedBN.addEdge(static_cast<uint8_t>(0), static_cast<uint8_t>(1));
  expectedBN.addEdge(static_cast<uint8_t>(1), static_cast<uint8_t>(0));
  expectedBN.addEdge(static_cast<uint8_t>(1), static_cast<uint8_t>(13));
  expectedBN.addEdge(static_cast<uint8_t>(2), static_cast<uint8_t>(4));
  expectedBN.addEdge(static_cast<uint8_t>(2), static_cast<uint8_t>(8));
  expectedBN.addEdge(static_cast<uint8_t>(3), static_cast<uint8_t>(9));
  expectedBN.addEdge(static_cast<uint8_t>(3), static_cast<uint8_t>(13));
  expectedBN.addEdge(static_cast<uint8_t>(3), static_cast<uint8_t>(18));
  expectedBN.addEdge(static_cast<uint8_t>(4), static_cast<uint8_t>(6));
  expectedBN.addEdge(static_cast<uint8_t>(4), static_cast<uint8_t>(11));
  expectedBN.addEdge(static_cast<uint8_t>(5), static_cast<uint8_t>(7));
  expectedBN.addEdge(static_cast<uint8_t>(7), static_cast<uint8_t>(5));
  expectedBN.addEdge(static_cast<uint8_t>(8), static_cast<uint8_t>(2));
  expectedBN.addEdge(static_cast<uint8_t>(8), static_cast<uint8_t>(6));
  expectedBN.addEdge(static_cast<uint8_t>(8), static_cast<uint8_t>(11));
  expectedBN.addEdge(static_cast<uint8_t>(12), static_cast<uint8_t>(9));
  expectedBN.addEdge(static_cast<uint8_t>(18), static_cast<uint8_t>(3));
  expectedBN.addEdge(static_cast<uint8_t>(23), static_cast<uint8_t>(24));
  expectedBN.addEdge(static_cast<uint8_t>(24), static_cast<uint8_t>(4));
  // Check false positives
  EXPECT_NE(expectedBN, computedBN);

  expectedBN.addEdge(static_cast<uint8_t>(24), static_cast<uint8_t>(23));
  EXPECT_EQ(expectedBN, computedBN);
}

INSTANTIATE_TEST_CASE_P(Sequential, InsuranceData_GS, testing::Values(false));
INSTANTIATE_TEST_CASE_P(Parallel, InsuranceData_GS, testing::Values(true));

class InsuranceData_MMPC : public InsuranceData<MMPC<DiscreteData<Counter, uint8_t>, uint8_t, UintSet<uint8_t>>>,
                           public testing::WithParamInterface<bool> {
};

TEST_P(InsuranceData_MMPC, DirectedNetwork) {
  auto computedBN = this->algo->getNetwork(true, GetParam());

  auto expectedBN = BayesianNetwork<uint8_t>(this->data->varNames());
  expectedBN.addEdge(static_cast<uint8_t>(0), static_cast<uint8_t>(1));
  expectedBN.addEdge(static_cast<uint8_t>(1), static_cast<uint8_t>(0));
  expectedBN.addEdge(static_cast<uint8_t>(1), static_cast<uint8_t>(13));
  expectedBN.addEdge(static_cast<uint8_t>(2), static_cast<uint8_t>(4));
  expectedBN.addEdge(static_cast<uint8_t>(2), static_cast<uint8_t>(8));
  expectedBN.addEdge(static_cast<uint8_t>(2), static_cast<uint8_t>(17));
  expectedBN.addEdge(static_cast<uint8_t>(2), static_cast<uint8_t>(18));
  expectedBN.addEdge(static_cast<uint8_t>(2), static_cast<uint8_t>(21));
  expectedBN.addEdge(static_cast<uint8_t>(3), static_cast<uint8_t>(9));
  expectedBN.addEdge(static_cast<uint8_t>(3), static_cast<uint8_t>(17));
  expectedBN.addEdge(static_cast<uint8_t>(3), static_cast<uint8_t>(18));
  expectedBN.addEdge(static_cast<uint8_t>(4), static_cast<uint8_t>(11));
  expectedBN.addEdge(static_cast<uint8_t>(4), static_cast<uint8_t>(16));
  expectedBN.addEdge(static_cast<uint8_t>(5), static_cast<uint8_t>(7));
  expectedBN.addEdge(static_cast<uint8_t>(5), static_cast<uint8_t>(14));
  expectedBN.addEdge(static_cast<uint8_t>(6), static_cast<uint8_t>(23));
  expectedBN.addEdge(static_cast<uint8_t>(7), static_cast<uint8_t>(5));
  expectedBN.addEdge(static_cast<uint8_t>(7), static_cast<uint8_t>(20));
  expectedBN.addEdge(static_cast<uint8_t>(8), static_cast<uint8_t>(2));
  expectedBN.addEdge(static_cast<uint8_t>(8), static_cast<uint8_t>(6));
  expectedBN.addEdge(static_cast<uint8_t>(8), static_cast<uint8_t>(16));
  expectedBN.addEdge(static_cast<uint8_t>(10), static_cast<uint8_t>(16));
  expectedBN.addEdge(static_cast<uint8_t>(11), static_cast<uint8_t>(8));
  expectedBN.addEdge(static_cast<uint8_t>(12), static_cast<uint8_t>(9));
  expectedBN.addEdge(static_cast<uint8_t>(12), static_cast<uint8_t>(26));
  expectedBN.addEdge(static_cast<uint8_t>(13), static_cast<uint8_t>(1));
  expectedBN.addEdge(static_cast<uint8_t>(14), static_cast<uint8_t>(5));
  expectedBN.addEdge(static_cast<uint8_t>(14), static_cast<uint8_t>(19));
  expectedBN.addEdge(static_cast<uint8_t>(20), static_cast<uint8_t>(7));
  expectedBN.addEdge(static_cast<uint8_t>(20), static_cast<uint8_t>(19));
  expectedBN.addEdge(static_cast<uint8_t>(21), static_cast<uint8_t>(2));
  expectedBN.addEdge(static_cast<uint8_t>(22), static_cast<uint8_t>(23));
  expectedBN.addEdge(static_cast<uint8_t>(24), static_cast<uint8_t>(4));
  expectedBN.addEdge(static_cast<uint8_t>(24), static_cast<uint8_t>(23));
  // Check false positives
  EXPECT_NE(expectedBN, computedBN);

  expectedBN.addEdge(static_cast<uint8_t>(26), static_cast<uint8_t>(12));
  EXPECT_EQ(expectedBN, computedBN);
}

INSTANTIATE_TEST_CASE_P(Sequential, InsuranceData_MMPC, testing::Values(false));
INSTANTIATE_TEST_CASE_P(Parallel, InsuranceData_MMPC, testing::Values(true));


class MildewData_GS : public MildewData<GS<DiscreteData<Counter, uint8_t>, uint8_t, UintSet<uint8_t>>>,
                      public testing::WithParamInterface<bool> {
};

TEST_P(MildewData_GS, DirectedNetwork) {
  auto computedBN = this->algo->getNetwork(true, GetParam());

  auto expectedBN = BayesianNetwork<uint8_t>(this->data->varNames());
  expectedBN.addEdge(static_cast<uint8_t>(0), static_cast<uint8_t>(1));
  expectedBN.addEdge(static_cast<uint8_t>(1), static_cast<uint8_t>(0));
  expectedBN.addEdge(static_cast<uint8_t>(3), static_cast<uint8_t>(25));
  expectedBN.addEdge(static_cast<uint8_t>(4), static_cast<uint8_t>(6));
  expectedBN.addEdge(static_cast<uint8_t>(6), static_cast<uint8_t>(4));
  expectedBN.addEdge(static_cast<uint8_t>(10), static_cast<uint8_t>(27));
  expectedBN.addEdge(static_cast<uint8_t>(16), static_cast<uint8_t>(29));
  expectedBN.addEdge(static_cast<uint8_t>(27), static_cast<uint8_t>(10));
  expectedBN.addEdge(static_cast<uint8_t>(27), static_cast<uint8_t>(32));
  expectedBN.addEdge(static_cast<uint8_t>(31), static_cast<uint8_t>(25));
  expectedBN.addEdge(static_cast<uint8_t>(32), static_cast<uint8_t>(27));
  // Check false positives
  EXPECT_NE(expectedBN, computedBN);

  expectedBN.addEdge(static_cast<uint8_t>(33), static_cast<uint8_t>(29));
  EXPECT_EQ(expectedBN, computedBN);
}

INSTANTIATE_TEST_CASE_P(Sequential, MildewData_GS, testing::Values(false));
INSTANTIATE_TEST_CASE_P(Parallel, MildewData_GS, testing::Values(true));

class MildewData_MMPC : public MildewData<MMPC<DiscreteData<Counter, uint8_t>, uint8_t, UintSet<uint8_t>>>,
                        public testing::WithParamInterface<bool> {
};

TEST_P(MildewData_MMPC, DirectedNetwork) {
  auto computedBN = this->algo->getNetwork(true, GetParam());

  auto expectedBN = BayesianNetwork<uint8_t>(this->data->varNames());
  expectedBN.addEdge(static_cast<uint8_t>(0), static_cast<uint8_t>(1));
  expectedBN.addEdge(static_cast<uint8_t>(1), static_cast<uint8_t>(0));
  expectedBN.addEdge(static_cast<uint8_t>(3), static_cast<uint8_t>(25));
  expectedBN.addEdge(static_cast<uint8_t>(4), static_cast<uint8_t>(11));
  expectedBN.addEdge(static_cast<uint8_t>(5), static_cast<uint8_t>(12));
  expectedBN.addEdge(static_cast<uint8_t>(8), static_cast<uint8_t>(9));
  expectedBN.addEdge(static_cast<uint8_t>(9), static_cast<uint8_t>(8));
  expectedBN.addEdge(static_cast<uint8_t>(10), static_cast<uint8_t>(27));
  expectedBN.addEdge(static_cast<uint8_t>(11), static_cast<uint8_t>(4));
  expectedBN.addEdge(static_cast<uint8_t>(16), static_cast<uint8_t>(29));
  expectedBN.addEdge(static_cast<uint8_t>(17), static_cast<uint8_t>(23));
  expectedBN.addEdge(static_cast<uint8_t>(18), static_cast<uint8_t>(28));
  expectedBN.addEdge(static_cast<uint8_t>(23), static_cast<uint8_t>(17));
  expectedBN.addEdge(static_cast<uint8_t>(24), static_cast<uint8_t>(30));
  expectedBN.addEdge(static_cast<uint8_t>(26), static_cast<uint8_t>(12));
  expectedBN.addEdge(static_cast<uint8_t>(27), static_cast<uint8_t>(10));
  expectedBN.addEdge(static_cast<uint8_t>(27), static_cast<uint8_t>(32));
  expectedBN.addEdge(static_cast<uint8_t>(28), static_cast<uint8_t>(18));
  expectedBN.addEdge(static_cast<uint8_t>(30), static_cast<uint8_t>(24));
  expectedBN.addEdge(static_cast<uint8_t>(31), static_cast<uint8_t>(25));
  expectedBN.addEdge(static_cast<uint8_t>(32), static_cast<uint8_t>(27));
  // Check false positives
  EXPECT_NE(expectedBN, computedBN);

  expectedBN.addEdge(static_cast<uint8_t>(33), static_cast<uint8_t>(29));
  EXPECT_EQ(expectedBN, computedBN);
}

INSTANTIATE_TEST_CASE_P(Sequential, MildewData_MMPC, testing::Values(false));
INSTANTIATE_TEST_CASE_P(Parallel, MildewData_MMPC, testing::Values(true));


class AlarmData_GS : public AlarmData<GS<DiscreteData<Counter, uint8_t>, uint8_t, UintSet<uint8_t>>>,
                     public testing::WithParamInterface<bool> {
};

TEST_P(AlarmData_GS, DirectedNetwork) {
  auto computedBN = this->algo->getNetwork(true, GetParam());

  auto expectedBN = BayesianNetwork<uint8_t>(this->data->varNames());
  expectedBN.addEdge(static_cast<uint8_t>(0), static_cast<uint8_t>(5));
  expectedBN.addEdge(static_cast<uint8_t>(1), static_cast<uint8_t>(4));
  expectedBN.addEdge(static_cast<uint8_t>(2), static_cast<uint8_t>(4));
  expectedBN.addEdge(static_cast<uint8_t>(3), static_cast<uint8_t>(4));
  expectedBN.addEdge(static_cast<uint8_t>(3), static_cast<uint8_t>(6));
  expectedBN.addEdge(static_cast<uint8_t>(4), static_cast<uint8_t>(1));
  expectedBN.addEdge(static_cast<uint8_t>(4), static_cast<uint8_t>(2));
  expectedBN.addEdge(static_cast<uint8_t>(4), static_cast<uint8_t>(3));
  expectedBN.addEdge(static_cast<uint8_t>(5), static_cast<uint8_t>(0));
  expectedBN.addEdge(static_cast<uint8_t>(5), static_cast<uint8_t>(6));
  expectedBN.addEdge(static_cast<uint8_t>(6), static_cast<uint8_t>(35));
  expectedBN.addEdge(static_cast<uint8_t>(7), static_cast<uint8_t>(8));
  expectedBN.addEdge(static_cast<uint8_t>(9), static_cast<uint8_t>(10));
  expectedBN.addEdge(static_cast<uint8_t>(11), static_cast<uint8_t>(10));
  expectedBN.addEdge(static_cast<uint8_t>(11), static_cast<uint8_t>(34));
  expectedBN.addEdge(static_cast<uint8_t>(14), static_cast<uint8_t>(13));
  expectedBN.addEdge(static_cast<uint8_t>(16), static_cast<uint8_t>(25));
  expectedBN.addEdge(static_cast<uint8_t>(17), static_cast<uint8_t>(15));
  expectedBN.addEdge(static_cast<uint8_t>(19), static_cast<uint8_t>(15));
  expectedBN.addEdge(static_cast<uint8_t>(19), static_cast<uint8_t>(17));
  expectedBN.addEdge(static_cast<uint8_t>(19), static_cast<uint8_t>(20));
  expectedBN.addEdge(static_cast<uint8_t>(20), static_cast<uint8_t>(19));
  expectedBN.addEdge(static_cast<uint8_t>(21), static_cast<uint8_t>(22));
  expectedBN.addEdge(static_cast<uint8_t>(22), static_cast<uint8_t>(21));
  expectedBN.addEdge(static_cast<uint8_t>(23), static_cast<uint8_t>(17));
  expectedBN.addEdge(static_cast<uint8_t>(24), static_cast<uint8_t>(15));
  expectedBN.addEdge(static_cast<uint8_t>(26), static_cast<uint8_t>(25));
  expectedBN.addEdge(static_cast<uint8_t>(27), static_cast<uint8_t>(28));
  expectedBN.addEdge(static_cast<uint8_t>(28), static_cast<uint8_t>(27));
  expectedBN.addEdge(static_cast<uint8_t>(33), static_cast<uint8_t>(14));
  expectedBN.addEdge(static_cast<uint8_t>(34), static_cast<uint8_t>(8));
  expectedBN.addEdge(static_cast<uint8_t>(34), static_cast<uint8_t>(11));
  // Check false positives
  EXPECT_NE(expectedBN, computedBN);

  expectedBN.addEdge(static_cast<uint8_t>(36), static_cast<uint8_t>(14));
  EXPECT_EQ(expectedBN, computedBN);
}

INSTANTIATE_TEST_CASE_P(Sequential, AlarmData_GS, testing::Values(false));
INSTANTIATE_TEST_CASE_P(Parallel, AlarmData_GS, testing::Values(true));

class AlarmData_MMPC : public AlarmData<MMPC<DiscreteData<Counter, uint8_t>, uint8_t, UintSet<uint8_t>>>,
                       public testing::WithParamInterface<bool> {
};

TEST_P(AlarmData_MMPC, DirectedNetwork) {
  auto computedBN = this->algo->getNetwork(true, GetParam());

  auto expectedBN = BayesianNetwork<uint8_t>(this->data->varNames());
  expectedBN.addEdge(static_cast<uint8_t>(0), static_cast<uint8_t>(5));
  expectedBN.addEdge(static_cast<uint8_t>(1), static_cast<uint8_t>(4));
  expectedBN.addEdge(static_cast<uint8_t>(2), static_cast<uint8_t>(4));
  expectedBN.addEdge(static_cast<uint8_t>(3), static_cast<uint8_t>(4));
  expectedBN.addEdge(static_cast<uint8_t>(3), static_cast<uint8_t>(6));
  expectedBN.addEdge(static_cast<uint8_t>(4), static_cast<uint8_t>(1));
  expectedBN.addEdge(static_cast<uint8_t>(4), static_cast<uint8_t>(2));
  expectedBN.addEdge(static_cast<uint8_t>(4), static_cast<uint8_t>(3));
  expectedBN.addEdge(static_cast<uint8_t>(5), static_cast<uint8_t>(0));
  expectedBN.addEdge(static_cast<uint8_t>(5), static_cast<uint8_t>(6));
  expectedBN.addEdge(static_cast<uint8_t>(6), static_cast<uint8_t>(35));
  expectedBN.addEdge(static_cast<uint8_t>(7), static_cast<uint8_t>(8));
  expectedBN.addEdge(static_cast<uint8_t>(10), static_cast<uint8_t>(9));
  expectedBN.addEdge(static_cast<uint8_t>(10), static_cast<uint8_t>(11));
  expectedBN.addEdge(static_cast<uint8_t>(13), static_cast<uint8_t>(14));
  expectedBN.addEdge(static_cast<uint8_t>(14), static_cast<uint8_t>(13));
  expectedBN.addEdge(static_cast<uint8_t>(14), static_cast<uint8_t>(33));
  expectedBN.addEdge(static_cast<uint8_t>(14), static_cast<uint8_t>(36));
  expectedBN.addEdge(static_cast<uint8_t>(15), static_cast<uint8_t>(30));
  expectedBN.addEdge(static_cast<uint8_t>(15), static_cast<uint8_t>(32));
  expectedBN.addEdge(static_cast<uint8_t>(16), static_cast<uint8_t>(25));
  expectedBN.addEdge(static_cast<uint8_t>(17), static_cast<uint8_t>(24));
  expectedBN.addEdge(static_cast<uint8_t>(17), static_cast<uint8_t>(30));
  expectedBN.addEdge(static_cast<uint8_t>(18), static_cast<uint8_t>(19));
  expectedBN.addEdge(static_cast<uint8_t>(19), static_cast<uint8_t>(18));
  expectedBN.addEdge(static_cast<uint8_t>(19), static_cast<uint8_t>(20));
  expectedBN.addEdge(static_cast<uint8_t>(21), static_cast<uint8_t>(22));
  expectedBN.addEdge(static_cast<uint8_t>(22), static_cast<uint8_t>(21));
  expectedBN.addEdge(static_cast<uint8_t>(22), static_cast<uint8_t>(23));
  expectedBN.addEdge(static_cast<uint8_t>(23), static_cast<uint8_t>(20));
  expectedBN.addEdge(static_cast<uint8_t>(24), static_cast<uint8_t>(23));
  expectedBN.addEdge(static_cast<uint8_t>(26), static_cast<uint8_t>(29));
  expectedBN.addEdge(static_cast<uint8_t>(27), static_cast<uint8_t>(28));
  expectedBN.addEdge(static_cast<uint8_t>(28), static_cast<uint8_t>(27));
  expectedBN.addEdge(static_cast<uint8_t>(28), static_cast<uint8_t>(29));
  expectedBN.addEdge(static_cast<uint8_t>(29), static_cast<uint8_t>(25));
  expectedBN.addEdge(static_cast<uint8_t>(30), static_cast<uint8_t>(15));
  expectedBN.addEdge(static_cast<uint8_t>(30), static_cast<uint8_t>(17));
  expectedBN.addEdge(static_cast<uint8_t>(31), static_cast<uint8_t>(24));
  expectedBN.addEdge(static_cast<uint8_t>(31), static_cast<uint8_t>(32));
  expectedBN.addEdge(static_cast<uint8_t>(33), static_cast<uint8_t>(14));
  expectedBN.addEdge(static_cast<uint8_t>(34), static_cast<uint8_t>(8));
  expectedBN.addEdge(static_cast<uint8_t>(34), static_cast<uint8_t>(9));
  expectedBN.addEdge(static_cast<uint8_t>(34), static_cast<uint8_t>(11));
  expectedBN.addEdge(static_cast<uint8_t>(34), static_cast<uint8_t>(35));
  // Check false positives
  EXPECT_NE(expectedBN, computedBN);

  expectedBN.addEdge(static_cast<uint8_t>(35), static_cast<uint8_t>(36));
  EXPECT_EQ(expectedBN, computedBN);
}

INSTANTIATE_TEST_CASE_P(Sequential, AlarmData_MMPC, testing::Values(false));
INSTANTIATE_TEST_CASE_P(Parallel, AlarmData_MMPC, testing::Values(true));

#endif // TEST_SEQUENTIALNETWORK_HPP_
