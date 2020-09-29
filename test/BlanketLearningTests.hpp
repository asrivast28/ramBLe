/**
 * @file BlanketLearningTests.hpp
 * @brief Unit tests for the blanket learning algorithms.
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
#ifndef TEST_BLANKETLEARNING_HPP_
#define TEST_BLANKETLEARNING_HPP_

#include "Environment.hpp"
#include "BlanketLearning.hpp"
#include "CTCounter.hpp"
#include "DiscreteData.hpp"
#include "UintSet.hpp"


using Counter = CTCounter<>;
// All the different blanket learning algorithms
using BlanketLearningAlgorithms = testing::Types<GS<DiscreteData<Counter, uint8_t>, uint8_t, UintSet<uint8_t>>,
                                                 IAMB<DiscreteData<Counter, uint8_t>, uint8_t, UintSet<uint8_t>>,
                                                 InterIAMB<DiscreteData<Counter, uint8_t>, uint8_t, UintSet<uint8_t>>>;

template <typename Algorithm>
class CoronaryBlanketLearning : public testing::Test { 
protected:
  void
  SetUp() override {
    auto n = CoronaryEnvironment::n;
    auto m = CoronaryEnvironment::m;
    auto reader = CoronaryEnvironment::reader;
    auto counter = Counter::create(n, m, std::begin(reader->data()));
    data = new DiscreteData<Counter, uint8_t>(counter, reader->varNames());
    algo = new Algorithm(comm, *data);
  }

  void
  TearDown() override {
    delete algo;
    delete data;
  }

  mxx::comm comm;
  DiscreteData<Counter, uint8_t>* data;
  Algorithm* algo;
};

TYPED_TEST_SUITE(CoronaryBlanketLearning, BlanketLearningAlgorithms);

TYPED_TEST(CoronaryBlanketLearning, MarkovBlanket) {
  auto target = this->data->varIndex("Smoking");
  auto trueSmokingMB = this->data->template varIndices<UintSet<uint8_t>>({"M. Work", "P. Work", "Pressure", "Proteins"});
  auto computedSmokingMB = this->algo->getMB(target);
  EXPECT_EQ(computedSmokingMB, trueSmokingMB);

  target = this->data->varIndex("M. Work");
  auto trueMWorkMB = this->data->template varIndices<UintSet<uint8_t>>({"Smoking", "P. Work", "Pressure", "Proteins", "Family"});
  auto computedMWorkMB = this->algo->getMB(target);
  EXPECT_EQ(computedMWorkMB, trueMWorkMB);

  target = this->data->varIndex("P. Work");
  auto truePWorkMB = this->data->template varIndices<UintSet<uint8_t>>({"Smoking", "M. Work", "Pressure", "Proteins"});
  auto computedPWorkMB = this->algo->getMB(target);
  EXPECT_EQ(computedPWorkMB, truePWorkMB);

  target = this->data->varIndex("Pressure");
  auto truePressureMB = this->data->template varIndices<UintSet<uint8_t>>({"Smoking", "M. Work", "P. Work", "Proteins"});
  auto computedPressureMB = this->algo->getMB(target);
  EXPECT_EQ(computedPressureMB, truePressureMB);

  target = this->data->varIndex("Proteins");
  auto trueProteinsMB = this->data->template varIndices<UintSet<uint8_t>>({"Smoking", "M. Work", "P. Work", "Pressure"});
  auto computedProteinsMB = this->algo->getMB(target);
  EXPECT_EQ(computedProteinsMB, trueProteinsMB);

  target = this->data->varIndex("Family");
  auto trueFamilyMB = this->data->template varIndices<UintSet<uint8_t>>({"M. Work"});
  auto computedFamilyMB = this->algo->getMB(target);
  EXPECT_EQ(computedFamilyMB, trueFamilyMB);
}

TYPED_TEST(CoronaryBlanketLearning, ParentsChildren) {
  auto target = this->data->varIndex("Smoking");
  auto trueSmokingPC = this->data->template varIndices<UintSet<uint8_t>>({"M. Work", "P. Work", "Pressure", "Proteins"});
  auto computedSmokingPC = this->algo->getPC(target);
  EXPECT_EQ(computedSmokingPC, trueSmokingPC);

  target = this->data->varIndex("M. Work");
  auto trueMWorkPC = this->data->template varIndices<UintSet<uint8_t>>({"Smoking", "P. Work", "Pressure", "Proteins", "Family"});
  auto computedMWorkPC = this->algo->getPC(target);
  EXPECT_EQ(computedMWorkPC, trueMWorkPC);

  target = this->data->varIndex("P. Work");
  auto truePWorkPC = this->data->template varIndices<UintSet<uint8_t>>({"Smoking", "M. Work"});
  auto computedPWorkPC = this->algo->getPC(target);
  EXPECT_EQ(computedPWorkPC, truePWorkPC);

  target = this->data->varIndex("Pressure");
  auto truePressurePC = this->data->template varIndices<UintSet<uint8_t>>({"Smoking", "M. Work", "Proteins"});
  auto computedPressurePC = this->algo->getPC(target);
  EXPECT_EQ(computedPressurePC, truePressurePC);

  target = this->data->varIndex("Proteins");
  auto trueProteinsPC = this->data->template varIndices<UintSet<uint8_t>>({"Smoking", "M. Work", "Pressure"});
  auto computedProteinsPC = this->algo->getPC(target);
  EXPECT_EQ(computedProteinsPC, trueProteinsPC);

  target = this->data->varIndex("Family");
  auto trueFamilyPC = this->data->template varIndices<UintSet<uint8_t>>({"M. Work"});
  auto computedFamilyPC = this->algo->getPC(target);
  EXPECT_EQ(computedFamilyPC, trueFamilyPC);
}

TYPED_TEST(CoronaryBlanketLearning, DirectedNetwork) {
  auto computedBN = this->algo->getNetwork(true, false);

  auto smoking = this->data->varIndex("Smoking");
  auto mWork = this->data->varIndex("M. Work");
  auto pWork  = this->data->varIndex("P. Work");
  auto pressure = this->data->varIndex("Pressure");
  auto proteins = this->data->varIndex("Proteins");
  auto family = this->data->varIndex("Family");

  auto expectedBN = BayesianNetwork<uint8_t>(this->data->varNames());
  expectedBN.addEdge(smoking, mWork);
  expectedBN.addEdge(smoking, pressure);
  expectedBN.addEdge(mWork, proteins);
  expectedBN.addEdge(mWork, family);
  expectedBN.addEdge(pWork, smoking);
  expectedBN.addEdge(pWork, mWork);
  expectedBN.addEdge(pressure, mWork);
  expectedBN.addEdge(proteins, smoking);
  expectedBN.addEdge(proteins, mWork);
  // Check for false positives
  EXPECT_NE(expectedBN, computedBN);

  expectedBN.addEdge(proteins, pressure);
  EXPECT_EQ(expectedBN, computedBN);
}


template <typename Algorithm>
class AsiaBlanketLearning : public testing::Test {
protected:
  void
  SetUp() override {
    auto n = AsiaEnvironment::n;
    auto m = AsiaEnvironment::m;
    auto reader = AsiaEnvironment::reader;
    auto counter = Counter::create(n, m, std::begin(reader->data()));
    data = new DiscreteData<Counter, uint8_t>(counter, reader->varNames());
    algo = new Algorithm(comm, *data);
  }

  void
  TearDown() override {
    delete algo;
    delete data;
  }

  mxx::comm comm;
  DiscreteData<Counter, uint8_t>* data;
  Algorithm* algo;
};

TYPED_TEST_SUITE(AsiaBlanketLearning, BlanketLearningAlgorithms);

TYPED_TEST(AsiaBlanketLearning, MarkovBlanket) {
  auto target = this->data->varIndex("asia");
  auto trueAsiaMB = this->data->template varIndices<UintSet<uint8_t>>({});
  auto computedAsiaMB = this->algo->getMB(target);
  EXPECT_EQ(computedAsiaMB, trueAsiaMB);

  target = this->data->varIndex("smoke");
  auto trueSmokeMB = this->data->template varIndices<UintSet<uint8_t>>({"bronc"});
  auto computedSmokeMB = this->algo->getMB(target);
  EXPECT_EQ(computedSmokeMB, trueSmokeMB);

  target = this->data->varIndex("tub");
  auto trueTubMB = this->data->template varIndices<UintSet<uint8_t>>({"either", "lung"});
  auto computedTubMB = this->algo->getMB(target);
  EXPECT_EQ(computedTubMB, trueTubMB);

  target = this->data->varIndex("lung");
  auto trueLungMB = this->data->template varIndices<UintSet<uint8_t>>({"either", "tub"});
  auto computedLungMB = this->algo->getMB(target);
  EXPECT_EQ(computedLungMB, trueLungMB);

  target = this->data->varIndex("bronc");
  auto trueBroncMB = this->data->template varIndices<UintSet<uint8_t>>({"dysp", "smoke"});
  auto computedBroncMB = this->algo->getMB(target);
  EXPECT_EQ(computedBroncMB, trueBroncMB);

  target = this->data->varIndex("either");
  auto trueEitherMB = this->data->template varIndices<UintSet<uint8_t>>({"lung", "tub"});
  auto computedEitherMB = this->algo->getMB(target);
  EXPECT_EQ(computedEitherMB, trueEitherMB);

  target = this->data->varIndex("xray");
  auto trueXrayMB = this->data->template varIndices<UintSet<uint8_t>>({});
  auto computedXrayMB = this->algo->getMB(target);
  EXPECT_EQ(computedXrayMB, trueXrayMB);

  target = this->data->varIndex("dysp");
  auto trueDyspMB = this->data->template varIndices<UintSet<uint8_t>>({"bronc"});
  auto computedDyspMB = this->algo->getMB(target);
  EXPECT_EQ(computedDyspMB, trueDyspMB);
}

TYPED_TEST(AsiaBlanketLearning, ParentsChildren) {
  auto target = this->data->varIndex("asia");
  auto trueAsiaPC = this->data->template varIndices<UintSet<uint8_t>>({});
  auto computedAsiaPC = this->algo->getPC(target);
  EXPECT_EQ(computedAsiaPC, trueAsiaPC);

  target = this->data->varIndex("smoke");
  auto trueSmokePC = this->data->template varIndices<UintSet<uint8_t>>({"bronc"});
  auto computedSmokePC = this->algo->getPC(target);
  EXPECT_EQ(computedSmokePC, trueSmokePC);

  target = this->data->varIndex("tub");
  auto trueTubPC = this->data->template varIndices<UintSet<uint8_t>>({"either"});
  auto computedTubPC = this->algo->getPC(target);
  EXPECT_EQ(computedTubPC, trueTubPC);

  target = this->data->varIndex("lung");
  auto trueLungPC = this->data->template varIndices<UintSet<uint8_t>>({"either"});
  auto computedLungPC = this->algo->getPC(target);
  EXPECT_EQ(computedLungPC, trueLungPC);

  target = this->data->varIndex("bronc");
  auto trueBroncPC = this->data->template varIndices<UintSet<uint8_t>>({"dysp", "smoke"});
  auto computedBroncPC = this->algo->getPC(target);
  EXPECT_EQ(computedBroncPC, trueBroncPC);

  target = this->data->varIndex("either");
  auto trueEitherPC = this->data->template varIndices<UintSet<uint8_t>>({"lung", "tub"});
  auto computedEitherPC = this->algo->getPC(target);
  EXPECT_EQ(computedEitherPC, trueEitherPC);

  target = this->data->varIndex("xray");
  auto trueXrayPC = this->data->template varIndices<UintSet<uint8_t>>({});
  auto computedXrayPC = this->algo->getPC(target);
  EXPECT_EQ(computedXrayPC, trueXrayPC);

  target = this->data->varIndex("dysp");
  auto trueDyspPC = this->data->template varIndices<UintSet<uint8_t>>({"bronc"});
  auto computedDyspPC = this->algo->getPC(target);
  EXPECT_EQ(computedDyspPC, trueDyspPC);
}

TYPED_TEST(AsiaBlanketLearning, DirectedNetwork) {
  auto computedBN = this->algo->getNetwork(true, false);

  auto smoke = this->data->varIndex("smoke");
  auto tub = this->data->varIndex("tub");
  auto lung = this->data->varIndex("lung");
  auto bronc = this->data->varIndex("bronc");
  auto either = this->data->varIndex("either");
  auto dysp = this->data->varIndex("dysp");

  auto expectedBN = BayesianNetwork<uint8_t>(this->data->varNames());
  expectedBN.addEdge(smoke, bronc);
  expectedBN.addEdge(tub, either);
  expectedBN.addEdge(lung, either);
  // Check for false positives
  EXPECT_NE(expectedBN, computedBN);

  expectedBN.addEdge(dysp, bronc);
  EXPECT_EQ(expectedBN, computedBN);
}

#endif // TEST_BLANKETLEARNING_HPP_
