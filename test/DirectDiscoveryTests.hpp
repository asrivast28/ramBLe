/**
 * @file DirectDiscoveryTests.hpp
 * @brief Unit tests for the Direct MB discovery algorithms.
 */
#ifndef TEST_DIRECTDISCOVERY_HPP_
#define TEST_DIRECTDISCOVERY_HPP_

#include "Environment.hpp"
#include "CTCounter.hpp"
#include "DataQuery.hpp"
#include "DirectDiscovery.hpp"
#include "UintSet.hpp"


using Counter = CTCounter<1>;
// All the different direct discovery algorithms
using DirectDiscoveryAlgorithms = testing::Types<GSMB<DataQuery<Counter, uint8_t>, uint8_t, UintSet<uint8_t>>,
                                                 IAMB<DataQuery<Counter, uint8_t>, uint8_t, UintSet<uint8_t>>,
                                                 InterIAMB<DataQuery<Counter, uint8_t>, uint8_t, UintSet<uint8_t>>>;

template <typename Algorithm>
class CoronaryDirectDiscovery : public testing::Test {
protected:
  void
  SetUp() override {
    auto n = CoronaryEnvironment::n;
    auto m = CoronaryEnvironment::m;
    auto reader = CoronaryEnvironment::reader;
    auto counter = Counter::create(n, m, std::begin(reader->data()));
    data = new DataQuery<Counter, uint8_t>(counter, reader->varNames());
    algo = new Algorithm(*data);
  }

  void
  TearDown() override {
    delete algo;
    delete data;
  }

  DataQuery<Counter, uint8_t>* data;
  Algorithm* algo;
};

TYPED_TEST_CASE(CoronaryDirectDiscovery, DirectDiscoveryAlgorithms);

TYPED_TEST(CoronaryDirectDiscovery, MarkovBlanket) {
  uint8_t target = this->data->varIndex("Smoking");
  UintSet<uint8_t> trueSmokingMB = this->data->template varIndices<UintSet<uint8_t>>({"M. Work", "P. Work", "Pressure", "Proteins"});
  UintSet<uint8_t> computedSmokingMB = this->algo->getMB(target);
  EXPECT_EQ(computedSmokingMB, trueSmokingMB);

  target = this->data->varIndex("M. Work");
  UintSet<uint8_t> trueMWorkMB = this->data->template varIndices<UintSet<uint8_t>>({"Smoking", "P. Work", "Pressure", "Proteins", "Family"});
  UintSet<uint8_t> computedMWorkMB = this->algo->getMB(target);
  EXPECT_EQ(computedMWorkMB, trueMWorkMB);

  target = this->data->varIndex("P. Work");
  UintSet<uint8_t> truePWorkMB = this->data->template varIndices<UintSet<uint8_t>>({"Smoking", "M. Work", "Pressure", "Proteins"});
  UintSet<uint8_t> computedPWorkMB = this->algo->getMB(target);
  EXPECT_EQ(computedPWorkMB, truePWorkMB);

  target = this->data->varIndex("Pressure");
  UintSet<uint8_t> truePressureMB = this->data->template varIndices<UintSet<uint8_t>>({"Smoking", "M. Work", "P. Work", "Proteins"});
  UintSet<uint8_t> computedPressureMB = this->algo->getMB(target);
  EXPECT_EQ(computedPressureMB, truePressureMB);

  target = this->data->varIndex("Proteins");
  UintSet<uint8_t> trueProteinsMB = this->data->template varIndices<UintSet<uint8_t>>({"Smoking", "M. Work", "P. Work", "Pressure"});
  UintSet<uint8_t> computedProteinsMB = this->algo->getMB(target);
  EXPECT_EQ(computedProteinsMB, trueProteinsMB);

  target = this->data->varIndex("Family");
  UintSet<uint8_t> trueFamilyMB = this->data->template varIndices<UintSet<uint8_t>>({"M. Work"});
  UintSet<uint8_t> computedFamilyMB = this->algo->getMB(target);
  EXPECT_EQ(computedFamilyMB, trueFamilyMB);
}

TYPED_TEST(CoronaryDirectDiscovery, ParentsChildren) {
  uint8_t target = this->data->varIndex("Smoking");
  UintSet<uint8_t> trueSmokingPC = this->data->template varIndices<UintSet<uint8_t>>({"M. Work", "P. Work", "Pressure", "Proteins"});
  UintSet<uint8_t> computedSmokingPC = this->algo->getPC(target);
  EXPECT_EQ(computedSmokingPC, trueSmokingPC);

  target = this->data->varIndex("M. Work");
  UintSet<uint8_t> trueMWorkPC = this->data->template varIndices<UintSet<uint8_t>>({"Smoking", "P. Work", "Pressure", "Proteins", "Family"});
  UintSet<uint8_t> computedMWorkPC = this->algo->getPC(target);
  EXPECT_EQ(computedMWorkPC, trueMWorkPC);

  target = this->data->varIndex("P. Work");
  UintSet<uint8_t> truePWorkPC = this->data->template varIndices<UintSet<uint8_t>>({"Smoking", "M. Work"});
  UintSet<uint8_t> computedPWorkPC = this->algo->getPC(target);
  EXPECT_EQ(computedPWorkPC, truePWorkPC);

  target = this->data->varIndex("Pressure");
  UintSet<uint8_t> truePressurePC = this->data->template varIndices<UintSet<uint8_t>>({"Smoking", "M. Work", "Proteins"});
  UintSet<uint8_t> computedPressurePC = this->algo->getPC(target);
  EXPECT_EQ(computedPressurePC, truePressurePC);

  target = this->data->varIndex("Proteins");
  UintSet<uint8_t> trueProteinsPC = this->data->template varIndices<UintSet<uint8_t>>({"Smoking", "M. Work", "Pressure"});
  UintSet<uint8_t> computedProteinsPC = this->algo->getPC(target);
  EXPECT_EQ(computedProteinsPC, trueProteinsPC);

  target = this->data->varIndex("Family");
  UintSet<uint8_t> trueFamilyPC = this->data->template varIndices<UintSet<uint8_t>>({"M. Work"});
  UintSet<uint8_t> computedFamilyPC = this->algo->getPC(target);
  EXPECT_EQ(computedFamilyPC, trueFamilyPC);
}


template <typename Algorithm>
class AsiaDirectDiscovery : public testing::Test {
protected:
  void
  SetUp() override {
    auto n = AsiaEnvironment::n;
    auto m = AsiaEnvironment::m;
    auto reader = AsiaEnvironment::reader;
    auto counter = Counter::create(n, m, std::begin(reader->data()));
    data = new DataQuery<Counter, uint8_t>(counter, reader->varNames());
    algo = new Algorithm(*data);
  }

  void
  TearDown() override {
    delete algo;
    delete data;
  }

  DataQuery<Counter, uint8_t>* data;
  Algorithm* algo;
};

TYPED_TEST_CASE(AsiaDirectDiscovery, DirectDiscoveryAlgorithms);

TYPED_TEST(AsiaDirectDiscovery, MarkovBlanket) {
  uint8_t target = this->data->varIndex("asia");
  UintSet<uint8_t> trueAsiaMB = this->data->template varIndices<UintSet<uint8_t>>({});
  UintSet<uint8_t> computedAsiaMB = this->algo->getMB(target);
  EXPECT_EQ(computedAsiaMB, trueAsiaMB);

  target = this->data->varIndex("smoke");
  UintSet<uint8_t> trueSmokeMB = this->data->template varIndices<UintSet<uint8_t>>({"bronc"});
  UintSet<uint8_t> computedSmokeMB = this->algo->getMB(target);
  EXPECT_EQ(computedSmokeMB, trueSmokeMB);

  target = this->data->varIndex("tub");
  UintSet<uint8_t> trueTubMB = this->data->template varIndices<UintSet<uint8_t>>({"either", "lung"});
  UintSet<uint8_t> computedTubMB = this->algo->getMB(target);
  EXPECT_EQ(computedTubMB, trueTubMB);

  target = this->data->varIndex("lung");
  UintSet<uint8_t> trueLungMB = this->data->template varIndices<UintSet<uint8_t>>({"either", "tub"});
  UintSet<uint8_t> computedLungMB = this->algo->getMB(target);
  EXPECT_EQ(computedLungMB, trueLungMB);

  target = this->data->varIndex("bronc");
  UintSet<uint8_t> trueBroncMB = this->data->template varIndices<UintSet<uint8_t>>({"dysp", "smoke"});
  UintSet<uint8_t> computedBroncMB = this->algo->getMB(target);
  EXPECT_EQ(computedBroncMB, trueBroncMB);

  target = this->data->varIndex("either");
  UintSet<uint8_t> trueEitherMB = this->data->template varIndices<UintSet<uint8_t>>({"lung", "tub"});
  UintSet<uint8_t> computedEitherMB = this->algo->getMB(target);
  EXPECT_EQ(computedEitherMB, trueEitherMB);

  target = this->data->varIndex("xray");
  UintSet<uint8_t> trueXrayMB = this->data->template varIndices<UintSet<uint8_t>>({});
  UintSet<uint8_t> computedXrayMB = this->algo->getMB(target);
  EXPECT_EQ(computedXrayMB, trueXrayMB);

  target = this->data->varIndex("dysp");
  UintSet<uint8_t> trueDyspMB = this->data->template varIndices<UintSet<uint8_t>>({"bronc"});
  UintSet<uint8_t> computedDyspMB = this->algo->getMB(target);
  EXPECT_EQ(computedDyspMB, trueDyspMB);
}

TYPED_TEST(AsiaDirectDiscovery, ParentsChildren) {
  uint8_t target = this->data->varIndex("asia");
  UintSet<uint8_t> trueAsiaPC = this->data->template varIndices<UintSet<uint8_t>>({});
  UintSet<uint8_t> computedAsiaPC = this->algo->getPC(target);
  EXPECT_EQ(computedAsiaPC, trueAsiaPC);

  target = this->data->varIndex("smoke");
  UintSet<uint8_t> trueSmokePC = this->data->template varIndices<UintSet<uint8_t>>({"bronc"});
  UintSet<uint8_t> computedSmokePC = this->algo->getPC(target);
  EXPECT_EQ(computedSmokePC, trueSmokePC);

  target = this->data->varIndex("tub");
  UintSet<uint8_t> trueTubPC = this->data->template varIndices<UintSet<uint8_t>>({"either"});
  UintSet<uint8_t> computedTubPC = this->algo->getPC(target);
  EXPECT_EQ(computedTubPC, trueTubPC);

  target = this->data->varIndex("lung");
  UintSet<uint8_t> trueLungPC = this->data->template varIndices<UintSet<uint8_t>>({"either"});
  UintSet<uint8_t> computedLungPC = this->algo->getPC(target);
  EXPECT_EQ(computedLungPC, trueLungPC);

  target = this->data->varIndex("bronc");
  UintSet<uint8_t> trueBroncPC = this->data->template varIndices<UintSet<uint8_t>>({"dysp", "smoke"});
  UintSet<uint8_t> computedBroncPC = this->algo->getPC(target);
  EXPECT_EQ(computedBroncPC, trueBroncPC);

  target = this->data->varIndex("either");
  UintSet<uint8_t> trueEitherPC = this->data->template varIndices<UintSet<uint8_t>>({"lung", "tub"});
  UintSet<uint8_t> computedEitherPC = this->algo->getPC(target);
  EXPECT_EQ(computedEitherPC, trueEitherPC);

  target = this->data->varIndex("xray");
  UintSet<uint8_t> trueXrayPC = this->data->template varIndices<UintSet<uint8_t>>({});
  UintSet<uint8_t> computedXrayPC = this->algo->getPC(target);
  EXPECT_EQ(computedXrayPC, trueXrayPC);

  target = this->data->varIndex("dysp");
  UintSet<uint8_t> trueDyspPC = this->data->template varIndices<UintSet<uint8_t>>({"bronc"});
  UintSet<uint8_t> computedDyspPC = this->algo->getPC(target);
  EXPECT_EQ(computedDyspPC, trueDyspPC);
}

#endif // TEST_DIRECTDISCOVERY_HPP_
