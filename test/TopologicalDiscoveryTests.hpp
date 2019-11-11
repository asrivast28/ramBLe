/**
 * @file TopologicalDiscoveryTests.hpp
 * @brief Unit tests for the Topological MB discovery algorithms.
 */
#ifndef TEST_TOPOLOGICALDISCOVERY_HPP_
#define TEST_TOPOLOGICALDISCOVERY_HPP_

#include "Environment.hpp"
#include "DiscreteData.hpp"
#include "TopologicalDiscovery.hpp"


using Counter = CTCounter<1>;
// All the different direct discovery algorithms
using TopologicalDiscoveryAlgorithms = testing::Types<MMPC<DiscreteData<Counter, uint8_t>, uint8_t, UintSet<uint8_t>>,
                                                      SemiInterleavedHITON<DiscreteData<Counter, uint8_t>, uint8_t, UintSet<uint8_t>>,
                                                      GetPC<DiscreteData<Counter, uint8_t>, uint8_t, UintSet<uint8_t>>>;

template <typename Algorithm>
class CoronaryTopologicalDiscovery : public testing::Test {
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

TYPED_TEST_CASE(CoronaryTopologicalDiscovery, TopologicalDiscoveryAlgorithms);

TYPED_TEST(CoronaryTopologicalDiscovery, ParentsChildren) {
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


template <typename Algorithm>
class AsiaTopologicalDiscovery : public testing::Test {
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

TYPED_TEST_CASE(AsiaTopologicalDiscovery, TopologicalDiscoveryAlgorithms);

TYPED_TEST(AsiaTopologicalDiscovery, ParentsChildren) {
  auto target = this->data->varIndex("asia");
  auto trueAsiaPC = this->data->template varIndices<UintSet<uint8_t>>({});
  auto computedAsiaPC = this->algo->getPC(target);
  EXPECT_EQ(computedAsiaPC, trueAsiaPC);

  target = this->data->varIndex("smoke");
  auto trueSmokePC = this->data->template varIndices<UintSet<uint8_t>>({"bronc", "lung"});
  auto computedSmokePC = this->algo->getPC(target);
  EXPECT_EQ(computedSmokePC, trueSmokePC);

  target = this->data->varIndex("tub");
  auto trueTubPC = this->data->template varIndices<UintSet<uint8_t>>({"either"});
  auto computedTubPC = this->algo->getPC(target);
  EXPECT_EQ(computedTubPC, trueTubPC);

  target = this->data->varIndex("lung");
  auto trueLungPC = this->data->template varIndices<UintSet<uint8_t>>({"either", "smoke"});
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

#endif // TEST_TOPOLOGICALDISCOVERY_HPP_
