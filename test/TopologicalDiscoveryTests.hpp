/**
 * @file TopologicalDiscoveryTests.hpp
 * @brief Unit tests for the Topological MB discovery algorithms.
 */
#ifndef TEST_TOPOLOGICALDISCOVERY_HPP_
#define TEST_TOPOLOGICALDISCOVERY_HPP_

#include "Common.hpp"
#include "TopologicalDiscovery.hpp"


TEST_F(CoronaryTest, TopologicalDiscovery_MMPC) {
  MMPC<decltype(data), uint8_t> mm(data);

  uint8_t target = data.varIndex("Smoking");
  std::set<uint8_t> trueSmokingPC = data.varIndices({"M. Work", "P. Work", "Pressure", "Proteins"});
  std::set<uint8_t> computedSmokingPC = mm.getCorrectPC(target);
  EXPECT_EQ(computedSmokingPC, trueSmokingPC);

  target = data.varIndex("M. Work");
  std::set<uint8_t> trueMWorkPC = data.varIndices({"Smoking", "P. Work", "Pressure", "Proteins", "Family"});
  std::set<uint8_t> computedMWorkPC = mm.getCorrectPC(target);
  EXPECT_EQ(computedMWorkPC, trueMWorkPC);

  target = data.varIndex("P. Work");
  std::set<uint8_t> truePWorkPC = data.varIndices({"Smoking", "M. Work"});
  std::set<uint8_t> computedPWorkPC = mm.getCorrectPC(target);
  EXPECT_EQ(computedPWorkPC, truePWorkPC);

  target = data.varIndex("Pressure");
  std::set<uint8_t> truePressurePC = data.varIndices({"Smoking", "M. Work", "Proteins"});
  std::set<uint8_t> computedPressurePC = mm.getCorrectPC(target);
  EXPECT_EQ(computedPressurePC, truePressurePC);

  target = data.varIndex("Proteins");
  std::set<uint8_t> trueProteinsPC = data.varIndices({"Smoking", "M. Work", "Pressure"});
  std::set<uint8_t> computedProteinsPC = mm.getCorrectPC(target);
  EXPECT_EQ(computedProteinsPC, trueProteinsPC);

  target = data.varIndex("Family");
  std::set<uint8_t> trueFamilyPC = data.varIndices({"M. Work"});
  std::set<uint8_t> computedFamilyPC = mm.getCorrectPC(target);
  EXPECT_EQ(computedFamilyPC, trueFamilyPC);
}

TEST_F(CoronaryTest, TopologicalDiscovery_HITON) {
  HITON<decltype(data), uint8_t> hiton(data);

  uint8_t target = data.varIndex("Smoking");
  std::set<uint8_t> trueSmokingPC = data.varIndices({"M. Work", "P. Work", "Pressure", "Proteins"});
  std::set<uint8_t> computedSmokingPC = hiton.getCorrectPC(target);
  EXPECT_EQ(computedSmokingPC, trueSmokingPC);

  target = data.varIndex("M. Work");
  std::set<uint8_t> trueMWorkPC = data.varIndices({"Smoking", "P. Work", "Pressure", "Proteins", "Family"});
  std::set<uint8_t> computedMWorkPC = hiton.getCorrectPC(target);
  EXPECT_EQ(computedMWorkPC, trueMWorkPC);

  target = data.varIndex("P. Work");
  std::set<uint8_t> truePWorkPC = data.varIndices({"Smoking", "M. Work"});
  std::set<uint8_t> computedPWorkPC = hiton.getCorrectPC(target);
  EXPECT_EQ(computedPWorkPC, truePWorkPC);

  target = data.varIndex("Pressure");
  std::set<uint8_t> truePressurePC = data.varIndices({"Smoking", "M. Work", "Proteins"});
  std::set<uint8_t> computedPressurePC = hiton.getCorrectPC(target);
  EXPECT_EQ(computedPressurePC, truePressurePC);

  target = data.varIndex("Proteins");
  std::set<uint8_t> trueProteinsPC = data.varIndices({"Smoking", "M. Work", "Pressure"});
  std::set<uint8_t> computedProteinsPC = hiton.getCorrectPC(target);
  EXPECT_EQ(computedProteinsPC, trueProteinsPC);

  target = data.varIndex("Family");
  std::set<uint8_t> trueFamilyPC = data.varIndices({"M. Work"});
  std::set<uint8_t> computedFamilyPC = hiton.getCorrectPC(target);
  EXPECT_EQ(computedFamilyPC, trueFamilyPC);
}

TEST_F(CoronaryTest, TopologicalDiscovery_GetPC) {
  GetPC<decltype(data), uint8_t> gmb(data);

  uint8_t target = data.varIndex("Smoking");
  std::set<uint8_t> trueSmokingPC = data.varIndices({"M. Work", "P. Work", "Pressure", "Proteins"});
  std::set<uint8_t> computedSmokingPC = gmb.getCorrectPC(target);
  EXPECT_EQ(computedSmokingPC, trueSmokingPC);

  target = data.varIndex("M. Work");
  std::set<uint8_t> trueMWorkPC = data.varIndices({"Smoking", "P. Work", "Pressure", "Proteins", "Family"});
  std::set<uint8_t> computedMWorkPC = gmb.getCorrectPC(target);
  EXPECT_EQ(computedMWorkPC, trueMWorkPC);

  target = data.varIndex("P. Work");
  std::set<uint8_t> truePWorkPC = data.varIndices({"Smoking", "M. Work"});
  std::set<uint8_t> computedPWorkPC = gmb.getCorrectPC(target);
  EXPECT_EQ(computedPWorkPC, truePWorkPC);

  target = data.varIndex("Pressure");
  std::set<uint8_t> truePressurePC = data.varIndices({"Smoking", "M. Work", "Proteins"});
  std::set<uint8_t> computedPressurePC = gmb.getCorrectPC(target);
  EXPECT_EQ(computedPressurePC, truePressurePC);

  target = data.varIndex("Proteins");
  std::set<uint8_t> trueProteinsPC = data.varIndices({"Smoking", "M. Work", "Pressure"});
  std::set<uint8_t> computedProteinsPC = gmb.getCorrectPC(target);
  EXPECT_EQ(computedProteinsPC, trueProteinsPC);

  target = data.varIndex("Family");
  std::set<uint8_t> trueFamilyPC = data.varIndices({"M. Work"});
  std::set<uint8_t> computedFamilyPC = gmb.getCorrectPC(target);
  EXPECT_EQ(computedFamilyPC, trueFamilyPC);
}

TEST_F(AsiaTest, TopologicalDiscovery_MMPC) {
  MMPC<decltype(data), uint8_t> mm(data);

  uint8_t target = data.varIndex("asia");
  std::set<uint8_t> trueAsiaPC = data.varIndices({});
  std::set<uint8_t> computedAsiaPC = mm.getCorrectPC(target);
  EXPECT_EQ(computedAsiaPC, trueAsiaPC);

  target = data.varIndex("smoke");
  std::set<uint8_t> trueSmokePC = data.varIndices({"bronc", "lung"});
  std::set<uint8_t> computedSmokePC = mm.getCorrectPC(target);
  EXPECT_EQ(computedSmokePC, trueSmokePC);

  target = data.varIndex("tub");
  std::set<uint8_t> trueTubPC = data.varIndices({"either"});
  std::set<uint8_t> computedTubPC = mm.getCorrectPC(target);
  EXPECT_EQ(computedTubPC, trueTubPC);

  target = data.varIndex("lung");
  std::set<uint8_t> trueLungPC = data.varIndices({"either", "smoke"});
  std::set<uint8_t> computedLungPC = mm.getCorrectPC(target);
  EXPECT_EQ(computedLungPC, trueLungPC);

  target = data.varIndex("bronc");
  std::set<uint8_t> trueBroncPC = data.varIndices({"dysp", "smoke"});
  std::set<uint8_t> computedBroncPC = mm.getCorrectPC(target);
  EXPECT_EQ(computedBroncPC, trueBroncPC);

  target = data.varIndex("either");
  std::set<uint8_t> trueEitherPC = data.varIndices({"lung", "tub"});
  std::set<uint8_t> computedEitherPC = mm.getCorrectPC(target);
  EXPECT_EQ(computedEitherPC, trueEitherPC);

  target = data.varIndex("xray");
  std::set<uint8_t> trueXrayPC = data.varIndices({});
  std::set<uint8_t> computedXrayPC = mm.getCorrectPC(target);
  EXPECT_EQ(computedXrayPC, trueXrayPC);

  target = data.varIndex("dysp");
  std::set<uint8_t> trueDyspPC = data.varIndices({"bronc"});
  std::set<uint8_t> computedDyspPC = mm.getCorrectPC(target);
  EXPECT_EQ(computedDyspPC, trueDyspPC);
}

TEST_F(AsiaTest, TopologicalDiscovery_GetPC) {
  GetPC<decltype(data), uint8_t> gmb(data);

  uint8_t target = data.varIndex("asia");
  std::set<uint8_t> trueAsiaPC = data.varIndices({});
  std::set<uint8_t> computedAsiaPC = gmb.getCorrectPC(target);
  EXPECT_EQ(computedAsiaPC, trueAsiaPC);

  target = data.varIndex("smoke");
  std::set<uint8_t> trueSmokePC = data.varIndices({"bronc", "lung"});
  std::set<uint8_t> computedSmokePC = gmb.getCorrectPC(target);
  EXPECT_EQ(computedSmokePC, trueSmokePC);

  target = data.varIndex("tub");
  std::set<uint8_t> trueTubPC = data.varIndices({"either"});
  std::set<uint8_t> computedTubPC = gmb.getCorrectPC(target);
  EXPECT_EQ(computedTubPC, trueTubPC);

  target = data.varIndex("lung");
  std::set<uint8_t> trueLungPC = data.varIndices({"either", "smoke"});
  std::set<uint8_t> computedLungPC = gmb.getCorrectPC(target);
  EXPECT_EQ(computedLungPC, trueLungPC);

  target = data.varIndex("bronc");
  std::set<uint8_t> trueBroncPC = data.varIndices({"dysp", "smoke"});
  std::set<uint8_t> computedBroncPC = gmb.getCorrectPC(target);
  EXPECT_EQ(computedBroncPC, trueBroncPC);

  target = data.varIndex("either");
  std::set<uint8_t> trueEitherPC = data.varIndices({"lung", "tub"});
  std::set<uint8_t> computedEitherPC = gmb.getCorrectPC(target);
  EXPECT_EQ(computedEitherPC, trueEitherPC);

  target = data.varIndex("xray");
  std::set<uint8_t> trueXrayPC = data.varIndices({});
  std::set<uint8_t> computedXrayPC = gmb.getCorrectPC(target);
  EXPECT_EQ(computedXrayPC, trueXrayPC);

  target = data.varIndex("dysp");
  std::set<uint8_t> trueDyspPC = data.varIndices({"bronc"});
  std::set<uint8_t> computedDyspPC = gmb.getCorrectPC(target);
  EXPECT_EQ(computedDyspPC, trueDyspPC);
}

#endif // TEST_TOPOLOGICALDISCOVERY_HPP_
