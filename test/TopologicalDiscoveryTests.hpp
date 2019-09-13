/**
 * @file TopologicalDiscoveryTests.hpp
 * @brief Unit tests for the Topological MB discovery algorithms.
 */
#ifndef TEST_TOPOLOGICALDISCOVERY_HPP_
#define TEST_TOPOLOGICALDISCOVERY_HPP_

#include "Common.hpp"
#include "TopologicalDiscovery.hpp"


TEST_F(CoronaryTest, TopologicalDiscovery_MMPC) {
  MMPC<decltype(data), uint8_t, UintSet<uint8_t>> mm(data);

  uint8_t target = data.varIndex("Smoking");
  UintSet<uint8_t> trueSmokingPC = data.varIndices<UintSet<uint8_t>>({"M. Work", "P. Work", "Pressure", "Proteins"});
  UintSet<uint8_t> computedSmokingPC = mm.getPC(target);
  EXPECT_EQ(computedSmokingPC, trueSmokingPC);

  target = data.varIndex("M. Work");
  UintSet<uint8_t> trueMWorkPC = data.varIndices<UintSet<uint8_t>>({"Smoking", "P. Work", "Pressure", "Proteins", "Family"});
  UintSet<uint8_t> computedMWorkPC = mm.getPC(target);
  EXPECT_EQ(computedMWorkPC, trueMWorkPC);

  target = data.varIndex("P. Work");
  UintSet<uint8_t> truePWorkPC = data.varIndices<UintSet<uint8_t>>({"Smoking", "M. Work"});
  UintSet<uint8_t> computedPWorkPC = mm.getPC(target);
  EXPECT_EQ(computedPWorkPC, truePWorkPC);

  target = data.varIndex("Pressure");
  UintSet<uint8_t> truePressurePC = data.varIndices<UintSet<uint8_t>>({"Smoking", "M. Work", "Proteins"});
  UintSet<uint8_t> computedPressurePC = mm.getPC(target);
  EXPECT_EQ(computedPressurePC, truePressurePC);

  target = data.varIndex("Proteins");
  UintSet<uint8_t> trueProteinsPC = data.varIndices<UintSet<uint8_t>>({"Smoking", "M. Work", "Pressure"});
  UintSet<uint8_t> computedProteinsPC = mm.getPC(target);
  EXPECT_EQ(computedProteinsPC, trueProteinsPC);

  target = data.varIndex("Family");
  UintSet<uint8_t> trueFamilyPC = data.varIndices<UintSet<uint8_t>>({"M. Work"});
  UintSet<uint8_t> computedFamilyPC = mm.getPC(target);
  EXPECT_EQ(computedFamilyPC, trueFamilyPC);
}

TEST_F(CoronaryTest, TopologicalDiscovery_SemiInterleavedHITON) {
  SemiInterleavedHITON<decltype(data), uint8_t, UintSet<uint8_t>> sih(data);

  uint8_t target = data.varIndex("Smoking");
  UintSet<uint8_t> trueSmokingPC = data.varIndices<UintSet<uint8_t>>({"M. Work", "P. Work", "Pressure", "Proteins"});
  UintSet<uint8_t> computedSmokingPC = sih.getPC(target);
  EXPECT_EQ(computedSmokingPC, trueSmokingPC);

  target = data.varIndex("M. Work");
  UintSet<uint8_t> trueMWorkPC = data.varIndices<UintSet<uint8_t>>({"Smoking", "P. Work", "Pressure", "Proteins", "Family"});
  UintSet<uint8_t> computedMWorkPC = sih.getPC(target);
  EXPECT_EQ(computedMWorkPC, trueMWorkPC);

  target = data.varIndex("P. Work");
  UintSet<uint8_t> truePWorkPC = data.varIndices<UintSet<uint8_t>>({"Smoking", "M. Work"});
  UintSet<uint8_t> computedPWorkPC = sih.getPC(target);
  EXPECT_EQ(computedPWorkPC, truePWorkPC);

  target = data.varIndex("Pressure");
  UintSet<uint8_t> truePressurePC = data.varIndices<UintSet<uint8_t>>({"Smoking", "M. Work", "Proteins"});
  UintSet<uint8_t> computedPressurePC = sih.getPC(target);
  EXPECT_EQ(computedPressurePC, truePressurePC);

  target = data.varIndex("Proteins");
  UintSet<uint8_t> trueProteinsPC = data.varIndices<UintSet<uint8_t>>({"Smoking", "M. Work", "Pressure"});
  UintSet<uint8_t> computedProteinsPC = sih.getPC(target);
  EXPECT_EQ(computedProteinsPC, trueProteinsPC);

  target = data.varIndex("Family");
  UintSet<uint8_t> trueFamilyPC = data.varIndices<UintSet<uint8_t>>({"M. Work"});
  UintSet<uint8_t> computedFamilyPC = sih.getPC(target);
  EXPECT_EQ(computedFamilyPC, trueFamilyPC);
}

TEST_F(CoronaryTest, TopologicalDiscovery_GetPC) {
  GetPC<decltype(data), uint8_t, UintSet<uint8_t>> gmb(data);

  uint8_t target = data.varIndex("Smoking");
  UintSet<uint8_t> trueSmokingPC = data.varIndices<UintSet<uint8_t>>({"M. Work", "P. Work", "Pressure", "Proteins"});
  UintSet<uint8_t> computedSmokingPC = gmb.getPC(target);
  EXPECT_EQ(computedSmokingPC, trueSmokingPC);

  target = data.varIndex("M. Work");
  UintSet<uint8_t> trueMWorkPC = data.varIndices<UintSet<uint8_t>>({"Smoking", "P. Work", "Pressure", "Proteins", "Family"});
  UintSet<uint8_t> computedMWorkPC = gmb.getPC(target);
  EXPECT_EQ(computedMWorkPC, trueMWorkPC);

  target = data.varIndex("P. Work");
  UintSet<uint8_t> truePWorkPC = data.varIndices<UintSet<uint8_t>>({"Smoking", "M. Work"});
  UintSet<uint8_t> computedPWorkPC = gmb.getPC(target);
  EXPECT_EQ(computedPWorkPC, truePWorkPC);

  target = data.varIndex("Pressure");
  UintSet<uint8_t> truePressurePC = data.varIndices<UintSet<uint8_t>>({"Smoking", "M. Work", "Proteins"});
  UintSet<uint8_t> computedPressurePC = gmb.getPC(target);
  EXPECT_EQ(computedPressurePC, truePressurePC);

  target = data.varIndex("Proteins");
  UintSet<uint8_t> trueProteinsPC = data.varIndices<UintSet<uint8_t>>({"Smoking", "M. Work", "Pressure"});
  UintSet<uint8_t> computedProteinsPC = gmb.getPC(target);
  EXPECT_EQ(computedProteinsPC, trueProteinsPC);

  target = data.varIndex("Family");
  UintSet<uint8_t> trueFamilyPC = data.varIndices<UintSet<uint8_t>>({"M. Work"});
  UintSet<uint8_t> computedFamilyPC = gmb.getPC(target);
  EXPECT_EQ(computedFamilyPC, trueFamilyPC);
}

TEST_F(AsiaTest, TopologicalDiscovery_MMPC) {
  MMPC<decltype(data), uint8_t, UintSet<uint8_t>> mm(data);

  uint8_t target = data.varIndex("asia");
  UintSet<uint8_t> trueAsiaPC = data.varIndices<UintSet<uint8_t>>({});
  UintSet<uint8_t> computedAsiaPC = mm.getPC(target);
  EXPECT_EQ(computedAsiaPC, trueAsiaPC);

  target = data.varIndex("smoke");
  UintSet<uint8_t> trueSmokePC = data.varIndices<UintSet<uint8_t>>({"bronc", "lung"});
  UintSet<uint8_t> computedSmokePC = mm.getPC(target);
  EXPECT_EQ(computedSmokePC, trueSmokePC);

  target = data.varIndex("tub");
  UintSet<uint8_t> trueTubPC = data.varIndices<UintSet<uint8_t>>({"either"});
  UintSet<uint8_t> computedTubPC = mm.getPC(target);
  EXPECT_EQ(computedTubPC, trueTubPC);

  target = data.varIndex("lung");
  UintSet<uint8_t> trueLungPC = data.varIndices<UintSet<uint8_t>>({"either", "smoke"});
  UintSet<uint8_t> computedLungPC = mm.getPC(target);
  EXPECT_EQ(computedLungPC, trueLungPC);

  target = data.varIndex("bronc");
  UintSet<uint8_t> trueBroncPC = data.varIndices<UintSet<uint8_t>>({"dysp", "smoke"});
  UintSet<uint8_t> computedBroncPC = mm.getPC(target);
  EXPECT_EQ(computedBroncPC, trueBroncPC);

  target = data.varIndex("either");
  UintSet<uint8_t> trueEitherPC = data.varIndices<UintSet<uint8_t>>({"lung", "tub"});
  UintSet<uint8_t> computedEitherPC = mm.getPC(target);
  EXPECT_EQ(computedEitherPC, trueEitherPC);

  target = data.varIndex("xray");
  UintSet<uint8_t> trueXrayPC = data.varIndices<UintSet<uint8_t>>({});
  UintSet<uint8_t> computedXrayPC = mm.getPC(target);
  EXPECT_EQ(computedXrayPC, trueXrayPC);

  target = data.varIndex("dysp");
  UintSet<uint8_t> trueDyspPC = data.varIndices<UintSet<uint8_t>>({"bronc"});
  UintSet<uint8_t> computedDyspPC = mm.getPC(target);
  EXPECT_EQ(computedDyspPC, trueDyspPC);
}

TEST_F(AsiaTest, TopologicalDiscovery_SemiInterleavedHITON) {
  SemiInterleavedHITON<decltype(data), uint8_t, UintSet<uint8_t>> sih(data);

  uint8_t target = data.varIndex("asia");
  UintSet<uint8_t> trueAsiaPC = data.varIndices<UintSet<uint8_t>>({});
  UintSet<uint8_t> computedAsiaPC = sih.getPC(target);
  EXPECT_EQ(computedAsiaPC, trueAsiaPC);

  target = data.varIndex("smoke");
  UintSet<uint8_t> trueSmokePC = data.varIndices<UintSet<uint8_t>>({"bronc", "lung"});
  UintSet<uint8_t> computedSmokePC = sih.getPC(target);
  EXPECT_EQ(computedSmokePC, trueSmokePC);

  target = data.varIndex("tub");
  UintSet<uint8_t> trueTubPC = data.varIndices<UintSet<uint8_t>>({"either"});
  UintSet<uint8_t> computedTubPC = sih.getPC(target);
  EXPECT_EQ(computedTubPC, trueTubPC);

  target = data.varIndex("lung");
  UintSet<uint8_t> trueLungPC = data.varIndices<UintSet<uint8_t>>({"either", "smoke"});
  UintSet<uint8_t> computedLungPC = sih.getPC(target);
  EXPECT_EQ(computedLungPC, trueLungPC);

  target = data.varIndex("bronc");
  UintSet<uint8_t> trueBroncPC = data.varIndices<UintSet<uint8_t>>({"dysp", "smoke"});
  UintSet<uint8_t> computedBroncPC = sih.getPC(target);
  EXPECT_EQ(computedBroncPC, trueBroncPC);

  target = data.varIndex("either");
  UintSet<uint8_t> trueEitherPC = data.varIndices<UintSet<uint8_t>>({"lung", "tub"});
  UintSet<uint8_t> computedEitherPC = sih.getPC(target);
  EXPECT_EQ(computedEitherPC, trueEitherPC);

  target = data.varIndex("xray");
  UintSet<uint8_t> trueXrayPC = data.varIndices<UintSet<uint8_t>>({});
  UintSet<uint8_t> computedXrayPC = sih.getPC(target);
  EXPECT_EQ(computedXrayPC, trueXrayPC);

  target = data.varIndex("dysp");
  UintSet<uint8_t> trueDyspPC = data.varIndices<UintSet<uint8_t>>({"bronc"});
  UintSet<uint8_t> computedDyspPC = sih.getPC(target);
  EXPECT_EQ(computedDyspPC, trueDyspPC);
}

TEST_F(AsiaTest, TopologicalDiscovery_GetPC) {
  GetPC<decltype(data), uint8_t, UintSet<uint8_t>> gmb(data);

  uint8_t target = data.varIndex("asia");
  UintSet<uint8_t> trueAsiaPC = data.varIndices<UintSet<uint8_t>>({});
  UintSet<uint8_t> computedAsiaPC = gmb.getPC(target);
  EXPECT_EQ(computedAsiaPC, trueAsiaPC);

  target = data.varIndex("smoke");
  UintSet<uint8_t> trueSmokePC = data.varIndices<UintSet<uint8_t>>({"bronc", "lung"});
  UintSet<uint8_t> computedSmokePC = gmb.getPC(target);
  EXPECT_EQ(computedSmokePC, trueSmokePC);

  target = data.varIndex("tub");
  UintSet<uint8_t> trueTubPC = data.varIndices<UintSet<uint8_t>>({"either"});
  UintSet<uint8_t> computedTubPC = gmb.getPC(target);
  EXPECT_EQ(computedTubPC, trueTubPC);

  target = data.varIndex("lung");
  UintSet<uint8_t> trueLungPC = data.varIndices<UintSet<uint8_t>>({"either", "smoke"});
  UintSet<uint8_t> computedLungPC = gmb.getPC(target);
  EXPECT_EQ(computedLungPC, trueLungPC);

  target = data.varIndex("bronc");
  UintSet<uint8_t> trueBroncPC = data.varIndices<UintSet<uint8_t>>({"dysp", "smoke"});
  UintSet<uint8_t> computedBroncPC = gmb.getPC(target);
  EXPECT_EQ(computedBroncPC, trueBroncPC);

  target = data.varIndex("either");
  UintSet<uint8_t> trueEitherPC = data.varIndices<UintSet<uint8_t>>({"lung", "tub"});
  UintSet<uint8_t> computedEitherPC = gmb.getPC(target);
  EXPECT_EQ(computedEitherPC, trueEitherPC);

  target = data.varIndex("xray");
  UintSet<uint8_t> trueXrayPC = data.varIndices<UintSet<uint8_t>>({});
  UintSet<uint8_t> computedXrayPC = gmb.getPC(target);
  EXPECT_EQ(computedXrayPC, trueXrayPC);

  target = data.varIndex("dysp");
  UintSet<uint8_t> trueDyspPC = data.varIndices<UintSet<uint8_t>>({"bronc"});
  UintSet<uint8_t> computedDyspPC = gmb.getPC(target);
  EXPECT_EQ(computedDyspPC, trueDyspPC);
}

#endif // TEST_TOPOLOGICALDISCOVERY_HPP_
