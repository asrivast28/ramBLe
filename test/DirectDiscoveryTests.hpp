/**
 * @file DirectDiscoveryTests.hpp
 * @brief Unit tests for the Direct MB discovery algorithms.
 */
#ifndef TEST_DIRECTDISCOVERY_HPP_
#define TEST_DIRECTDISCOVERY_HPP_

#include "Common.hpp"
#include "DirectDiscovery.hpp"


TEST_F(CoronaryTest, DirectDiscovery_GSMB) {
  GSMB<decltype(data), uint8_t> gs(data);
  std::set<uint8_t> exclude;

  uint8_t target = data.varIndex("Smoking");
  std::set<uint8_t> trueSmokingMB = data.varIndices({"M. Work", "P. Work", "Pressure", "Proteins"});
  std::set<uint8_t> computedSmokingMB = gs.getMB(target);
  EXPECT_EQ(computedSmokingMB, trueSmokingMB);
  std::set<uint8_t> trueSmokingPC = data.varIndices({"M. Work", "P. Work", "Pressure", "Proteins"});
  std::set<uint8_t> computedSmokingPC = gs.getPC(target);
  EXPECT_EQ(computedSmokingPC, trueSmokingPC);

  target = data.varIndex("M. Work");
  std::set<uint8_t> trueMWorkMB = data.varIndices({"Smoking", "P. Work", "Pressure", "Proteins", "Family"});
  std::set<uint8_t> computedMWorkMB = gs.getMB(target);
  EXPECT_EQ(computedMWorkMB, trueMWorkMB);
  std::set<uint8_t> trueMWorkPC = data.varIndices({"Smoking", "P. Work", "Pressure", "Proteins", "Family"});
  std::set<uint8_t> computedMWorkPC = gs.getPC(target);
  EXPECT_EQ(computedMWorkPC, trueMWorkPC);

  target = data.varIndex("P. Work");
  std::set<uint8_t> truePWorkMB = data.varIndices({"Smoking", "M. Work", "Pressure", "Proteins"});
  std::set<uint8_t> computedPWorkMB = gs.getMB(target);
  EXPECT_EQ(computedPWorkMB, truePWorkMB);
  std::set<uint8_t> truePWorkPC = data.varIndices({"Smoking", "M. Work"});
  std::set<uint8_t> computedPWorkPC = gs.getPC(target);
  EXPECT_EQ(computedPWorkPC, truePWorkPC);

  target = data.varIndex("Pressure");
  std::set<uint8_t> truePressureMB = data.varIndices({"Smoking", "M. Work", "P. Work", "Proteins"});
  std::set<uint8_t> computedPressureMB = gs.getMB(target);
  EXPECT_EQ(computedPressureMB, truePressureMB);
  std::set<uint8_t> truePressurePC = data.varIndices({"Smoking", "M. Work", "Proteins"});
  std::set<uint8_t> computedPressurePC = gs.getPC(target);
  EXPECT_EQ(computedPressurePC, truePressurePC);

  target = data.varIndex("Proteins");
  std::set<uint8_t> trueProteinsMB = data.varIndices({"Smoking", "M. Work", "P. Work", "Pressure"});
  std::set<uint8_t> computedProteinsMB = gs.getMB(target);
  EXPECT_EQ(computedProteinsMB, trueProteinsMB);
  std::set<uint8_t> trueProteinsPC = data.varIndices({"Smoking", "M. Work", "Pressure"});
  std::set<uint8_t> computedProteinsPC = gs.getPC(target);
  EXPECT_EQ(computedProteinsPC, trueProteinsPC);

  target = data.varIndex("Family");
  std::set<uint8_t> trueFamilyMB = data.varIndices({"M. Work"});
  std::set<uint8_t> computedFamilyMB = gs.getMB(target);
  EXPECT_EQ(computedFamilyMB, trueFamilyMB);
  std::set<uint8_t> trueFamilyPC = data.varIndices({"M. Work"});
  std::set<uint8_t> computedFamilyPC = gs.getPC(target);
  EXPECT_EQ(computedFamilyPC, trueFamilyPC);
}

TEST_F(CoronaryTest, DirectDiscovery_IAMB) {
  IAMB<decltype(data), uint8_t> ia(data);
  std::set<uint8_t> exclude;

  uint8_t target = data.varIndex("Smoking");
  std::set<uint8_t> trueSmokingMB = data.varIndices({"M. Work", "P. Work", "Pressure", "Proteins"});
  std::set<uint8_t> computedSmokingMB = ia.getMB(target);
  EXPECT_EQ(computedSmokingMB, trueSmokingMB);
  std::set<uint8_t> trueSmokingPC = data.varIndices({"M. Work", "P. Work", "Pressure", "Proteins"});
  std::set<uint8_t> computedSmokingPC = ia.getPC(target);
  EXPECT_EQ(computedSmokingPC, trueSmokingPC);

  target = data.varIndex("M. Work");
  std::set<uint8_t> trueMWorkMB = data.varIndices({"Smoking", "P. Work", "Pressure", "Proteins", "Family"});
  std::set<uint8_t> computedMWorkMB = ia.getMB(target);
  EXPECT_EQ(computedMWorkMB, trueMWorkMB);
  std::set<uint8_t> trueMWorkPC = data.varIndices({"Smoking", "P. Work", "Pressure", "Proteins", "Family"});
  std::set<uint8_t> computedMWorkPC = ia.getPC(target);
  EXPECT_EQ(computedMWorkPC, trueMWorkPC);

  target = data.varIndex("P. Work");
  std::set<uint8_t> truePWorkMB = data.varIndices({"Smoking", "M. Work", "Pressure", "Proteins"});
  std::set<uint8_t> computedPWorkMB = ia.getMB(target);
  EXPECT_EQ(computedPWorkMB, truePWorkMB);
  std::set<uint8_t> truePWorkPC = data.varIndices({"Smoking", "M. Work"});
  std::set<uint8_t> computedPWorkPC = ia.getPC(target);
  EXPECT_EQ(computedPWorkPC, truePWorkPC);

  target = data.varIndex("Pressure");
  std::set<uint8_t> truePressureMB = data.varIndices({"Smoking", "M. Work", "P. Work", "Proteins"});
  std::set<uint8_t> computedPressureMB = ia.getMB(target);
  EXPECT_EQ(computedPressureMB, truePressureMB);
  std::set<uint8_t> truePressurePC = data.varIndices({"Smoking", "M. Work", "Proteins"});
  std::set<uint8_t> computedPressurePC = ia.getPC(target);
  EXPECT_EQ(computedPressurePC, truePressurePC);

  target = data.varIndex("Proteins");
  std::set<uint8_t> trueProteinsMB = data.varIndices({"Smoking", "M. Work", "P. Work", "Pressure"});
  std::set<uint8_t> computedProteinsMB = ia.getMB(target);
  EXPECT_EQ(computedProteinsMB, trueProteinsMB);
  std::set<uint8_t> trueProteinsPC = data.varIndices({"Smoking", "M. Work", "Pressure"});
  std::set<uint8_t> computedProteinsPC = ia.getPC(target);
  EXPECT_EQ(computedProteinsPC, trueProteinsPC);

  target = data.varIndex("Family");
  std::set<uint8_t> trueFamilyMB = data.varIndices({"M. Work"});
  std::set<uint8_t> computedFamilyMB = ia.getMB(target);
  EXPECT_EQ(computedFamilyMB, trueFamilyMB);
  std::set<uint8_t> trueFamilyPC = data.varIndices({"M. Work"});
  std::set<uint8_t> computedFamilyPC = ia.getPC(target);
  EXPECT_EQ(computedFamilyPC, trueFamilyPC);
}

TEST_F(CoronaryTest, DirectDiscovery_InterIAMB) {
  InterIAMB<decltype(data), uint8_t> inter(data);
  std::set<uint8_t> exclude;

  uint8_t target = data.varIndex("Smoking");
  std::set<uint8_t> trueSmokingMB = data.varIndices({"M. Work", "P. Work", "Pressure", "Proteins"});
  std::set<uint8_t> computedSmokingMB = inter.getMB(target);
  EXPECT_EQ(computedSmokingMB, trueSmokingMB);
  std::set<uint8_t> trueSmokingPC = data.varIndices({"M. Work", "P. Work", "Pressure", "Proteins"});
  std::set<uint8_t> computedSmokingPC = inter.getPC(target);
  EXPECT_EQ(computedSmokingPC, trueSmokingPC);

  target = data.varIndex("M. Work");
  std::set<uint8_t> trueMWorkMB = data.varIndices({"Smoking", "P. Work", "Pressure", "Proteins", "Family"});
  std::set<uint8_t> computedMWorkMB = inter.getMB(target);
  EXPECT_EQ(computedMWorkMB, trueMWorkMB);
  std::set<uint8_t> trueMWorkPC = data.varIndices({"Smoking", "P. Work", "Pressure", "Proteins", "Family"});
  std::set<uint8_t> computedMWorkPC = inter.getPC(target);
  EXPECT_EQ(computedMWorkPC, trueMWorkPC);

  target = data.varIndex("P. Work");
  std::set<uint8_t> truePWorkMB = data.varIndices({"Smoking", "M. Work", "Pressure", "Proteins"});
  std::set<uint8_t> computedPWorkMB = inter.getMB(target);
  EXPECT_EQ(computedPWorkMB, truePWorkMB);
  std::set<uint8_t> truePWorkPC = data.varIndices({"Smoking", "M. Work"});
  std::set<uint8_t> computedPWorkPC = inter.getPC(target);
  EXPECT_EQ(computedPWorkPC, truePWorkPC);

  target = data.varIndex("Pressure");
  std::set<uint8_t> truePressureMB = data.varIndices({"Smoking", "M. Work", "P. Work", "Proteins"});
  std::set<uint8_t> computedPressureMB = inter.getMB(target);
  EXPECT_EQ(computedPressureMB, truePressureMB);
  std::set<uint8_t> truePressurePC = data.varIndices({"Smoking", "M. Work", "Proteins"});
  std::set<uint8_t> computedPressurePC = inter.getPC(target);
  EXPECT_EQ(computedPressurePC, truePressurePC);

  target = data.varIndex("Proteins");
  std::set<uint8_t> trueProteinsMB = data.varIndices({"Smoking", "M. Work", "P. Work", "Pressure"});
  std::set<uint8_t> computedProteinsMB = inter.getMB(target);
  EXPECT_EQ(computedProteinsMB, trueProteinsMB);
  std::set<uint8_t> trueProteinsPC = data.varIndices({"Smoking", "M. Work", "Pressure"});
  std::set<uint8_t> computedProteinsPC = inter.getPC(target);
  EXPECT_EQ(computedProteinsPC, trueProteinsPC);

  target = data.varIndex("Family");
  std::set<uint8_t> trueFamilyMB = data.varIndices({"M. Work"});
  std::set<uint8_t> computedFamilyMB = inter.getMB(target);
  EXPECT_EQ(computedFamilyMB, trueFamilyMB);
  std::set<uint8_t> trueFamilyPC = data.varIndices({"M. Work"});
  std::set<uint8_t> computedFamilyPC = inter.getPC(target);
  EXPECT_EQ(computedFamilyPC, trueFamilyPC);
}

TEST_F(AsiaTest, DirectDiscovery_GSMB) {
  GSMB<decltype(data), uint8_t> gs(data);

  uint8_t target = data.varIndex("asia");
  std::set<uint8_t> trueAsiaMB = data.varIndices({});
  std::set<uint8_t> computedAsiaMB = gs.getMB(target);
  EXPECT_EQ(computedAsiaMB, trueAsiaMB);
  std::set<uint8_t> trueAsiaPC = data.varIndices({});
  std::set<uint8_t> computedAsiaPC = gs.getPC(target);
  EXPECT_EQ(computedAsiaPC, trueAsiaPC);

  target = data.varIndex("smoke");
  std::set<uint8_t> trueSmokeMB = data.varIndices({"bronc"});
  std::set<uint8_t> computedSmokeMB = gs.getMB(target);
  EXPECT_EQ(computedSmokeMB, trueSmokeMB);
  std::set<uint8_t> trueSmokePC = data.varIndices({"bronc"});
  std::set<uint8_t> computedSmokePC = gs.getPC(target);
  EXPECT_EQ(computedSmokePC, trueSmokePC);

  target = data.varIndex("tub");
  std::set<uint8_t> trueTubMB = data.varIndices({"either", "lung"});
  std::set<uint8_t> computedTubMB = gs.getMB(target);
  EXPECT_EQ(computedTubMB, trueTubMB);
  std::set<uint8_t> trueTubPC = data.varIndices({"either"});
  std::set<uint8_t> computedTubPC = gs.getPC(target);
  EXPECT_EQ(computedTubPC, trueTubPC);

  target = data.varIndex("lung");
  std::set<uint8_t> trueLungMB = data.varIndices({"either", "tub"});
  std::set<uint8_t> computedLungMB = gs.getMB(target);
  EXPECT_EQ(computedLungMB, trueLungMB);
  std::set<uint8_t> trueLungPC = data.varIndices({"either"});
  std::set<uint8_t> computedLungPC = gs.getPC(target);
  EXPECT_EQ(computedLungPC, trueLungPC);

  target = data.varIndex("bronc");
  std::set<uint8_t> trueBroncMB = data.varIndices({"dysp", "smoke"});
  std::set<uint8_t> computedBroncMB = gs.getMB(target);
  EXPECT_EQ(computedBroncMB, trueBroncMB);
  std::set<uint8_t> trueBroncPC = data.varIndices({"dysp", "smoke"});
  std::set<uint8_t> computedBroncPC = gs.getPC(target);
  EXPECT_EQ(computedBroncPC, trueBroncPC);

  target = data.varIndex("either");
  std::set<uint8_t> trueEitherMB = data.varIndices({"lung", "tub"});
  std::set<uint8_t> computedEitherMB = gs.getMB(target);
  EXPECT_EQ(computedEitherMB, trueEitherMB);
  std::set<uint8_t> trueEitherPC = data.varIndices({"lung", "tub"});
  std::set<uint8_t> computedEitherPC = gs.getPC(target);
  EXPECT_EQ(computedEitherPC, trueEitherPC);

  target = data.varIndex("xray");
  std::set<uint8_t> trueXrayMB = data.varIndices({});
  std::set<uint8_t> computedXrayMB = gs.getMB(target);
  EXPECT_EQ(computedXrayMB, trueXrayMB);
  std::set<uint8_t> trueXrayPC = data.varIndices({});
  std::set<uint8_t> computedXrayPC = gs.getPC(target);
  EXPECT_EQ(computedXrayPC, trueXrayPC);

  target = data.varIndex("dysp");
  std::set<uint8_t> trueDyspMB = data.varIndices({"bronc"});
  std::set<uint8_t> computedDyspMB = gs.getMB(target);
  EXPECT_EQ(computedDyspMB, trueDyspMB);
  std::set<uint8_t> trueDyspPC = data.varIndices({"bronc"});
  std::set<uint8_t> computedDyspPC = gs.getPC(target);
  EXPECT_EQ(computedDyspPC, trueDyspPC);
}

TEST_F(AsiaTest, DirectDiscovery_IAMB) {
  IAMB<decltype(data), uint8_t> ia(data);

  uint8_t target = data.varIndex("asia");
  std::set<uint8_t> trueAsiaMB = data.varIndices({});
  std::set<uint8_t> computedAsiaMB = ia.getMB(target);
  EXPECT_EQ(computedAsiaMB, trueAsiaMB);
  std::set<uint8_t> trueAsiaPC = data.varIndices({});
  std::set<uint8_t> computedAsiaPC = ia.getPC(target);
  EXPECT_EQ(computedAsiaPC, trueAsiaPC);

  target = data.varIndex("smoke");
  std::set<uint8_t> trueSmokeMB = data.varIndices({"bronc"});
  std::set<uint8_t> computedSmokeMB = ia.getMB(target);
  EXPECT_EQ(computedSmokeMB, trueSmokeMB);
  std::set<uint8_t> trueSmokePC = data.varIndices({"bronc"});
  std::set<uint8_t> computedSmokePC = ia.getPC(target);
  EXPECT_EQ(computedSmokePC, trueSmokePC);

  target = data.varIndex("tub");
  std::set<uint8_t> trueTubMB = data.varIndices({"either", "lung"});
  std::set<uint8_t> computedTubMB = ia.getMB(target);
  EXPECT_EQ(computedTubMB, trueTubMB);
  std::set<uint8_t> trueTubPC = data.varIndices({"either"});
  std::set<uint8_t> computedTubPC = ia.getPC(target);
  EXPECT_EQ(computedTubPC, trueTubPC);

  target = data.varIndex("lung");
  std::set<uint8_t> trueLungMB = data.varIndices({"either", "tub"});
  std::set<uint8_t> computedLungMB = ia.getMB(target);
  EXPECT_EQ(computedLungMB, trueLungMB);
  std::set<uint8_t> trueLungPC = data.varIndices({"either"});
  std::set<uint8_t> computedLungPC = ia.getPC(target);
  EXPECT_EQ(computedLungPC, trueLungPC);

  target = data.varIndex("bronc");
  std::set<uint8_t> trueBroncMB = data.varIndices({"dysp", "smoke"});
  std::set<uint8_t> computedBroncMB = ia.getMB(target);
  EXPECT_EQ(computedBroncMB, trueBroncMB);
  std::set<uint8_t> trueBroncPC = data.varIndices({"dysp", "smoke"});
  std::set<uint8_t> computedBroncPC = ia.getPC(target);
  EXPECT_EQ(computedBroncPC, trueBroncPC);

  target = data.varIndex("either");
  std::set<uint8_t> trueEitherMB = data.varIndices({"lung", "tub"});
  std::set<uint8_t> computedEitherMB = ia.getMB(target);
  EXPECT_EQ(computedEitherMB, trueEitherMB);
  std::set<uint8_t> trueEitherPC = data.varIndices({"lung", "tub"});
  std::set<uint8_t> computedEitherPC = ia.getPC(target);
  EXPECT_EQ(computedEitherPC, trueEitherPC);

  target = data.varIndex("xray");
  std::set<uint8_t> trueXrayMB = data.varIndices({});
  std::set<uint8_t> computedXrayMB = ia.getMB(target);
  EXPECT_EQ(computedXrayMB, trueXrayMB);
  std::set<uint8_t> trueXrayPC = data.varIndices({});
  std::set<uint8_t> computedXrayPC = ia.getPC(target);
  EXPECT_EQ(computedXrayPC, trueXrayPC);

  target = data.varIndex("dysp");
  std::set<uint8_t> trueDyspMB = data.varIndices({"bronc"});
  std::set<uint8_t> computedDyspMB = ia.getMB(target);
  EXPECT_EQ(computedDyspMB, trueDyspMB);
  std::set<uint8_t> trueDyspPC = data.varIndices({"bronc"});
  std::set<uint8_t> computedDyspPC = ia.getPC(target);
  EXPECT_EQ(computedDyspPC, trueDyspPC);
}

TEST_F(AsiaTest, DirectDiscovery_InterIAMB) {
  InterIAMB<decltype(data), uint8_t> inter(data);

  uint8_t target = data.varIndex("asia");
  std::set<uint8_t> trueAsiaMB = data.varIndices({});
  std::set<uint8_t> computedAsiaMB = inter.getMB(target);
  EXPECT_EQ(computedAsiaMB, trueAsiaMB);
  std::set<uint8_t> trueAsiaPC = data.varIndices({});
  std::set<uint8_t> computedAsiaPC = inter.getPC(target);
  EXPECT_EQ(computedAsiaPC, trueAsiaPC);

  target = data.varIndex("smoke");
  std::set<uint8_t> trueSmokeMB = data.varIndices({"bronc"});
  std::set<uint8_t> computedSmokeMB = inter.getMB(target);
  EXPECT_EQ(computedSmokeMB, trueSmokeMB);
  std::set<uint8_t> trueSmokePC = data.varIndices({"bronc"});
  std::set<uint8_t> computedSmokePC = inter.getPC(target);
  EXPECT_EQ(computedSmokePC, trueSmokePC);

  target = data.varIndex("tub");
  std::set<uint8_t> trueTubMB = data.varIndices({"either", "lung"});
  std::set<uint8_t> computedTubMB = inter.getMB(target);
  EXPECT_EQ(computedTubMB, trueTubMB);
  std::set<uint8_t> trueTubPC = data.varIndices({"either"});
  std::set<uint8_t> computedTubPC = inter.getPC(target);
  EXPECT_EQ(computedTubPC, trueTubPC);

  target = data.varIndex("lung");
  std::set<uint8_t> trueLungMB = data.varIndices({"either", "tub"});
  std::set<uint8_t> computedLungMB = inter.getMB(target);
  EXPECT_EQ(computedLungMB, trueLungMB);
  std::set<uint8_t> trueLungPC = data.varIndices({"either"});
  std::set<uint8_t> computedLungPC = inter.getPC(target);
  EXPECT_EQ(computedLungPC, trueLungPC);

  target = data.varIndex("bronc");
  std::set<uint8_t> trueBroncMB = data.varIndices({"dysp", "smoke"});
  std::set<uint8_t> computedBroncMB = inter.getMB(target);
  EXPECT_EQ(computedBroncMB, trueBroncMB);
  std::set<uint8_t> trueBroncPC = data.varIndices({"dysp", "smoke"});
  std::set<uint8_t> computedBroncPC = inter.getPC(target);
  EXPECT_EQ(computedBroncPC, trueBroncPC);

  target = data.varIndex("either");
  std::set<uint8_t> trueEitherMB = data.varIndices({"lung", "tub"});
  std::set<uint8_t> computedEitherMB = inter.getMB(target);
  EXPECT_EQ(computedEitherMB, trueEitherMB);
  std::set<uint8_t> trueEitherPC = data.varIndices({"lung", "tub"});
  std::set<uint8_t> computedEitherPC = inter.getPC(target);
  EXPECT_EQ(computedEitherPC, trueEitherPC);

  target = data.varIndex("xray");
  std::set<uint8_t> trueXrayMB = data.varIndices({});
  std::set<uint8_t> computedXrayMB = inter.getMB(target);
  EXPECT_EQ(computedXrayMB, trueXrayMB);
  std::set<uint8_t> trueXrayPC = data.varIndices({});
  std::set<uint8_t> computedXrayPC = inter.getPC(target);
  EXPECT_EQ(computedXrayPC, trueXrayPC);

  target = data.varIndex("dysp");
  std::set<uint8_t> trueDyspMB = data.varIndices({"bronc"});
  std::set<uint8_t> computedDyspMB = inter.getMB(target);
  EXPECT_EQ(computedDyspMB, trueDyspMB);
  std::set<uint8_t> trueDyspPC = data.varIndices({"bronc"});
  std::set<uint8_t> computedDyspPC = inter.getPC(target);
  EXPECT_EQ(computedDyspPC, trueDyspPC);
}

#endif // TEST_DIRECTDISCOVERY_HPP_
