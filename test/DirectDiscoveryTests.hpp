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

  target = data.varIndex("M. Work");
  std::set<uint8_t> trueMWorkMB = data.varIndices({"Smoking", "P. Work", "Pressure", "Proteins", "Family"});
  std::set<uint8_t> computedMWorkMB = gs.getMB(target);
  EXPECT_EQ(computedMWorkMB, trueMWorkMB);

  target = data.varIndex("P. Work");
  std::set<uint8_t> truePWorkMB = data.varIndices({"Smoking", "M. Work", "Pressure", "Proteins"});
  std::set<uint8_t> computedPWorkMB = gs.getMB(target);
  EXPECT_EQ(computedPWorkMB, truePWorkMB);

  target = data.varIndex("Pressure");
  std::set<uint8_t> truePressureMB = data.varIndices({"Smoking", "M. Work", "P. Work", "Proteins"});
  std::set<uint8_t> computedPressureMB = gs.getMB(target);
  EXPECT_EQ(computedPressureMB, truePressureMB);

  target = data.varIndex("Proteins");
  std::set<uint8_t> trueProteinsMB = data.varIndices({"Smoking", "M. Work", "P. Work", "Pressure"});
  std::set<uint8_t> computedProteinsMB = gs.getMB(target);
  EXPECT_EQ(computedProteinsMB, trueProteinsMB);

  target = data.varIndex("Family");
  std::set<uint8_t> trueFamilyMB = data.varIndices({"M. Work"});
  std::set<uint8_t> computedFamilyMB = gs.getMB(target);
  EXPECT_EQ(computedFamilyMB, trueFamilyMB);
}

TEST_F(CoronaryTest, DirectDiscovery_IAMB) {
  IAMB<decltype(data), uint8_t> ia(data);
  std::set<uint8_t> exclude;

  uint8_t target = data.varIndex("Smoking");
  std::set<uint8_t> trueSmokingMB = data.varIndices({"M. Work", "P. Work", "Pressure", "Proteins"});
  std::set<uint8_t> computedSmokingMB = ia.getMB(target);
  EXPECT_EQ(computedSmokingMB, trueSmokingMB);

  target = data.varIndex("M. Work");
  std::set<uint8_t> trueMWorkMB = data.varIndices({"Smoking", "P. Work", "Pressure", "Proteins", "Family"});
  std::set<uint8_t> computedMWorkMB = ia.getMB(target);
  EXPECT_EQ(computedMWorkMB, trueMWorkMB);

  target = data.varIndex("P. Work");
  std::set<uint8_t> truePWorkMB = data.varIndices({"Smoking", "M. Work", "Pressure", "Proteins"});
  std::set<uint8_t> computedPWorkMB = ia.getMB(target);
  EXPECT_EQ(computedPWorkMB, truePWorkMB);

  target = data.varIndex("Pressure");
  std::set<uint8_t> truePressureMB = data.varIndices({"Smoking", "M. Work", "P. Work", "Proteins"});
  std::set<uint8_t> computedPressureMB = ia.getMB(target);
  EXPECT_EQ(computedPressureMB, truePressureMB);

  target = data.varIndex("Proteins");
  std::set<uint8_t> trueProteinsMB = data.varIndices({"Smoking", "M. Work", "P. Work", "Pressure"});
  std::set<uint8_t> computedProteinsMB = ia.getMB(target);
  EXPECT_EQ(computedProteinsMB, trueProteinsMB);

  target = data.varIndex("Family");
  std::set<uint8_t> trueFamilyMB = data.varIndices({"M. Work"});
  std::set<uint8_t> computedFamilyMB = ia.getMB(target);
  EXPECT_EQ(computedFamilyMB, trueFamilyMB);
}

TEST_F(CoronaryTest, DirectDiscovery_InterIAMB) {
  InterIAMB<decltype(data), uint8_t> inter(data);
  std::set<uint8_t> exclude;

  uint8_t target = data.varIndex("Smoking");
  std::set<uint8_t> trueSmokingMB = data.varIndices({"M. Work", "P. Work", "Pressure", "Proteins"});
  std::set<uint8_t> computedSmokingMB = inter.getMB(target);
  EXPECT_EQ(computedSmokingMB, trueSmokingMB);

  target = data.varIndex("M. Work");
  std::set<uint8_t> trueMWorkMB = data.varIndices({"Smoking", "P. Work", "Pressure", "Proteins", "Family"});
  std::set<uint8_t> computedMWorkMB = inter.getMB(target);
  EXPECT_EQ(computedMWorkMB, trueMWorkMB);

  target = data.varIndex("P. Work");
  std::set<uint8_t> truePWorkMB = data.varIndices({"Smoking", "M. Work", "Pressure", "Proteins"});
  std::set<uint8_t> computedPWorkMB = inter.getMB(target);
  EXPECT_EQ(computedPWorkMB, truePWorkMB);

  target = data.varIndex("Pressure");
  std::set<uint8_t> truePressureMB = data.varIndices({"Smoking", "M. Work", "P. Work", "Proteins"});
  std::set<uint8_t> computedPressureMB = inter.getMB(target);
  EXPECT_EQ(computedPressureMB, truePressureMB);

  target = data.varIndex("Proteins");
  std::set<uint8_t> trueProteinsMB = data.varIndices({"Smoking", "M. Work", "P. Work", "Pressure"});
  std::set<uint8_t> computedProteinsMB = inter.getMB(target);
  EXPECT_EQ(computedProteinsMB, trueProteinsMB);

  target = data.varIndex("Family");
  std::set<uint8_t> trueFamilyMB = data.varIndices({"M. Work"});
  std::set<uint8_t> computedFamilyMB = inter.getMB(target);
  EXPECT_EQ(computedFamilyMB, trueFamilyMB);
}

#endif // TEST_DIRECTDISCOVERY_HPP_
