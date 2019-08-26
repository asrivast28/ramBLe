/**
 * @file TopologicalDiscoveryTests.hpp
 * @brief Unit tests for the Topological MB discovery algorithms.
 */
#ifndef TEST_TOPOLOGICALDISCOVERY_HPP_
#define TEST_TOPOLOGICALDISCOVERY_HPP_

#include "Common.hpp"
#include "TopologicalDiscovery.hpp"


TEST_F(CoronaryTest, TopologicalDiscovery_MMPC) {
  MMPC<decltype(data), uint8_t> mmpc(data);

  uint8_t target = data.varIndex("Smoking");
  std::set<uint8_t> trueSmokingPC = data.varIndices({"M. Work", "P. Work", "Pressure", "Proteins"});
  std::set<uint8_t> computedSmokingPC = mmpc.getCorrectPC(target);
  EXPECT_EQ(computedSmokingPC, trueSmokingPC);

  target = data.varIndex("M. Work");
  std::set<uint8_t> trueMWorkPC = data.varIndices({"Smoking", "P. Work", "Pressure", "Proteins", "Family"});
  std::set<uint8_t> computedMWorkPC = mmpc.getCorrectPC(target);
  EXPECT_EQ(computedMWorkPC, trueMWorkPC);

  target = data.varIndex("P. Work");
  std::set<uint8_t> truePWorkPC = data.varIndices({"Smoking", "M. Work"});
  std::set<uint8_t> computedPWorkPC = mmpc.getCorrectPC(target);
  EXPECT_EQ(computedPWorkPC, truePWorkPC);

  target = data.varIndex("Pressure");
  std::set<uint8_t> truePressurePC = data.varIndices({"Smoking", "M. Work", "Proteins"});
  std::set<uint8_t> computedPressurePC = mmpc.getCorrectPC(target);
  EXPECT_EQ(computedPressurePC, truePressurePC);

  target = data.varIndex("Proteins");
  std::set<uint8_t> trueProteinsPC = data.varIndices({"Smoking", "M. Work", "Pressure"});
  std::set<uint8_t> computedProteinsPC = mmpc.getCorrectPC(target);
  EXPECT_EQ(computedProteinsPC, trueProteinsPC);

  target = data.varIndex("Family");
  std::set<uint8_t> trueFamilyPC = data.varIndices({"M. Work"});
  std::set<uint8_t> computedFamilyPC = mmpc.getCorrectPC(target);
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

#endif // TEST_TOPOLOGICALDISCOVERY_HPP_
