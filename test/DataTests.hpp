/**
 * @file DataTests.hpp
 * @brief Unit tests for the dataset queries.
 */
#ifndef TEST_DATA_HPP_
#define TEST_DATA_HPP_

#include "Common.hpp"


TEST_F(NeapolitanTest, MarginalPValue) {
  double neapolitan1 = data[0].pValue(0, 1);
  EXPECT_NEAR(neapolitan1, 0.001, 0.0002);

  double neapolitan2 = data[1].pValue(0, 1);
  EXPECT_NEAR(neapolitan2, 1, 0.0001);

  double neapolitan3 = data[2].pValue(0, 1);
  EXPECT_NEAR(neapolitan3, 0.46, 0.002);
}

TEST_F(LizardsTest, MarginalPValue) {
  double spec_diam  = data.pValue(0, 1);
  EXPECT_NEAR(spec_diam, 0.0003845, 0.0001);

  double spec_hght = data.pValue(0, 2);
  EXPECT_NEAR(spec_hght, 0.001257, 0.0001);

  double diam_hght = data.pValue(1, 2);
  EXPECT_NEAR(diam_hght, 0.4357, 0.0001);
}

TEST_F(LizardsTest, ConditionalPValue) {
  std::set<uint8_t> given{2};
  double spec_diam  = data.pValue(0, 1, given);
  EXPECT_NEAR(spec_diam, 0.0009009, 0.0001);

  given = std::set<uint8_t>{1};
  double spec_hght = data.pValue(0, 2, given);
  EXPECT_NEAR(spec_hght, 0.002708, 0.0001);

  given = std::set<uint8_t>{0};
  double diam_hght = data.pValue(1, 2, given);
  EXPECT_NEAR(diam_hght, 0.3632, 0.0001);
}


TEST_F(CoronaryTest, MarginalPValue) {
  double mwork_pwork  = data.pValue(1, 2);
  EXPECT_NEAR(mwork_pwork, 0, 0.0001);

  double pwork_press = data.pValue(2, 3);
  EXPECT_NEAR(pwork_press, 0.7574, 0.0001);

  double mwork_press = data.pValue(1, 3);
  EXPECT_NEAR(mwork_press, 0.0004921, 0.0001);

  double smoke_press = data.pValue(0, 3);
  EXPECT_NEAR(smoke_press, 0.0008954, 0.0001);

  double press_family = data.pValue(3, 5);
  EXPECT_NEAR(press_family, 0.2891, 0.0001);
}

TEST_F(CoronaryTest, ConditionalPValue) {
  std::set<uint8_t> given{0};
  double press_family  = data.pValue(3, 5, given);
  EXPECT_NEAR(press_family, 0.3172, 0.0001);

  given.insert(1);
  press_family = data.pValue(3, 5, given);
  EXPECT_NEAR(press_family, 0.2186, 0.0001);

  given.insert(2);
  press_family  = data.pValue(3, 5, given);
  EXPECT_NEAR(press_family, 0.2452, 0.0001);

  given.insert(4);
  press_family  = data.pValue(3, 5, given);
  EXPECT_NEAR(press_family, 0.0002556, 0.0001);
}

#endif // TEST_DATA_HPP_
