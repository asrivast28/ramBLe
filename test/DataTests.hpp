/**
 * @file data->hpp
 * @brief Unit tests for the dataset queries.
 */
#ifndef TEST_DATA_HPP_
#define TEST_DATA_HPP_

#include "Environment.hpp"
#include "CTCounter.hpp"
#include "DiscreteData.hpp"
#include "UintSet.hpp"


// All the different counter implementations
using AllCounters = testing::Types<CTCounter<>>;

template <typename Counter>
class NeapolitanData : public testing::Test {
protected:
  void
  SetUp() override {
    auto n = NeapolitanEnvironment::n;
    auto m = NeapolitanEnvironment::m;
    auto reader = NeapolitanEnvironment::reader;
    data.resize(reader.size());
    for (auto i = 0u; i < data.size(); ++i) {
      auto counter = Counter::create(n, m, std::begin(reader[i]->data()));
      data[i] = new DiscreteData<Counter, uint8_t>(counter, reader[i]->varNames());
    }
  }

  void
  TearDown() override {
    for (auto i = 0u; i < data.size(); ++i) {
      delete data[i];
    }
  }

  std::vector<DiscreteData<Counter, uint8_t>*> data;
}; // class NeapolitanData


TYPED_TEST_CASE(NeapolitanData, AllCounters);

TYPED_TEST(NeapolitanData, MarginalPValue) {
  auto neapolitan1 = this->data[0]->pValue(0, 1);
  EXPECT_NEAR(neapolitan1, 0.001, 0.0002);

  auto neapolitan2 = this->data[1]->pValue(0, 1);
  EXPECT_NEAR(neapolitan2, 1, 0.0001);

  auto neapolitan3 = this->data[2]->pValue(0, 1);
  EXPECT_NEAR(neapolitan3, 0.46, 0.002);
}


template <typename Counter>
class LizardsData : public testing::Test {
protected:
  void
  SetUp() override {
    auto n = LizardsEnvironment::n;
    auto m = LizardsEnvironment::m;
    auto reader = LizardsEnvironment::reader;
    auto counter = Counter::create(n, m, std::begin(reader->data()));
    data = new DiscreteData<Counter, uint8_t>(counter, reader->varNames());
  }

  void
  TearDown() override {
    delete data;
  }

  DiscreteData<Counter, uint8_t>* data;
}; // class LizardsData

TYPED_TEST_CASE(LizardsData, AllCounters);

TYPED_TEST(LizardsData, MarginalPValue) {
  auto spec_diam  = this->data->pValue(0, 1);
  EXPECT_NEAR(spec_diam, 0.0003845, 0.0001);

  auto spec_hght = this->data->pValue(0, 2);
  EXPECT_NEAR(spec_hght, 0.001257, 0.0001);

  auto diam_hght = this->data->pValue(1, 2);
  EXPECT_NEAR(diam_hght, 0.4357, 0.0001);
}

TYPED_TEST(LizardsData, ConditionalPValue) {
  UintSet<uint8_t> given{2};
  auto spec_diam  = this->data->pValue(0, 1, given);
  EXPECT_NEAR(spec_diam, 0.0009009, 0.0001);

  given = UintSet<uint8_t>{1};
  auto spec_hght = this->data->pValue(0, 2, given);
  EXPECT_NEAR(spec_hght, 0.002708, 0.0001);

  given = UintSet<uint8_t>{0};
  auto diam_hght = this->data->pValue(1, 2, given);
  EXPECT_NEAR(diam_hght, 0.3632, 0.0001);
}


template <typename Counter>
class CoronaryData : public testing::Test {
protected:
  void
  SetUp() override {
    auto n = CoronaryEnvironment::n;
    auto m = CoronaryEnvironment::m;
    auto reader = CoronaryEnvironment::reader;
    auto counter = Counter::create(n, m, std::begin(reader->data()));
    data = new DiscreteData<Counter, uint8_t>(counter, reader->varNames());
  }

  void
  TearDown() override {
    delete data;
  }

  DiscreteData<Counter, uint8_t>* data;
}; // class CoronaryData

TYPED_TEST_CASE(CoronaryData, AllCounters);

TYPED_TEST(CoronaryData, MarginalPValue) {
  auto mwork_pwork  = this->data->pValue(1, 2);
  EXPECT_NEAR(mwork_pwork, 0, 0.0001);

  auto pwork_press = this->data->pValue(2, 3);
  EXPECT_NEAR(pwork_press, 0.7574, 0.0001);

  auto mwork_press = this->data->pValue(1, 3);
  EXPECT_NEAR(mwork_press, 0.0004921, 0.0001);

  auto smoke_press = this->data->pValue(0, 3);
  EXPECT_NEAR(smoke_press, 0.0008954, 0.0001);

  auto press_family = this->data->pValue(3, 5);
  EXPECT_NEAR(press_family, 0.2891, 0.0001);
}

TYPED_TEST(CoronaryData, ConditionalPValue) {
  UintSet<uint8_t> given{0};
  auto press_family  = this->data->pValue(3, 5, given);
  EXPECT_NEAR(press_family, 0.3172, 0.0001);

  given.insert(1);
  press_family = this->data->pValue(3, 5, given);
  EXPECT_NEAR(press_family, 0.2186, 0.0001);

  given.insert(2);
  press_family  = this->data->pValue(3, 5, given);
  EXPECT_NEAR(press_family, 0.2452, 0.0001);

  given.insert(4);
  press_family  = this->data->pValue(3, 5, given);
  EXPECT_NEAR(press_family, 0.0002556, 0.0001);
}

#endif // TEST_DATA_HPP_
