/**
 * @file DirectedNetworkTests.hpp
 * @brief Unit tests for the directed network.
 */
#ifndef TEST_DIRECTEDNETWORK_HPP_
#define TEST_DIRECTEDNETWORK_HPP_

#include "CTCounter.hpp"
#include "DiscreteData.hpp"
#include "DirectDiscovery.hpp"
#include "UintSet.hpp"


using Counter = CTCounter<1>;
using Algorithm = GSMB<DiscreteData<Counter, uint8_t>, uint8_t, UintSet<uint8_t>>;

class ChildData : public testing::Test {
protected:
  void
  SetUp() override {
    uint32_t n = 20;
    uint32_t m = 10000;
    ColumnObservationReader<uint8_t> reader("child.csv", n, m, ' ', false, false, true);
    auto counter = Counter::create(n, m, std::begin(reader.data()));
    data = new DiscreteData<Counter, uint8_t>(counter, reader.varNames());
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

TEST_F(ChildData, DirectedNetwork) {
  auto computedBN = this->algo->getNetwork(true);

  auto expectedBN = BayesianNetwork<uint8_t>(this->data->varNames());
  expectedBN.addEdge(static_cast<uint8_t>(1), static_cast<uint8_t>(2));
  expectedBN.addEdge(static_cast<uint8_t>(1), static_cast<uint8_t>(7));
  expectedBN.addEdge(static_cast<uint8_t>(2), static_cast<uint8_t>(7));
  expectedBN.addEdge(static_cast<uint8_t>(2), static_cast<uint8_t>(8));
  expectedBN.addEdge(static_cast<uint8_t>(3), static_cast<uint8_t>(4));
  expectedBN.addEdge(static_cast<uint8_t>(3), static_cast<uint8_t>(9));
  expectedBN.addEdge(static_cast<uint8_t>(4), static_cast<uint8_t>(2));
  expectedBN.addEdge(static_cast<uint8_t>(4), static_cast<uint8_t>(3));
  expectedBN.addEdge(static_cast<uint8_t>(4), static_cast<uint8_t>(10));
  expectedBN.addEdge(static_cast<uint8_t>(5), static_cast<uint8_t>(12));
  expectedBN.addEdge(static_cast<uint8_t>(5), static_cast<uint8_t>(17));
  expectedBN.addEdge(static_cast<uint8_t>(6), static_cast<uint8_t>(14));
  expectedBN.addEdge(static_cast<uint8_t>(9), static_cast<uint8_t>(3));
  expectedBN.addEdge(static_cast<uint8_t>(10), static_cast<uint8_t>(4));
  expectedBN.addEdge(static_cast<uint8_t>(11), static_cast<uint8_t>(15));
  expectedBN.addEdge(static_cast<uint8_t>(12), static_cast<uint8_t>(5));
  expectedBN.addEdge(static_cast<uint8_t>(14), static_cast<uint8_t>(6));
  expectedBN.addEdge(static_cast<uint8_t>(15), static_cast<uint8_t>(11));
  // Check false positives
  EXPECT_NE(expectedBN, computedBN);

  expectedBN.addEdge(static_cast<uint8_t>(17), static_cast<uint8_t>(5));
  EXPECT_EQ(expectedBN, computedBN);
}


class InsuranceData : public testing::Test {
protected:
  void
  SetUp() override {
    uint32_t n = 27;
    uint32_t m = 10000;
    ColumnObservationReader<uint8_t> reader("insurance.csv", n, m, ' ', false, false, true);
    auto counter = Counter::create(n, m, std::begin(reader.data()));
    data = new DiscreteData<Counter, uint8_t>(counter, reader.varNames());
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

TEST_F(InsuranceData, DirectedNetwork) {
  auto computedBN = this->algo->getNetwork(true);

  auto expectedBN = BayesianNetwork<uint8_t>(this->data->varNames());
  expectedBN.addEdge(static_cast<uint8_t>(0), static_cast<uint8_t>(1));
  expectedBN.addEdge(static_cast<uint8_t>(1), static_cast<uint8_t>(0));
  expectedBN.addEdge(static_cast<uint8_t>(1), static_cast<uint8_t>(13));
  expectedBN.addEdge(static_cast<uint8_t>(2), static_cast<uint8_t>(4));
  expectedBN.addEdge(static_cast<uint8_t>(2), static_cast<uint8_t>(8));
  expectedBN.addEdge(static_cast<uint8_t>(3), static_cast<uint8_t>(18));
  expectedBN.addEdge(static_cast<uint8_t>(5), static_cast<uint8_t>(7));
  expectedBN.addEdge(static_cast<uint8_t>(6), static_cast<uint8_t>(4));
  expectedBN.addEdge(static_cast<uint8_t>(6), static_cast<uint8_t>(8));
  expectedBN.addEdge(static_cast<uint8_t>(7), static_cast<uint8_t>(5));
  expectedBN.addEdge(static_cast<uint8_t>(8), static_cast<uint8_t>(2));
  expectedBN.addEdge(static_cast<uint8_t>(8), static_cast<uint8_t>(6));
  expectedBN.addEdge(static_cast<uint8_t>(8), static_cast<uint8_t>(11));
  expectedBN.addEdge(static_cast<uint8_t>(9), static_cast<uint8_t>(3));
  expectedBN.addEdge(static_cast<uint8_t>(9), static_cast<uint8_t>(12));
  expectedBN.addEdge(static_cast<uint8_t>(11), static_cast<uint8_t>(4));
  expectedBN.addEdge(static_cast<uint8_t>(11), static_cast<uint8_t>(8));
  expectedBN.addEdge(static_cast<uint8_t>(12), static_cast<uint8_t>(9));
  expectedBN.addEdge(static_cast<uint8_t>(13), static_cast<uint8_t>(1));
  expectedBN.addEdge(static_cast<uint8_t>(13), static_cast<uint8_t>(3));
  expectedBN.addEdge(static_cast<uint8_t>(23), static_cast<uint8_t>(24));
  expectedBN.addEdge(static_cast<uint8_t>(24), static_cast<uint8_t>(4));
  // Check false positives
  EXPECT_NE(expectedBN, computedBN);

  expectedBN.addEdge(static_cast<uint8_t>(24), static_cast<uint8_t>(23));
  EXPECT_EQ(expectedBN, computedBN);
}


class MildewData : public testing::Test {
protected:
  void
  SetUp() override {
    uint32_t n = 35;
    uint32_t m = 10000;
    ColumnObservationReader<uint8_t> reader("mildew.csv", n, m, ' ', false, false, true);
    auto counter = Counter::create(n, m, std::begin(reader.data()));
    data = new DiscreteData<Counter, uint8_t>(counter, reader.varNames());
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

TEST_F(MildewData, DirectedNetwork) {
  auto computedBN = this->algo->getNetwork(true);

  auto expectedBN = BayesianNetwork<uint8_t>(this->data->varNames());
  expectedBN.addEdge(static_cast<uint8_t>(0), static_cast<uint8_t>(1));
  expectedBN.addEdge(static_cast<uint8_t>(1), static_cast<uint8_t>(0));
  expectedBN.addEdge(static_cast<uint8_t>(3), static_cast<uint8_t>(25));
  expectedBN.addEdge(static_cast<uint8_t>(4), static_cast<uint8_t>(6));
  expectedBN.addEdge(static_cast<uint8_t>(6), static_cast<uint8_t>(4));
  expectedBN.addEdge(static_cast<uint8_t>(10), static_cast<uint8_t>(27));
  expectedBN.addEdge(static_cast<uint8_t>(16), static_cast<uint8_t>(29));
  expectedBN.addEdge(static_cast<uint8_t>(27), static_cast<uint8_t>(10));
  expectedBN.addEdge(static_cast<uint8_t>(27), static_cast<uint8_t>(32));
  expectedBN.addEdge(static_cast<uint8_t>(31), static_cast<uint8_t>(25));
  expectedBN.addEdge(static_cast<uint8_t>(32), static_cast<uint8_t>(27));
  // Check false positives
  EXPECT_NE(expectedBN, computedBN);

  expectedBN.addEdge(static_cast<uint8_t>(33), static_cast<uint8_t>(29));
  EXPECT_EQ(expectedBN, computedBN);
}


class AlarmData : public testing::Test {
protected:
  void
  SetUp() override {
    uint32_t n = 37;
    uint32_t m = 10000;
    ColumnObservationReader<uint8_t> reader("alarm.csv", n, m, ' ', false, false, true);
    auto counter = Counter::create(n, m, std::begin(reader.data()));
    data = new DiscreteData<Counter, uint8_t>(counter, reader.varNames());
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

TEST_F(AlarmData, DirectedNetwork) {
  auto computedBN = this->algo->getNetwork(true);

  auto expectedBN = BayesianNetwork<uint8_t>(this->data->varNames());
  expectedBN.addEdge(static_cast<uint8_t>(0), static_cast<uint8_t>(5));
  expectedBN.addEdge(static_cast<uint8_t>(1), static_cast<uint8_t>(4));
  expectedBN.addEdge(static_cast<uint8_t>(2), static_cast<uint8_t>(4));
  expectedBN.addEdge(static_cast<uint8_t>(4), static_cast<uint8_t>(1));
  expectedBN.addEdge(static_cast<uint8_t>(4), static_cast<uint8_t>(2));
  expectedBN.addEdge(static_cast<uint8_t>(4), static_cast<uint8_t>(3));
  expectedBN.addEdge(static_cast<uint8_t>(5), static_cast<uint8_t>(0));
  expectedBN.addEdge(static_cast<uint8_t>(5), static_cast<uint8_t>(6));
  expectedBN.addEdge(static_cast<uint8_t>(6), static_cast<uint8_t>(3));
  expectedBN.addEdge(static_cast<uint8_t>(6), static_cast<uint8_t>(5));
  expectedBN.addEdge(static_cast<uint8_t>(6), static_cast<uint8_t>(35));
  expectedBN.addEdge(static_cast<uint8_t>(7), static_cast<uint8_t>(8));
  expectedBN.addEdge(static_cast<uint8_t>(9), static_cast<uint8_t>(10));
  expectedBN.addEdge(static_cast<uint8_t>(11), static_cast<uint8_t>(10));
  expectedBN.addEdge(static_cast<uint8_t>(11), static_cast<uint8_t>(34));
  expectedBN.addEdge(static_cast<uint8_t>(14), static_cast<uint8_t>(13));
  expectedBN.addEdge(static_cast<uint8_t>(16), static_cast<uint8_t>(25));
  expectedBN.addEdge(static_cast<uint8_t>(17), static_cast<uint8_t>(15));
  expectedBN.addEdge(static_cast<uint8_t>(19), static_cast<uint8_t>(15));
  expectedBN.addEdge(static_cast<uint8_t>(19), static_cast<uint8_t>(17));
  expectedBN.addEdge(static_cast<uint8_t>(19), static_cast<uint8_t>(20));
  expectedBN.addEdge(static_cast<uint8_t>(20), static_cast<uint8_t>(19));
  expectedBN.addEdge(static_cast<uint8_t>(21), static_cast<uint8_t>(22));
  expectedBN.addEdge(static_cast<uint8_t>(22), static_cast<uint8_t>(21));
  expectedBN.addEdge(static_cast<uint8_t>(23), static_cast<uint8_t>(17));
  expectedBN.addEdge(static_cast<uint8_t>(24), static_cast<uint8_t>(15));
  expectedBN.addEdge(static_cast<uint8_t>(26), static_cast<uint8_t>(25));
  expectedBN.addEdge(static_cast<uint8_t>(27), static_cast<uint8_t>(28));
  expectedBN.addEdge(static_cast<uint8_t>(28), static_cast<uint8_t>(27));
  expectedBN.addEdge(static_cast<uint8_t>(33), static_cast<uint8_t>(14));
  expectedBN.addEdge(static_cast<uint8_t>(34), static_cast<uint8_t>(8));
  expectedBN.addEdge(static_cast<uint8_t>(34), static_cast<uint8_t>(11));
  expectedBN.addEdge(static_cast<uint8_t>(35), static_cast<uint8_t>(6));
  // Check false positives
  EXPECT_NE(expectedBN, computedBN);

  expectedBN.addEdge(static_cast<uint8_t>(36), static_cast<uint8_t>(14));
  EXPECT_EQ(expectedBN, computedBN);
}

#endif // TEST_DIRECTEDNETWORK_HPP_
