/**
 * @file PCStableNetworks.hpp
 * @brief Definitions of expected networks learned using PC-stable algorithm.
 * @author Ankit Srivastava <asrivast@gatech.edu>
 *
 * Copyright 2020 Georgia Institute of Technology
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
#ifndef TEST_PCSTABLENETWORKS_HPP_
#define TEST_PCSTABLENETWORKS_HPP_

#include "BlanketLearning.hpp"
#include "NetworkData.hpp"


using PCStableAlgorithm = PCStable<DiscreteData<Counter, uint8_t>, uint8_t, UintSet<uint8_t>>;

static
BayesianNetwork<uint8_t>*
childNetwork(
  const DiscreteData<Counter, uint8_t>* const data
)
{
  auto* expectedBN = new BayesianNetwork<uint8_t>(data->varNames());
  expectedBN->addEdge(data->varIndex("V4"), data->varIndex("V10"), true);
  expectedBN->addEdge(data->varIndex("V7"), data->varIndex("V15"), true);
  expectedBN->addEdge(data->varIndex("V2"), data->varIndex("V8"), false);
  expectedBN->addEdge(data->varIndex("V16"), data->varIndex("V2"), false);
  expectedBN->addEdge(data->varIndex("V17"), data->varIndex("V2"), false);
  expectedBN->addEdge(data->varIndex("V5"), data->varIndex("V11"), false);
  expectedBN->addEdge(data->varIndex("V19"), data->varIndex("V5"), false);
  expectedBN->addEdge(data->varIndex("V6"), data->varIndex("V13"), true);
  expectedBN->addEdge(data->varIndex("V6"), data->varIndex("V20"), false);
  expectedBN->addEdge(data->varIndex("V14"), data->varIndex("V20"), false);
  expectedBN->addEdge(data->varIndex("V3"), data->varIndex("V8"), false);
  expectedBN->addEdge(data->varIndex("V3"), data->varIndex("V9"), true);
  expectedBN->addEdge(data->varIndex("V3"), data->varIndex("V17"), true);
  expectedBN->addEdge(data->varIndex("V4"), data->varIndex("V18"), true);
  expectedBN->addEdge(data->varIndex("V18"), data->varIndex("V5"), false);
  expectedBN->addEdge(data->varIndex("V6"), data->varIndex("V18"), true);
  expectedBN->addEdge(data->varIndex("V12"), data->varIndex("V15"), true);
  expectedBN->addEdge(data->varIndex("V12"), data->varIndex("V16"), true);
  expectedBN->addEdge(data->varIndex("V12"), data->varIndex("V17"), true);
  expectedBN->addEdge(data->varIndex("V12"), data->varIndex("V19"), true);
  return expectedBN;
}

template <>
BayesianNetwork<uint8_t>*
ChildData<PCStableAlgorithm>::expectedBN(
)
{
  return childNetwork(this->data);
}

static
BayesianNetwork<uint8_t>*
insuranceNetwork(
  const DiscreteData<Counter, uint8_t>* const data
)
{
  auto* expectedBN = new BayesianNetwork<uint8_t>(data->varNames());
  expectedBN->addEdge(data->varIndex("V6"), data->varIndex("V8"), true);
  expectedBN->addEdge(data->varIndex("V6"), data->varIndex("V15"), false);
  expectedBN->addEdge(data->varIndex("V7"), data->varIndex("V24"), false);
  expectedBN->addEdge(data->varIndex("V8"), data->varIndex("V10"), false);
  expectedBN->addEdge(data->varIndex("V8"), data->varIndex("V21"), true);
  expectedBN->addEdge(data->varIndex("V13"), data->varIndex("V10"), false);
  expectedBN->addEdge(data->varIndex("V13"), data->varIndex("V27"), false);
  expectedBN->addEdge(data->varIndex("V15"), data->varIndex("V20"), false);
  expectedBN->addEdge(data->varIndex("V21"), data->varIndex("V20"), false);
  expectedBN->addEdge(data->varIndex("V24"), data->varIndex("V23"), false);
  expectedBN->addEdge(data->varIndex("V24"), data->varIndex("V25"), false);
  expectedBN->addEdge(data->varIndex("V5"), data->varIndex("V7"), false);
  expectedBN->addEdge(data->varIndex("V5"), data->varIndex("V12"), false);
  expectedBN->addEdge(data->varIndex("V5"), data->varIndex("V17"), false);
  expectedBN->addEdge(data->varIndex("V5"), data->varIndex("V25"), false);
  expectedBN->addEdge(data->varIndex("V11"), data->varIndex("V17"), false);
  expectedBN->addEdge(data->varIndex("V17"), data->varIndex("V15"), false);
  expectedBN->addEdge(data->varIndex("V1"), data->varIndex("V2"), true);
  expectedBN->addEdge(data->varIndex("V2"), data->varIndex("V14"), true);
  expectedBN->addEdge(data->varIndex("V3"), data->varIndex("V9"), true);
  expectedBN->addEdge(data->varIndex("V3"), data->varIndex("V18"), false);
  expectedBN->addEdge(data->varIndex("V3"), data->varIndex("V19"), false);
  expectedBN->addEdge(data->varIndex("V3"), data->varIndex("V22"), true);
  expectedBN->addEdge(data->varIndex("V9"), data->varIndex("V7"), false);
  expectedBN->addEdge(data->varIndex("V9"), data->varIndex("V12"), false);
  expectedBN->addEdge(data->varIndex("V9"), data->varIndex("V17"), false);
  expectedBN->addEdge(data->varIndex("V9"), data->varIndex("V25"), false);
  expectedBN->addEdge(data->varIndex("V4"), data->varIndex("V10"), false);
  expectedBN->addEdge(data->varIndex("V4"), data->varIndex("V18"), false);
  expectedBN->addEdge(data->varIndex("V4"), data->varIndex("V19"), false);
  expectedBN->addEdge(data->varIndex("V4"), data->varIndex("V27"), false);
  return expectedBN;
}

template <>
BayesianNetwork<uint8_t>*
InsuranceData<PCStableAlgorithm>::expectedBN(
)
{
  return insuranceNetwork(this->data);
}

static
BayesianNetwork<uint8_t>*
mildewNetwork(
  const DiscreteData<Counter, uint8_t>* const data
)
{
  auto* expectedBN = new BayesianNetwork<uint8_t>(data->varNames());
  expectedBN->addEdge(data->varIndex("V5"), data->varIndex("V7"), true);
  expectedBN->addEdge(data->varIndex("V9"), data->varIndex("V10"), true);
  expectedBN->addEdge(data->varIndex("V11"), data->varIndex("V28"), false);
  expectedBN->addEdge(data->varIndex("V17"), data->varIndex("V30"), false);
  expectedBN->addEdge(data->varIndex("V18"), data->varIndex("V24"), true);
  expectedBN->addEdge(data->varIndex("V19"), data->varIndex("V29"), true);
  expectedBN->addEdge(data->varIndex("V25"), data->varIndex("V31"), true);
  expectedBN->addEdge(data->varIndex("V33"), data->varIndex("V28"), false);
  expectedBN->addEdge(data->varIndex("V34"), data->varIndex("V30"), false);
  expectedBN->addEdge(data->varIndex("V4"), data->varIndex("V26"), false);
  expectedBN->addEdge(data->varIndex("V5"), data->varIndex("V12"), true);
  expectedBN->addEdge(data->varIndex("V6"), data->varIndex("V13"), false);
  expectedBN->addEdge(data->varIndex("V12"), data->varIndex("V18"), true);
  expectedBN->addEdge(data->varIndex("V26"), data->varIndex("V13"), false);
  expectedBN->addEdge(data->varIndex("V27"), data->varIndex("V13"), false);
  expectedBN->addEdge(data->varIndex("V32"), data->varIndex("V26"), false);
  return expectedBN;
}

template <>
BayesianNetwork<uint8_t>*
MildewData<PCStableAlgorithm>::expectedBN(
)
{
  return mildewNetwork(this->data);
}

static
BayesianNetwork<uint8_t>*
alarmNetwork(
  const DiscreteData<Counter, uint8_t>* const data
)
{
  auto* expectedBN = new BayesianNetwork<uint8_t>(data->varNames());
  expectedBN->addEdge(data->varIndex("V8"), data->varIndex("V9"), false);
  expectedBN->addEdge(data->varIndex("V22"), data->varIndex("V23"), true);
  expectedBN->addEdge(data->varIndex("V28"), data->varIndex("V29"), true);
  expectedBN->addEdge(data->varIndex("V1"), data->varIndex("V6"), true);
  expectedBN->addEdge(data->varIndex("V4"), data->varIndex("V7"), false);
  expectedBN->addEdge(data->varIndex("V6"), data->varIndex("V7"), false);
  expectedBN->addEdge(data->varIndex("V7"), data->varIndex("V36"), false);
  expectedBN->addEdge(data->varIndex("V11"), data->varIndex("V10"), false);
  expectedBN->addEdge(data->varIndex("V11"), data->varIndex("V12"), false);
  expectedBN->addEdge(data->varIndex("V14"), data->varIndex("V15"), true);
  expectedBN->addEdge(data->varIndex("V15"), data->varIndex("V34"), true);
  expectedBN->addEdge(data->varIndex("V15"), data->varIndex("V37"), false);
  expectedBN->addEdge(data->varIndex("V16"), data->varIndex("V33"), false);
  expectedBN->addEdge(data->varIndex("V17"), data->varIndex("V26"), false);
  expectedBN->addEdge(data->varIndex("V19"), data->varIndex("V20"), false);
  expectedBN->addEdge(data->varIndex("V20"), data->varIndex("V21"), false);
  expectedBN->addEdge(data->varIndex("V24"), data->varIndex("V21"), false);
  expectedBN->addEdge(data->varIndex("V23"), data->varIndex("V24"), true);
  expectedBN->addEdge(data->varIndex("V34"), data->varIndex("V33"), false);
  expectedBN->addEdge(data->varIndex("V36"), data->varIndex("V37"), false);
  expectedBN->addEdge(data->varIndex("V5"), data->varIndex("V2"), false);
  expectedBN->addEdge(data->varIndex("V5"), data->varIndex("V3"), false);
  expectedBN->addEdge(data->varIndex("V4"), data->varIndex("V5"), false);
  expectedBN->addEdge(data->varIndex("V6"), data->varIndex("V5"), false);
  expectedBN->addEdge(data->varIndex("V16"), data->varIndex("V31"), true);
  expectedBN->addEdge(data->varIndex("V18"), data->varIndex("V25"), false);
  expectedBN->addEdge(data->varIndex("V32"), data->varIndex("V20"), false);
  expectedBN->addEdge(data->varIndex("V31"), data->varIndex("V25"), false);
  expectedBN->addEdge(data->varIndex("V25"), data->varIndex("V32"), false);
  expectedBN->addEdge(data->varIndex("V30"), data->varIndex("V26"), false);
  expectedBN->addEdge(data->varIndex("V27"), data->varIndex("V30"), false);
  expectedBN->addEdge(data->varIndex("V29"), data->varIndex("V30"), false);
  expectedBN->addEdge(data->varIndex("V31"), data->varIndex("V32"), false);
  expectedBN->addEdge(data->varIndex("V33"), data->varIndex("V32"), false);
  expectedBN->addEdge(data->varIndex("V35"), data->varIndex("V9"), false);
  expectedBN->addEdge(data->varIndex("V35"), data->varIndex("V10"), false);
  expectedBN->addEdge(data->varIndex("V35"), data->varIndex("V12"), false);
  expectedBN->addEdge(data->varIndex("V35"), data->varIndex("V36"), false);
  return expectedBN;
}

template <>
BayesianNetwork<uint8_t>*
AlarmData<PCStableAlgorithm>::expectedBN(
)
{
  return alarmNetwork(this->data);
}

using PCStable2Algorithm = PCStable2<DiscreteData<Counter, uint8_t>, uint8_t, UintSet<uint8_t>>;

template <>
BayesianNetwork<uint8_t>*
ChildData<PCStable2Algorithm>::expectedBN(
)
{
  return childNetwork(this->data);
}

template <>
BayesianNetwork<uint8_t>*
InsuranceData<PCStable2Algorithm>::expectedBN(
)
{
  return insuranceNetwork(this->data);
}

template <>
BayesianNetwork<uint8_t>*
MildewData<PCStable2Algorithm>::expectedBN(
)
{
  return mildewNetwork(this->data);
}

template <>
BayesianNetwork<uint8_t>*
AlarmData<PCStable2Algorithm>::expectedBN(
)
{
  return alarmNetwork(this->data);
}

#endif // TEST_PCSTABLENETWORKS_HPP_
