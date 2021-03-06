/**
 * @file GSNetworks.hpp
 * @brief Definitions of expected networks learned using GS algorithm.
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
#ifndef TEST_GSNETWORKS_HPP_
#define TEST_GSNETWORKS_HPP_

#include "BlanketLearning.hpp"
#include "NetworkData.hpp"


using GSAlgorithm = GS<DiscreteData<Counter, uint8_t>, uint8_t, UintSet<uint8_t>>;

template <>
BayesianNetwork<uint8_t>*
ChildData<GSAlgorithm>::expectedBN(
)
{
  auto* expectedBN = new BayesianNetwork<uint8_t>(this->data->varNames());
  expectedBN->addEdge(this->data->varIndex("V2"), this->data->varIndex("V3"), false);
  expectedBN->addEdge(this->data->varIndex("V2"), this->data->varIndex("V8"), false);
  expectedBN->addEdge(this->data->varIndex("V5"), this->data->varIndex("V3"), false);
  expectedBN->addEdge(this->data->varIndex("V3"), this->data->varIndex("V8"), false);
  expectedBN->addEdge(this->data->varIndex("V3"), this->data->varIndex("V9"), false);
  expectedBN->addEdge(this->data->varIndex("V4"), this->data->varIndex("V5"), true);
  expectedBN->addEdge(this->data->varIndex("V4"), this->data->varIndex("V10"), true);
  expectedBN->addEdge(this->data->varIndex("V5"), this->data->varIndex("V11"), true);
  expectedBN->addEdge(this->data->varIndex("V6"), this->data->varIndex("V13"), true);
  expectedBN->addEdge(this->data->varIndex("V6"), this->data->varIndex("V18"), true);
  expectedBN->addEdge(this->data->varIndex("V7"), this->data->varIndex("V15"), true);
  expectedBN->addEdge(this->data->varIndex("V12"), this->data->varIndex("V16"), true);
  return expectedBN;
}

template <>
BayesianNetwork<uint8_t>*
InsuranceData<GSAlgorithm>::expectedBN(
)
{
  auto* expectedBN = new BayesianNetwork<uint8_t>(this->data->varNames());
  expectedBN->addEdge(this->data->varIndex("V1"), this->data->varIndex("V2"), true);
  expectedBN->addEdge(this->data->varIndex("V2"), this->data->varIndex("V14"), false);
  expectedBN->addEdge(this->data->varIndex("V3"), this->data->varIndex("V5"), false);
  expectedBN->addEdge(this->data->varIndex("V3"), this->data->varIndex("V9"), true);
  expectedBN->addEdge(this->data->varIndex("V4"), this->data->varIndex("V10"), false);
  expectedBN->addEdge(this->data->varIndex("V4"), this->data->varIndex("V14"), false);
  expectedBN->addEdge(this->data->varIndex("V4"), this->data->varIndex("V19"), true);
  expectedBN->addEdge(this->data->varIndex("V5"), this->data->varIndex("V7"), false);
  expectedBN->addEdge(this->data->varIndex("V5"), this->data->varIndex("V12"), false);
  expectedBN->addEdge(this->data->varIndex("V25"), this->data->varIndex("V5"), false);
  expectedBN->addEdge(this->data->varIndex("V6"), this->data->varIndex("V8"), true);
  expectedBN->addEdge(this->data->varIndex("V9"), this->data->varIndex("V7"), false);
  expectedBN->addEdge(this->data->varIndex("V9"), this->data->varIndex("V12"), false);
  expectedBN->addEdge(this->data->varIndex("V13"), this->data->varIndex("V10"), false);
  expectedBN->addEdge(this->data->varIndex("V24"), this->data->varIndex("V25"), true);
  return expectedBN;
}

template <>
BayesianNetwork<uint8_t>*
MildewData<GSAlgorithm>::expectedBN(
)
{
  auto* expectedBN = new BayesianNetwork<uint8_t>(this->data->varNames());
  expectedBN->addEdge(this->data->varIndex("V1"), this->data->varIndex("V2"), true);
  expectedBN->addEdge(this->data->varIndex("V4"), this->data->varIndex("V26"), false);
  expectedBN->addEdge(this->data->varIndex("V5"), this->data->varIndex("V7"), true);
  expectedBN->addEdge(this->data->varIndex("V11"), this->data->varIndex("V28"), true);
  expectedBN->addEdge(this->data->varIndex("V17"), this->data->varIndex("V30"), false);
  expectedBN->addEdge(this->data->varIndex("V32"), this->data->varIndex("V26"), false);
  expectedBN->addEdge(this->data->varIndex("V28"), this->data->varIndex("V33"), true);
  expectedBN->addEdge(this->data->varIndex("V34"), this->data->varIndex("V30"), false);
  return expectedBN;
}

template <>
BayesianNetwork<uint8_t>*
AlarmData<GSAlgorithm>::expectedBN(
)
{
  auto* expectedBN = new BayesianNetwork<uint8_t>(this->data->varNames());
  expectedBN->addEdge(this->data->varIndex("V1"), this->data->varIndex("V6"), true);
  expectedBN->addEdge(this->data->varIndex("V2"), this->data->varIndex("V5"), true);
  expectedBN->addEdge(this->data->varIndex("V3"), this->data->varIndex("V5"), true);
  expectedBN->addEdge(this->data->varIndex("V4"), this->data->varIndex("V5"), true);
  expectedBN->addEdge(this->data->varIndex("V4"), this->data->varIndex("V7"), false);
  expectedBN->addEdge(this->data->varIndex("V6"), this->data->varIndex("V7"), false);
  expectedBN->addEdge(this->data->varIndex("V7"), this->data->varIndex("V36"), false);
  expectedBN->addEdge(this->data->varIndex("V8"), this->data->varIndex("V9"), false);
  expectedBN->addEdge(this->data->varIndex("V35"), this->data->varIndex("V9"), false);
  expectedBN->addEdge(this->data->varIndex("V10"), this->data->varIndex("V11"), false);
  expectedBN->addEdge(this->data->varIndex("V12"), this->data->varIndex("V11"), false);
  expectedBN->addEdge(this->data->varIndex("V12"), this->data->varIndex("V35"), true);
  expectedBN->addEdge(this->data->varIndex("V15"), this->data->varIndex("V14"), false);
  expectedBN->addEdge(this->data->varIndex("V34"), this->data->varIndex("V15"), false);
  expectedBN->addEdge(this->data->varIndex("V37"), this->data->varIndex("V15"), false);
  expectedBN->addEdge(this->data->varIndex("V18"), this->data->varIndex("V16"), false);
  expectedBN->addEdge(this->data->varIndex("V20"), this->data->varIndex("V16"), false);
  expectedBN->addEdge(this->data->varIndex("V25"), this->data->varIndex("V16"), false);
  expectedBN->addEdge(this->data->varIndex("V17"), this->data->varIndex("V26"), false);
  expectedBN->addEdge(this->data->varIndex("V20"), this->data->varIndex("V18"), false);
  expectedBN->addEdge(this->data->varIndex("V24"), this->data->varIndex("V18"), false);
  expectedBN->addEdge(this->data->varIndex("V20"), this->data->varIndex("V21"), true);
  expectedBN->addEdge(this->data->varIndex("V22"), this->data->varIndex("V23"), true);
  expectedBN->addEdge(this->data->varIndex("V27"), this->data->varIndex("V26"), false);
  expectedBN->addEdge(this->data->varIndex("V28"), this->data->varIndex("V29"), true);
  return expectedBN;
}

#endif // TEST_GSNETWORKS_HPP_
