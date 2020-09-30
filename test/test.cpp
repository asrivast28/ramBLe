/**
 * @file test.cpp
 * @brief The main file for running all the unit tests.
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

#include "SetUtilsTests.hpp"
#include "Environment.hpp"
#include "DataTests.hpp"
#include "BlanketLearningTests.hpp"
#include "DirectLearningTests.hpp"
#include "SequentialNetworkTests.hpp"
#include "ParallelNetworkTests.hpp"

#include "mxx/comm.hpp"
#include "mxx/env.hpp"
#include "mxx_eventlistener.hpp"
#include "utils/Logging.hpp"


int
main(
  int argc,
  char** argv
)
{
  testing::InitGoogleTest(&argc, argv);
  std::string logLevel = "error";
  if (argc >= 2) {
    logLevel = argv[1];
  }
  mxx::env e(argc, argv);
  mxx::env::set_exception_on_error();
  mxx::comm comm;
  INIT_LOGGING(logLevel);
  int result = 0;
  if (comm.size() > 1) {
    // set up wrapped test listener
    testing::TestEventListeners& listeners = testing::UnitTest::GetInstance()->listeners();
    testing::TestEventListener* default_listener = listeners.Release(listeners.default_result_printer());
    listeners.Append(new mxx_gtest::MpiTestEventListener(comm.rank(), default_listener));
    testing::GTEST_FLAG(color) = "yes";
    if (testing::GTEST_FLAG(filter) == "*") {
      // This is the default test filter
      // Reset it to only run parallel unit tests
      testing::GTEST_FLAG(filter) = "*ParallelNetwork*";
    }
    result = RUN_ALL_TESTS();
    if (!comm.is_first()) {
      result = 0;
    }
  }
  else {
    testing::AddGlobalTestEnvironment(new NeapolitanEnvironment);
    testing::AddGlobalTestEnvironment(new LizardsEnvironment);
    testing::AddGlobalTestEnvironment(new CoronaryEnvironment);
    testing::AddGlobalTestEnvironment(new AsiaEnvironment);
    if (testing::GTEST_FLAG(filter) == "*") {
      // This is the default test filter
      // Reset it to exclude parallel unit tests
      testing::GTEST_FLAG(filter) = "-*ParallelNetwork*";
    }
    result = RUN_ALL_TESTS();
  }

  return result;
}
