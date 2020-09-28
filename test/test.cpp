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
#include "DirectedNetworkTests.hpp"

#include "mxx/comm.hpp"
#include "mxx/env.hpp"
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
  if (comm.size() > 1) {
    if (comm.is_first()) {
      std::cerr << "ERROR: Multi-processor testing is not supported yet." << std::endl;
    }
    return 1;
  }
  INIT_LOGGING(logLevel);
  testing::AddGlobalTestEnvironment(new NeapolitanEnvironment);
  testing::AddGlobalTestEnvironment(new LizardsEnvironment);
  testing::AddGlobalTestEnvironment(new CoronaryEnvironment);
  testing::AddGlobalTestEnvironment(new AsiaEnvironment);
  auto result = RUN_ALL_TESTS();

  return result;
}
