/**
 * @file test.cpp
 * @brief The main file for running all the unit tests.
 */

#include "SetUtilsTests.hpp"
#include "Environment.hpp"
#include "DataTests.hpp"
#include "BlanketLearningTests.hpp"
#include "DirectLearningTests.hpp"
#include "DirectedNetworkTests.hpp"

#include "mxx/comm.hpp"
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
  MPI_Init(&argc, &argv);
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
  MPI_Finalize();

  return result;
}
