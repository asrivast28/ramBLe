/**
 * @file test.cpp
 * @brief The main file for running all the unit tests.
 */

#include "SetUtilsTests.hpp"
#include "Environment.hpp"
#include "DataTests.hpp"
#include "DirectDiscoveryTests.hpp"
#include "TopologicalDiscoveryTests.hpp"
#include "DirectedNetworkTests.hpp"

#include "utils/Logging.hpp"

#include <mpi.h>


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
  // Get communicator size and my rank
  MPI_Comm comm = MPI_COMM_WORLD;
  int p = 0, rank = -1;
  MPI_Comm_size(comm, &p);
  MPI_Comm_rank(comm, &rank);
  if (p > 1) {
    if (rank == 0) {
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
