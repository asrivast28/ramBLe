/**
 * @file test.cpp
 * @brief The main file for running all the unit tests.
 */

#include "Common.hpp"
#include "SetUtilsTests.hpp"
#include "DataTests.hpp"
#include "DirectDiscoveryTests.hpp"
#include "TopologicalDiscoveryTests.hpp"

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
  INIT_LOGGING(logLevel);
  return RUN_ALL_TESTS();
}
