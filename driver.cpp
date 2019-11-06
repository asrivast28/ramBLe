/**
 * @file driver.cpp
 * @brief The implementation of the main function, and other
 *        functions that drive the program execution.
 */
#include "DiscreteData.hpp"
#include "DataReader.hpp"
#include "DirectDiscovery.hpp"
#include "ProgramOptions.hpp"
#include "TopologicalDiscovery.hpp"
#include "UintSet.hpp"

#include "utils/Logging.hpp"
#include "utils/Timer.hpp"
#include "BVCounter.hpp"
#include "CTCounter.hpp"
#include "RadCounter.hpp"

#include <mpi.h>

#include <iostream>
#include <memory>
#include <vector>


/**
 * @brief Gets a pointer to the object of the required MB discovery algorithm.
 *
 * @tparam Var Type of the variables (expected to be an integral type).
 * @tparam Set Type of set container.
 * @tparam Data Type of the object which is used for querying data.
 * @param algoName The name of the algorithm.
 * @param data The object which is used for querying data.
 *
 * @return unique_ptr to the object of the given algorithm.
 *         The unique_ptr points to a nullptr if the algorithm is not found.
 */
template <typename Var, typename Set, typename Data>
std::unique_ptr<ConstraintBasedDiscovery<Data, Var, Set>>
getAlgorithm(
  const std::string& algoName,
  const Data& data
)
{
  std::stringstream ss;
  if (algoName.compare("gs") == 0) {
    return std::make_unique<GSMB<Data, Var, Set>>(data);
  }
  ss << "gs,";
  if (algoName.compare("iamb") == 0) {
    return std::make_unique<IAMB<Data, Var, Set>>(data);
  }
  ss << "iamb,";
  if (algoName.compare("inter.iamb") == 0) {
    return std::make_unique<InterIAMB<Data, Var, Set>>(data);
  }
  ss << "inter.iamb,";
  if (algoName.compare("mmpc") == 0) {
    return std::make_unique<MMPC<Data, Var, Set>>(data);
  }
  ss << "mmpc,";
  if (algoName.compare("hiton") == 0) {
    return std::make_unique<HITON<Data, Var, Set>>(data);
  }
  ss << "hiton,";
  if (algoName.compare("si.hiton.pc") == 0) {
    return std::make_unique<SemiInterleavedHITON<Data, Var, Set>>(data);
  }
  ss << "si.hiton.pc,";
  if (algoName.compare("getpc") == 0) {
    return std::make_unique<GetPC<Data, Var, Set>>(data);
  }
  ss << "getpc";
  throw std::runtime_error("Requested algorithm not found. Supported algorithms are: {" + ss.str() + "}");
  return std::unique_ptr<ConstraintBasedDiscovery<Data, Var, Set>>();
}

/**
 * @brief Gets the neighborhood for the given target variable.
 *
 * @tparam Var Type of the variables (expected to be an integral type).
 * @tparam Counter Type of the object that provides counting queries.
 * @param counter Object that executes counting queries.
 * @param varNames Names of all the variables.
 * @param options Program options provider.
 *
 * @return The list of labels of the variables in the neighborhood.
 */
template <typename Var, typename Counter>
std::vector<std::string>
getNeighborhood(
  const Counter& counter,
  const std::vector<std::string>& varNames,
  const ProgramOptions& options
)
{
  DiscreteData<Counter, Var> data(counter, varNames);
  auto algo = getAlgorithm<Var, UintSet<Var>>(options.algoName(), data);
  std::vector<std::string> neighborhoodVars;
  if (!options.targetVar().empty()) {
    TIMER_DECLARE(tNeighborhood);
    auto target = data.varIndex(options.targetVar());
    if (target == varNames.size()) {
      throw std::runtime_error("Target variable not found.");
    }
    if (options.discoverMB()) {
      neighborhoodVars = data.varNames(algo->getMB(target));
    }
    else {
      neighborhoodVars = data.varNames(algo->getPC(target));
    }
    TIMER_ELAPSED("Time taken in getting the neighborhood: ", tNeighborhood);
  }
  if (options.learnNetwork() || !options.outputFile().empty()) {
    TIMER_DECLARE(tNetwork);
    auto g = algo->getNetwork(options.directEdges());
    TIMER_ELAPSED("Time taken in getting the network: ", tNetwork);
    if (!options.outputFile().empty()) {
      TIMER_DECLARE(tWrite);
      g.writeGraphviz(options.outputFile());
      TIMER_ELAPSED("Time taken in writing the network: ", tNetwork);
    }
  }
  return neighborhoodVars;
}

/**
 * @brief Gets the neighborhood for the given target variable.
 *
 * @tparam CounterType Type of the counter to be used.
 * @tparam FileType Type of the file to be read.
 * @param n The total number of variables.
 * @param m The total number of observations.
 * @param dataFile File data provider.
 * @param options Program options provider.
 *
 * @return The list of labels of the variables in the neighborhood.
 */
template <template <int, typename...> class CounterType, typename FileType>
std::vector<std::string>
getNeighborhood(
  const uint32_t n,
  const uint32_t m,
  std::unique_ptr<FileType>&& dataFile,
  const ProgramOptions& options
)
{
  std::vector<std::string> varNames(dataFile->varNames());
  std::vector<std::string> nbrVars;
  if ((n - 1) <= UintSet<uint8_t>::capacity()) {
    constexpr int N = maxN<uint8_t>();
    auto counter = CounterType<N>::create(n, m, std::begin(dataFile->data()));
    dataFile.reset();
    nbrVars = getNeighborhood<uint8_t>(counter, varNames, options);
  }
  else if ((n - 1) <= UintSet<uint16_t>::capacity()) {
    constexpr int N = maxN<uint16_t>();
    auto counter = CounterType<N>::create(n, m, std::begin(dataFile->data()));
    dataFile.reset();
    nbrVars = getNeighborhood<uint16_t>(counter, varNames, options);
  }
  else {
    throw std::runtime_error("The given number of variables is not supported.");
  }
  return nbrVars;
}

void
my_mpi_errorhandler(
  MPI_Comm* comm,
  int* errorcode,
  ...
)
{
    // throw exception, enables gdb stack trace analysis
    throw std::runtime_error("MPI Error");
}

int
main(
  int argc,
  char** argv
)
{
  // Set up MPI
  MPI_Init(&argc, &argv);

  // Get communicator size and my rank
  MPI_Comm comm = MPI_COMM_WORLD;
  int p = 0, rank = -1;
  MPI_Comm_size(comm, &p);
  MPI_Comm_rank(comm, &rank);
  if (p > 1) {
    if (rank == 0) {
      std::cerr << "ERROR: Multi-processor execution is not supported yet." << std::endl;
    }
    return 1;
  }

  // Set custom error handler (for debugging with working stack-trace on gdb)
  MPI_Errhandler errhandler;
  MPI_Errhandler_create(&my_mpi_errorhandler, &errhandler);
  MPI_Errhandler_set(comm, errhandler);

  ProgramOptions options;
  try {
    options.parse(argc, argv);
  }
  catch (const po::error& pe) {
    if (rank == 0) {
      std::cerr << pe.what() << std::endl;
    }
    return 1;
  }

  try {
    INIT_LOGGING(options.logLevel());
    uint32_t n = options.numVars();
    uint32_t m = options.numObs();
    TIMER_DECLARE(tRead);
    std::unique_ptr<DataReader<uint8_t>> dataFile;
    if (options.colObs()) {
      dataFile.reset(new ColumnObservationReader<uint8_t>(options.fileName(), n, m, options.separator(), options.varNames(), options.obsIndices(), true));
    }
    else {
      dataFile.reset(new RowObservationReader<uint8_t>(options.fileName(), n, m, options.separator(), options.varNames(), options.obsIndices(), true));
    }
    TIMER_ELAPSED("Time taken in reading the file: ", tRead);

    bool counterFound = false;
    std::stringstream ss;
    std::vector<std::string> nbrVars;
    if (options.counterType().compare("ct") == 0) {
      nbrVars = getNeighborhood<CTCounter>(n, m, std::move(dataFile), options);
      counterFound = true;
    }
    ss << "ct,";
    if (options.counterType().compare("bv") == 0) {
      nbrVars = getNeighborhood<BVCounter>(n, m, std::move(dataFile), options);
      counterFound = true;
    }
    ss << "bv,";
    if (options.counterType().compare("rad") == 0) {
      nbrVars = getNeighborhood<RadCounter>(n, m, std::move(dataFile), options);
      counterFound = true;
    }
    ss << "rad";
    if (!counterFound) {
      throw std::runtime_error("Requested counter not found. Supported counter types are: {" + ss.str() + "}");
    }

    for (const auto var: nbrVars) {
      std::cout << var << ",";
    }
    std::cout << std::endl;
  }
  catch (const std::runtime_error& e) {
    std::cerr << "Encountered runtime error during execution:" << std::endl;
    std::cerr << e.what() << std::endl;
    std::cerr << "Aborting." << std::endl;
  }

  // Finalize MPI
  MPI_Finalize();
  return 0;
}
