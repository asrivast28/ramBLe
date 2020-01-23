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

#include "mxx/comm.hpp"
#include "utils/Logging.hpp"
#include "utils/Timer.hpp"
#include "CTCounter.hpp"

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
  const mxx::comm& comm,
  const Data& data,
  const Var maxConditioning
)
{
  std::stringstream ss;
  if (algoName.compare("gs") == 0) {
    return std::make_unique<GSMB<Data, Var, Set>>(comm, data, maxConditioning);
  }
  ss << "gs,";
  if (algoName.compare("iamb") == 0) {
    return std::make_unique<IAMB<Data, Var, Set>>(comm, data, maxConditioning);
  }
  ss << "iamb,";
  if (algoName.compare("inter.iamb") == 0) {
    return std::make_unique<InterIAMB<Data, Var, Set>>(comm, data, maxConditioning);
  }
  ss << "inter.iamb,";
  if (algoName.compare("mmpc") == 0) {
    return std::make_unique<MMPC<Data, Var, Set>>(comm, data, maxConditioning);
  }
  ss << "mmpc,";
  if (algoName.compare("hiton") == 0) {
    return std::make_unique<HITON<Data, Var, Set>>(comm, data, maxConditioning);
  }
  ss << "hiton,";
  if (algoName.compare("si.hiton.pc") == 0) {
    return std::make_unique<SemiInterleavedHITON<Data, Var, Set>>(comm, data, maxConditioning);
  }
  ss << "si.hiton.pc,";
  if (algoName.compare("getpc") == 0) {
    return std::make_unique<GetPC<Data, Var, Set>>(comm, data, maxConditioning);
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
template <typename Var, typename Size, typename Counter>
std::vector<std::string>
getNeighborhood(
  const Counter& counter,
  const std::vector<std::string>& varNames,
  const ProgramOptions& options
)
{
  mxx::comm comm;
  DiscreteData<Counter, Var> data(counter, varNames, options.alpha());
  Var maxConditioning = static_cast<Var>(std::min(options.numVars(), options.maxConditioning()));
  auto algo = getAlgorithm<Var, UintSet<Var, Size>>(options.algoName(), comm, data, maxConditioning);
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
    if (comm.is_first()) {
      TIMER_ELAPSED("Time taken in getting the neighborhood: ", tNeighborhood);
    }
  }
  if (options.learnNetwork() || !options.outputFile().empty()) {
    TIMER_DECLARE(tNetwork);
    auto g = algo->getNetwork(options.directEdges());
    if (comm.is_first()) {
      TIMER_ELAPSED("Time taken in getting the network: ", tNetwork);
    }
    if ((comm.is_first()) && !options.outputFile().empty()) {
      TIMER_DECLARE(tWrite);
      g.writeGraphviz(options.outputFile(), options.directEdges());
      TIMER_ELAPSED("Time taken in writing the network: ", tWrite);
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
 * @param reader File data reader.
 * @param options Program options provider.
 *
 * @return The list of labels of the variables in the neighborhood.
 */
template <template <typename...> class CounterType, typename FileType>
std::vector<std::string>
getNeighborhood(
  const uint32_t n,
  const uint32_t m,
  std::unique_ptr<FileType>&& reader,
  const ProgramOptions& options
)
{
  std::vector<std::string> varNames(reader->varNames());
  auto counter = CounterType<>::create(n, m, std::begin(reader->data()));
  reader.reset();
  std::vector<std::string> nbrVars;
  if ((n - 1) <= UintSet<uint8_t, std::integral_constant<int, (maxSize<uint8_t>() >> 2)>>::capacity()) {
    nbrVars = getNeighborhood<uint8_t, std::integral_constant<int, (maxSize<uint8_t>() >> 2)>>(counter, varNames, options);
  }
  else if ((n - 1) <= UintSet<uint8_t, std::integral_constant<int, (maxSize<uint8_t>() >> 1)>>::capacity()) {
    nbrVars = getNeighborhood<uint8_t, std::integral_constant<int, (maxSize<uint8_t>() >> 1)>>(counter, varNames, options);
  }
  else if ((n - 1) <= UintSet<uint8_t>::capacity()) {
    nbrVars = getNeighborhood<uint8_t, std::integral_constant<int, maxSize<uint8_t>()>>(counter, varNames, options);
  }
  else if ((n - 1) <= UintSet<uint16_t, std::integral_constant<int, (maxSize<uint16_t>() >> 7)>>::capacity()) {
    nbrVars = getNeighborhood<uint16_t, std::integral_constant<int, (maxSize<uint16_t>() >> 7)>>(counter, varNames, options);
  }
  else if ((n - 1) <= UintSet<uint16_t, std::integral_constant<int, (maxSize<uint16_t>() >> 6)>>::capacity()) {
    nbrVars = getNeighborhood<uint16_t, std::integral_constant<int, (maxSize<uint16_t>() >> 6)>>(counter, varNames, options);
  }
  else if ((n - 1) <= UintSet<uint16_t, std::integral_constant<int, (maxSize<uint16_t>() >> 5)>>::capacity()) {
    nbrVars = getNeighborhood<uint16_t, std::integral_constant<int, (maxSize<uint16_t>() >> 5)>>(counter, varNames, options);
  }
  else if ((n - 1) <= UintSet<uint16_t, std::integral_constant<int, (maxSize<uint16_t>() >> 4)>>::capacity()) {
    nbrVars = getNeighborhood<uint16_t, std::integral_constant<int, (maxSize<uint16_t>() >> 4)>>(counter, varNames, options);
  }
  else if ((n - 1) <= UintSet<uint16_t, std::integral_constant<int, (maxSize<uint16_t>() >> 3)>>::capacity()) {
    nbrVars = getNeighborhood<uint16_t, std::integral_constant<int, (maxSize<uint16_t>() >> 3)>>(counter, varNames, options);
  }
  else if ((n - 1) <= UintSet<uint16_t, std::integral_constant<int, (maxSize<uint16_t>() >> 2)>>::capacity()) {
    nbrVars = getNeighborhood<uint16_t, std::integral_constant<int, (maxSize<uint16_t>() >> 2)>>(counter, varNames, options);
  }
  else if ((n - 1) <= UintSet<uint16_t, std::integral_constant<int, (maxSize<uint16_t>() >> 1)>>::capacity()) {
    nbrVars = getNeighborhood<uint16_t, std::integral_constant<int, (maxSize<uint16_t>() >> 1)>>(counter, varNames, options);
  }
  else if ((n - 1) <= UintSet<uint16_t, std::integral_constant<int, maxSize<uint16_t>()>>::capacity()) {
    nbrVars = getNeighborhood<uint16_t, std::integral_constant<int, maxSize<uint16_t>()>>(counter, varNames, options);
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

  // Set custom error handler (for debugging with working stack-trace on gdb)
  MPI_Errhandler errhandler;
  MPI_Errhandler_create(&my_mpi_errorhandler, &errhandler);
  MPI_Errhandler_set(MPI_COMM_WORLD, errhandler);

  mxx::comm comm;
  if (comm.size() > 1) {
    if (comm.is_first()) {
      std::cerr << "ERROR: Multi-processor execution is not supported yet." << std::endl;
    }
    return 1;
  }

  ProgramOptions options;
  try {
    options.parse(argc, argv);
  }
  catch (const po::error& pe) {
    if (comm.is_first()) {
      std::cerr << pe.what() << std::endl;
    }
    return 1;
  }

  try {
    INIT_LOGGING(options.logLevel());
    uint32_t n = options.numVars();
    uint32_t m = options.numObs();
    TIMER_DECLARE(tRead);
    std::unique_ptr<DataReader<uint8_t>> reader;
    auto varMajor = true;
    if (options.colObs()) {
      reader.reset(new ColumnObservationReader<uint8_t>(options.fileName(), n, m, options.separator(), options.varNames(), options.obsIndices(), varMajor));
    }
    else {
      reader.reset(new RowObservationReader<uint8_t>(options.fileName(), n, m, options.separator(), options.varNames(), options.obsIndices(), varMajor));
    }
    if (comm.is_first()) {
      TIMER_ELAPSED("Time taken in reading the file: ", tRead);
    }

    bool counterFound = false;
    std::stringstream ss;
    std::vector<std::string> nbrVars;
    if (options.counterType().compare("ct") == 0) {
      nbrVars = getNeighborhood<CTCounter>(n, m, std::move(reader), options);
      counterFound = true;
    }
    ss << "ct";
    if (!counterFound) {
      throw std::runtime_error("Requested counter not found. Supported counter types are: {" + ss.str() + "}");
    }

    if (comm.is_first()) {
      for (const auto var: nbrVars) {
        std::cout << var << ",";
      }
      std::cout << std::endl;
    }
  }
  catch (const std::runtime_error& e) {
    std::cerr << "Encountered runtime error during execution:" << std::endl;
    std::cerr << e.what() << std::endl;
    std::cerr << "Aborting." << std::endl;
    return 1;
  }

  // Finalize MPI
  MPI_Finalize();
  return 0;
}
