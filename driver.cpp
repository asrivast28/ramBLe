/**
 * @file driver.cpp
 * @brief The implementation of the main function, and other
 *        functions that drive the program execution.
 */
#include "Data.hpp"
#include "DataFile.hpp"
#include "DirectDiscovery.hpp"
#include "ProgramOptions.hpp"
#include "TopologicalDiscovery.hpp"
#include "UintSet.hpp"

#include "utils/Logging.hpp"
#include "utils/Timer.hpp"
#include "BVCounter.hpp"
#include "RadCounter.hpp"

#include <iostream>
#include <memory>
#include <vector>


/**
 * @brief Gets a pointer to the object of the required MB discovery algorithm.
 *
 * @tparam VarType Type of the variables (expected to be an integral type).
 * @tparam SetType Type of set container.
 * @tparam DataType Type of the object which is used for querying data.
 * @param algoName The name of the algorithm.
 * @param data The object which is used for querying data.
 *
 * @return unique_ptr to the object of the given algorithm.
 *         The unique_ptr points to a nullptr if the algorithm is not found.
 */
template <typename VarType, typename SetType, typename DataType>
std::unique_ptr<ConstraintBasedDiscovery<DataType, VarType, SetType>>
getAlgorithm(
  const std::string& algoName,
  const DataType& data
)
{
  std::stringstream ss;
  if (algoName.compare("gs") == 0) {
    return std::make_unique<GSMB<DataType, VarType, SetType>>(data);
  }
  ss << "gs,";
  if (algoName.compare("iamb") == 0) {
    return std::make_unique<IAMB<DataType, VarType, SetType>>(data);
  }
  ss << "iamb,";
  if (algoName.compare("inter.iamb") == 0) {
    return std::make_unique<InterIAMB<DataType, VarType, SetType>>(data);
  }
  ss << "inter.iamb,";
  if (algoName.compare("mmpc") == 0) {
    return std::make_unique<MMPC<DataType, VarType, SetType>>(data);
  }
  ss << "mmpc,";
  if (algoName.compare("hiton") == 0) {
    return std::make_unique<HITON<DataType, VarType, SetType>>(data);
  }
  ss << "hiton,";
  if (algoName.compare("si.hiton.pc") == 0) {
    return std::make_unique<SemiInterleavedHITON<DataType, VarType, SetType>>(data);
  }
  ss << "si.hiton.pc,";
  if (algoName.compare("getpc") == 0) {
    return std::make_unique<GetPC<DataType, VarType, SetType>>(data);
  }
  ss << "getpc";
  throw std::runtime_error("Requested algorithm not found. Supported algorithms are: {" + ss.str() + "}");
  return std::unique_ptr<ConstraintBasedDiscovery<DataType, VarType, SetType>>();
}

/**
 * @brief Gets the neighborhood for the given target variable.
 *
 * @tparam VarType Type of the variables (expected to be an integral type).
 * @tparam Counter Type of the object that provides counting queries.
 * @param counter Object that executes counting queries.
 * @param varNames Names of all the variables.
 * @param options Program options provider.
 *
 * @return The list of labels of the variables in the neighborhood.
 */
template <typename VarType, typename Counter>
std::vector<std::string>
getNeighborhood(
  const Counter& counter,
  const std::vector<std::string>& varNames,
  const ProgramOptions& options
)
{
  Data<Counter, VarType> data(counter, varNames);
  auto algo = getAlgorithm<VarType, UintSet<VarType>>(options.algoName(), data);
  std::vector<std::string> neighborhoodVars;
  if (!options.targetVar().empty()) {
    Timer tNeighborhood;
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
    if (options.wallTime()) {
      std::cout << "Time taken in getting the neighborhood: " << tNeighborhood.elapsed() << " sec" << std::endl;
    }
  }
  if (options.learnNetwork() || !options.outputFile().empty()) {
    Timer tNetwork;
    auto g = algo->getNetwork(options.directEdges());
    if (options.wallTime()) {
      std::cout << "Time taken in getting the network: " << tNetwork.elapsed() << " sec" << std::endl;
    }
    if (!options.outputFile().empty()) {
      Timer tWrite;
      g.writeGraphviz(options.outputFile());
      if (options.wallTime()) {
        std::cout << "Time taken in writing the network: " << tWrite.elapsed() << " sec" << std::endl;
      }
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
  const uint32_t& n,
  const uint64_t& m,
  std::unique_ptr<FileType>&& dataFile,
  const ProgramOptions& options
)
{
  std::vector<std::string> varNames(dataFile->varNames());
  std::vector<std::string> nbrVars;
  if (n <= UintSet<uint8_t>::capacity()) {
    constexpr int N = UintTypeTrait<uint8_t>::N;
    auto counter = CounterType<N>::create(n, m, std::begin(dataFile->data()));
    dataFile.reset();
    nbrVars = getNeighborhood<uint8_t>(counter, varNames, options);
  }
  else if (n <= UintSet<uint16_t>::capacity()) {
    constexpr int N = UintTypeTrait<uint16_t>::N;
    auto counter = CounterType<N>::create(n, m, std::begin(dataFile->data()));
    dataFile.reset();
    nbrVars = getNeighborhood<uint16_t>(counter, varNames, options);
  }
  else if (n < UintSet<uint32_t>::capacity()) {
    constexpr int N = UintTypeTrait<uint32_t>::N;
    auto counter = CounterType<N>::create(n, m, std::begin(dataFile->data()));
    dataFile.reset();
    nbrVars = getNeighborhood<uint32_t>(counter, varNames, options);
  }
  else {
    throw std::runtime_error("The given number of variables is not supported.");
  }
  return nbrVars;
}

int
main(
  int argc,
  char** argv
)
{
  ProgramOptions options;
  try {
    options.parse(argc, argv);
  }
  catch (const po::error& pe) {
    std::cerr << pe.what() << std::endl;
    return 1;
  }

  try {
    INIT_LOGGING(options.logLevel());
    uint32_t n = options.numVars();
    uint32_t m = options.numObs();
    Timer tRead;
    std::unique_ptr<DataFile<uint8_t>> dataFile;
    if (options.colObs()) {
      dataFile.reset(new ColumnObservationFile<uint8_t>(options.fileName(), n, m, options.separator(), options.varNames(), options.obsIndices(), true));
    }
    else {
      dataFile.reset(new RowObservationFile<uint8_t>(options.fileName(), n, m, options.separator(), options.varNames(), options.obsIndices(), true));
    }
    if (options.wallTime()) {
      std::cout << "Time taken in reading the file: " << tRead.elapsed() << " sec" << std::endl;
    }

    bool counterFound = false;
    std::stringstream ss;
    std::vector<std::string> nbrVars;
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
  return 0;
}
