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
 * @tparam CounterType Type of the object that provides counting queries.
 * @param counter Object that executes counting queries.
 * @param varNames Names of all the variables.
 * @param options Program options provider.
 *
 * @return The list of labels of the variables in the neighborhood.
 */
template <typename CounterType, typename VarType>
std::vector<std::string>
getNeighborhood(
  const CounterType& counter,
  const std::vector<std::string>& varNames,
  const ProgramOptions& options
)
{
  Data<CounterType, VarType> data(counter, varNames);
  auto algo = getAlgorithm<VarType, UintSet<VarType>>(options.algoName(), data);
  std::vector<std::string> neighborhoodVars;
  if (!options.targetVar().empty()) {
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
  }
  if (!options.outputFile().empty()) {
    auto g = algo->getNetwork(options.directEdges());
    g.writeGraphviz(options.outputFile());
  }
  return neighborhoodVars;
}

/**
 * @brief Gets the neighborhood for the given target variable.
 *
 * @tparam DataType Type of the data observed.
 * @param n The total number of variables.
 * @param m The total number of observations.
 * @param dataFile File data provider.
 * @param options Program options provider.
 *
 * @return The list of labels of the variables in the neighborhood.
 */
template <typename FileType>
std::vector<std::string>
getNeighborhood(
  const uint32_t& n,
  const uint64_t& m,
  const FileType& dataFile,
  const ProgramOptions& options
)
{
  std::vector<std::string> nbrVars;
  if (n <= UintSet<uint8_t>::capacity()) {
    constexpr int N = UintTypeTrait<uint8_t>::N;
    auto bvc = create_BVCounter<N>(n, m, std::begin(dataFile.data()));
    nbrVars = getNeighborhood<BVCounter<N>, uint8_t>(bvc, dataFile.header(), options);
  }
  else if (n <= UintSet<uint16_t>::capacity()) {
    constexpr int N = UintTypeTrait<uint16_t>::N;
    auto bvc = create_BVCounter<N>(n, m, std::begin(dataFile.data()));
    nbrVars = getNeighborhood<BVCounter<N>, uint16_t>(bvc, dataFile.header(), options);
  }
  // TODO: Investigate the compiler warning in the following case.
  //else if (n < UintSet<4>::capacity()) {
    //auto bvc = create_BVCounter<4>(n, m, std::begin(dataFile.data()));
    //nbrVars = getNeighborhood<BVCounter<4>, 4>(bvc, dataFile.header(), options);
  //}
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
    uint32_t m = options.numRows();
    SeparatedFile<uint8_t> dataFile(options.fileName(), n, m, ',', true, true);
    auto nbrVars = getNeighborhood(n, m, dataFile, options);
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
