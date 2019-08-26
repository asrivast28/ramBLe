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

#include "utils/Logging.hpp"
#include "BVCounter.hpp"
#include "RadCounter.hpp"

#include <bitset>
#include <iostream>
#include <memory>
#include <vector>


/**
 * @brief Gets a pointer to the object of the required MB discovery algorithm.
 *
 * @tparam VarType Type of the variables (expected to be an integral type).
 * @tparam DataType Type of the object which is used for querying data.
 * @param algoName The name of the algorithm.
 * @param data The object which is used for querying data.
 *
 * @return unique_ptr to the object of the given algorithm.
 *         The unique_ptr points to a nullptr if the algorithm is not found.
 */
template <typename VarType, typename DataType>
std::unique_ptr<MBDiscovery<DataType, VarType>>
getAlgorithm(
  const std::string& algoName,
  const DataType& data
)
{
  std::stringstream ss;
  if (algoName.compare("gs") == 0) {
    return std::make_unique<GSMB<DataType, VarType>>(data);
  }
  ss << "gs,";
  if (algoName.compare("iamb") == 0) {
    return std::make_unique<IAMB<DataType, VarType>>(data);
  }
  ss << "iamb,";
  if (algoName.compare("inter.iamb") == 0) {
    return std::make_unique<InterIAMB<DataType, VarType>>(data);
  }
  ss << "inter.iamb,";
  if (algoName.compare("mmpc") == 0) {
    return std::make_unique<MMPC<DataType, VarType>>(data);
  }
  ss << "mmpc,";
  if (algoName.compare("hiton") == 0) {
    return std::make_unique<HITON<DataType, VarType>>(data);
  }
  ss << "hiton,";
  if (algoName.compare("getpc") == 0) {
    return std::make_unique<GetPC<DataType, VarType>>(data);
  }
  ss << "getpc";
  throw std::runtime_error("Requested algorithm not found. Supported algorithms are: {" + ss.str() + "}");
  return std::unique_ptr<MBDiscovery<DataType, VarType>>();
}

/**
 * @brief Gets the Markov Blanket for the given target variable.
 *
 * @tparam VarType Type of the variables (expected to be an integral type).
 * @tparam CounterType Type of the object that provides counting queries.
 * @param counter Object that executes counting queries.
 * @param varNames Names of all the variables.
 * @param algoName The name of the algorithm to be used for MB discovery.
 * @param targetVar The name of the variable for which the MB is to be discovered.
 *
 * @return The list of names in the MB of the given target variable.
 */
template <typename VarType, typename CounterType>
std::vector<std::string>
getMB(
  const CounterType& counter,
  const std::vector<std::string>& varNames,
  const std::string& algoName,
  const std::string& targetVar
)
{
  std::set<uint32_t> discovered;
  Data<CounterType, VarType> data(counter, varNames);
  auto algo = getAlgorithm<VarType>(algoName, data);
  auto target = data.varIndex(targetVar);
  if (target == varNames.size()) {
    throw std::runtime_error("Target variable not found.");
  }
  return data.varNames(algo->getMB(target));
}

/**
 * @brief Gets the Markov Blanket for the given target variable.
 *
 * @tparam DataType Type of the data observed.
 * @param n The total number of variables.
 * @param m The total number of observations.
 * @param data The observed data.
 * @param varNames The names of all the variables.
 * @param algoName The name of the algorithm to be used for MB discovery.
 * @param targetVar The name of the variable for which the MB is to be discovered.
 *
 * @return The list of names in the MB of the given target variable.
 */
template <typename DataType>
std::vector<std::string>
getMBVars(
  const uint32_t& n,
  const uint64_t& m,
  const std::vector<DataType>& data,
  const std::vector<std::string>& varNames,
  const std::string& algoName,
  const std::string& targetVar
)
{
  std::vector<std::string> mbVars;
  if (n <= static_cast<uint32_t>(std::bitset<8>().set().to_ulong())) {
    auto bvc = create_BVCounter<1>(n, m, std::begin(data));
    mbVars = getMB<uint8_t>(bvc, varNames, algoName, targetVar);
  }
  else if (n <= static_cast<uint32_t>(std::bitset<16>().set().to_ulong())) {
    auto bvc = create_BVCounter<2>(n, m, std::begin(data));
    mbVars = getMB<uint16_t>(bvc, varNames, algoName, targetVar);
  }
  // TODO: Investigate the compiler warning in the following case.
  //else if (n <= static_cast<uint32_t>(std::bitset<32>().set().to_ulong())) {
    //auto bvc = create_BVCounter<4>(n, m, std::begin(data));
    //mbVars = getMB<uint32_t>(bvc, varNames, algoName, targetVar);
  //}
  else {
    throw std::runtime_error("The given number of variables is not supported.");
  }
  return mbVars;
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
    auto mbVars = getMBVars(n, m, dataFile.data(), dataFile.header(), options.algoName(), options.targetVar());
    for (const auto var: mbVars) {
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
