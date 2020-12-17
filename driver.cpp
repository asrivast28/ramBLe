/**
 * @file driver.cpp
 * @brief The implementation of the main function, and other
 *        functions that drive the program execution.
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
#include "DataReader.hpp"
#include "BlanketLearning.hpp"
#include "DirectLearning.hpp"
#include "DiscreteData.hpp"
#include "GlobalLearning.hpp"
#include "ProgramOptions.hpp"
#include "UintSet.hpp"

#include "mxx/comm.hpp"
#include "mxx/env.hpp"
#include "utils/Logging.hpp"
#include "utils/Timer.hpp"
#include "BVCounter.hpp"
#include "CTCounter.hpp"
#include "RadCounter.hpp"

#include <boost/asio/ip/host_name.hpp>

#include <iostream>
#include <memory>
#include <vector>


/**
 * @brief Gets a pointer to the object of the required constraint-based algorithm.
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
std::unique_ptr<ConstraintBasedLearning<Data, Var, Set>>
getAlgorithm(
  const std::string& algoName,
  const mxx::comm& comm,
  const Data& data,
  const Var maxConditioning
)
{
  std::stringstream ss;
  if (algoName.compare("gs") == 0) {
    return std::make_unique<GS<Data, Var, Set>>(comm, data, maxConditioning);
  }
  ss << "gs";
  if (algoName.compare("iamb") == 0) {
    return std::make_unique<IAMB<Data, Var, Set>>(comm, data, maxConditioning);
  }
  ss << ",iamb";
  if (algoName.compare("inter.iamb") == 0) {
    return std::make_unique<InterIAMB<Data, Var, Set>>(comm, data, maxConditioning);
  }
  ss << ",inter.iamb";
  if (algoName.compare("mmpc") == 0) {
    return std::make_unique<MMPC<Data, Var, Set>>(comm, data, maxConditioning);
  }
  ss << ",mmpc";
  if (algoName.compare("hiton") == 0) {
    return std::make_unique<HITON<Data, Var, Set>>(comm, data, maxConditioning);
  }
  ss << ",hiton";
  if (algoName.compare("si.hiton.pc") == 0) {
    return std::make_unique<SemiInterleavedHITON<Data, Var, Set>>(comm, data, maxConditioning);
  }
  ss << ",si.hiton.pc";
  if (algoName.compare("getpc") == 0) {
    return std::make_unique<GetPC<Data, Var, Set>>(comm, data, maxConditioning);
  }
  ss << ",getpc";
  if (algoName.compare("pc.stable") == 0) {
    return std::make_unique<PCStable<Data, Var, Set>>(comm, data, maxConditioning);
  }
  ss << ",pc.stable";
  throw std::runtime_error("Requested algorithm not found. Supported algorithms are: {" + ss.str() + "}");
  return std::unique_ptr<ConstraintBasedLearning<Data, Var, Set>>();
}

/**
 * @brief Learns the BN using the given data counter.
 *
 * @tparam Var Type of the variables (expected to be an integral type).
 * @tparam Counter Type of the object that provides counting queries.
 * @param counter Object that executes counting queries.
 * @param varNames Names of all the variables.
 * @param options Program options provider.
 */
template <typename Var, typename Size, typename Counter>
void
learnNetwork(
  const Counter& counter,
  const std::vector<std::string>& varNames,
  const ProgramOptions& options,
  const mxx::comm& comm
)
{
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
      if (options.directEdges()) {
        for (const auto& vs : algo->findVStructures(target)) {
          std::cout << varNames[std::get<1>(vs)] << " -> " << varNames[std::get<2>(vs)] << " <- " <<
                       varNames[std::get<3>(vs)] << std::endl;
        }
      }
    }
    if (comm.is_first()) {
      for (const auto var : neighborhoodVars) {
        std::cout << var << ",";
      }
      std::cout << std::endl;
      TIMER_ELAPSED("Time taken in getting the neighborhood: ", tNeighborhood);
    }
  }
  if (options.learnNetwork() || !options.outputFile().empty()) {
    comm.barrier();
    TIMER_DECLARE(tNetwork);
    auto g = algo->getNetwork(options.directEdges(), (comm.size() > 1) || options.forceParallel(), options.imbalanceThreshold());
    comm.barrier();
    if (comm.is_first()) {
      TIMER_ELAPSED("Time taken in getting the network: ", tNetwork);
    }
    if ((comm.is_first()) && !options.outputFile().empty()) {
      TIMER_DECLARE(tWrite);
      g.writeGraphviz(options.outputFile());
      TIMER_ELAPSED("Time taken in writing the network: ", tWrite);
    }
  }
}

/**
 * @brief Creates the contingency table counter.
 */
template<template <typename...> class CounterType, typename Size, typename Iter>
typename std::enable_if<
  std::is_same<CounterType<>, CTCounter<>>::value,
  CounterType<>
>::type
createCounter(
  const uint32_t n,
  const uint32_t m,
  Iter it
)
{
  return CounterType<>::create(n, m, it);
}

/**
 * @brief Creates the SABNAtk library counters.
 */
template<template <typename...> class CounterType, typename Size, typename Iter>
typename std::enable_if<
  std::is_same<CounterType<Size>, BVCounter<Size>>::value ||
  std::is_same<CounterType<Size>, RadCounter<Size>>::value,
  CounterType<Size>
>::type
createCounter(
  const uint32_t n,
  const uint32_t m,
  Iter it
)
{
  return CounterType<Size>::create(n, m, it);
}

/**
 * @brief Learns the BN from the data in the given file.
 *
 * @tparam CounterType Type of the counter to be used.
 * @tparam FileType Type of the file to be read.
 * @param n The total number of variables.
 * @param m The total number of observations.
 * @param reader File data reader.
 * @param options Program options provider.
 */
template <template <typename...> class CounterType, typename FileType>
void
learnNetwork(
  const uint32_t n,
  const uint32_t m,
  std::unique_ptr<FileType>&& reader,
  const ProgramOptions& options,
  const mxx::comm& comm
)
{
  std::vector<std::string> varNames(reader->varNames());
  std::vector<std::string> nbrVars;
  if ((n - 1) <= UintSet<uint8_t, std::integral_constant<int, (maxSize<uint8_t>() >> 2)>>::capacity()) {
    auto counter = createCounter<CounterType, std::integral_constant<int, (maxSize<uint8_t>() >> 2)>>(n, m, std::begin(reader->data()));
    learnNetwork<uint8_t, std::integral_constant<int, (maxSize<uint8_t>() >> 2)>>(counter, varNames, options, comm);
  }
  else if ((n - 1) <= UintSet<uint8_t, std::integral_constant<int, (maxSize<uint8_t>() >> 1)>>::capacity()) {
    auto counter = createCounter<CounterType, std::integral_constant<int, (maxSize<uint8_t>() >> 1)>>(n, m, std::begin(reader->data()));
    learnNetwork<uint8_t, std::integral_constant<int, (maxSize<uint8_t>() >> 1)>>(counter, varNames, options, comm);
  }
  else if ((n - 1) <= UintSet<uint8_t>::capacity()) {
    auto counter = createCounter<CounterType, std::integral_constant<int, maxSize<uint8_t>()>>(n, m, std::begin(reader->data()));
    learnNetwork<uint8_t, std::integral_constant<int, maxSize<uint8_t>()>>(counter, varNames, options, comm);
  }
  else if ((n - 1) <= UintSet<uint16_t, std::integral_constant<int, (maxSize<uint16_t>() >> 7)>>::capacity()) {
    auto counter = createCounter<CounterType, std::integral_constant<int, (maxSize<uint16_t>() >> 7)>>(n, m, std::begin(reader->data()));
    learnNetwork<uint16_t, std::integral_constant<int, (maxSize<uint16_t>() >> 7)>>(counter, varNames, options, comm);
  }
  else if ((n - 1) <= UintSet<uint16_t, std::integral_constant<int, (maxSize<uint16_t>() >> 6)>>::capacity()) {
    auto counter = createCounter<CounterType, std::integral_constant<int, (maxSize<uint16_t>() >> 6)>>(n, m, std::begin(reader->data()));
    learnNetwork<uint16_t, std::integral_constant<int, (maxSize<uint16_t>() >> 6)>>(counter, varNames, options, comm);
  }
  else if ((n - 1) <= UintSet<uint16_t, std::integral_constant<int, (maxSize<uint16_t>() >> 5)>>::capacity()) {
    auto counter = createCounter<CounterType, std::integral_constant<int, (maxSize<uint16_t>() >> 5)>>(n, m, std::begin(reader->data()));
    learnNetwork<uint16_t, std::integral_constant<int, (maxSize<uint16_t>() >> 5)>>(counter, varNames, options, comm);
  }
  else if ((n - 1) <= UintSet<uint16_t, std::integral_constant<int, (maxSize<uint16_t>() >> 4)>>::capacity()) {
    auto counter = createCounter<CounterType, std::integral_constant<int, (maxSize<uint16_t>() >> 4)>>(n, m, std::begin(reader->data()));
    learnNetwork<uint16_t, std::integral_constant<int, (maxSize<uint16_t>() >> 4)>>(counter, varNames, options, comm);
  }
  else if ((n - 1) <= UintSet<uint16_t, std::integral_constant<int, (maxSize<uint16_t>() >> 3)>>::capacity()) {
    auto counter = createCounter<CounterType, std::integral_constant<int, (maxSize<uint16_t>() >> 3)>>(n, m, std::begin(reader->data()));
    learnNetwork<uint16_t, std::integral_constant<int, (maxSize<uint16_t>() >> 3)>>(counter, varNames, options, comm);
  }
  else if ((n - 1) <= UintSet<uint16_t, std::integral_constant<int, (maxSize<uint16_t>() >> 2)>>::capacity()) {
    auto counter = createCounter<CounterType, std::integral_constant<int, (maxSize<uint16_t>() >> 2)>>(n, m, std::begin(reader->data()));
    learnNetwork<uint16_t, std::integral_constant<int, (maxSize<uint16_t>() >> 2)>>(counter, varNames, options, comm);
  }
  else if ((n - 1) <= UintSet<uint16_t, std::integral_constant<int, (maxSize<uint16_t>() >> 1)>>::capacity()) {
    auto counter = createCounter<CounterType, std::integral_constant<int, (maxSize<uint16_t>() >> 1)>>(n, m, std::begin(reader->data()));
    learnNetwork<uint16_t, std::integral_constant<int, (maxSize<uint16_t>() >> 1)>>(counter, varNames, options, comm);
  }
  else if ((n - 1) <= UintSet<uint16_t, std::integral_constant<int, maxSize<uint16_t>()>>::capacity()) {
    auto counter = createCounter<CounterType, std::integral_constant<int, maxSize<uint16_t>()>>(n, m, std::begin(reader->data()));
    learnNetwork<uint16_t, std::integral_constant<int, maxSize<uint16_t>()>>(counter, varNames, options, comm);
  }
  else {
    throw std::runtime_error("The given number of variables is not supported.");
  }
}

void
warmupMPI(
  const mxx::comm& comm
)
{
  std::vector<uint8_t> send(comm.size());
  std::vector<uint8_t> recv(comm.size());
  // First, warmup Alltoall of size 1
  mxx::all2all(&send[0], 1, &recv[0], comm);
  // Then, warmup Alltoallv of size 1
  std::vector<size_t> sendSizes(comm.size(), 1);
  std::vector<size_t> sendDispls(comm.size());
  std::iota(sendDispls.begin(), sendDispls.end(), 0);
  std::vector<size_t> recvSizes(comm.size(), 1);
  std::vector<size_t> recvDispls(sendDispls);
  mxx::all2allv(&send[0], sendSizes, sendDispls, &recv[0], recvSizes, recvDispls, comm);
}

int
main(
  int argc,
  char** argv
)
{
  // Set up MPI
  TIMER_DECLARE(tInit);

  mxx::env e(argc, argv);
  mxx::env::set_exception_on_error();
  mxx::comm comm;
  comm.barrier();
  if (comm.is_first()) {
    TIMER_ELAPSED("Time taken in initializing MPI: ", tInit);
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

  if (options.hostNames()) {
    auto name = boost::asio::ip::host_name();
    if (comm.is_first()) {
      std::cout << std::endl << "*** Host names ***" << std::endl;
      std::cout << comm.rank() << ": " << name << std::endl;
    }
    for (int i = 1; i < comm.size(); ++i) {
      if (comm.rank() == i) {
        comm.send(name, 0, i);
      }
      if (comm.is_first()) {
        name = comm.recv<std::string>(i, i);
        std::cout << i << ": " << name << std::endl;
      }
    }
    if (comm.is_first()) {
      std::cout << "******" << std::endl;
    }
  }

  if ((comm.size() > 1) && options.warmupMPI()) {
    comm.barrier();
    TIMER_DECLARE(tWarmup);
    warmupMPI(comm);
    comm.barrier();
    if (comm.is_first()) {
      TIMER_ELAPSED("Time taken in warming up MPI: ", tWarmup);
    }
  }

  try {
    std::string logFile = options.logFile();
    if (!logFile.empty() && (comm.size() > 1)) {
      logFile += ".p" + std::to_string(comm.rank());
    }
    INIT_LOGGING(logFile, comm.rank(), options.logLevel());
    uint32_t n = options.numVars();
    uint32_t m = options.numObs();
    if (static_cast<double>(m) >= std::sqrt(std::numeric_limits<uint32_t>::max())) {
      // Warn the user if the number of observations is too big to be handled by uint32_t
      // We use sqrt here because we never multiply more than two observation counts without handling the consequences
      std::cerr << "WARNING: The given number of observations is possibly too big to be handled by 32-bit unsigned integer" << std::endl;
      std::cerr << "         This may result in silent errors because of overflow" << std::endl;
    }
    TIMER_DECLARE(tRead);
    std::unique_ptr<DataReader<uint8_t>> reader;
    constexpr auto varMajor = true;
    if (options.colObs()) {
      reader.reset(new ColumnObservationReader<uint8_t>(options.dataFile(), n, m, options.separator(),
                                                        options.varNames(), options.obsIndices(), varMajor, options.parallelRead()));
    }
    else {
      reader.reset(new RowObservationReader<uint8_t>(options.dataFile(), n, m, options.separator(),
                                                     options.varNames(), options.obsIndices(), varMajor, options.parallelRead()));
    }
    comm.barrier();
    if (comm.is_first()) {
      TIMER_ELAPSED("Time taken in reading the file: ", tRead);
    }

    bool counterFound = false;
    std::stringstream ss;
    std::vector<std::string> nbrVars;
    if (options.counterType().compare("ct") == 0) {
      learnNetwork<CTCounter>(n, m, std::move(reader), options, comm);
      counterFound = true;
    }
    ss << "ct";
    if (options.counterType().compare("bv") == 0) {
      learnNetwork<BVCounter>(n, m, std::move(reader), options, comm);
      counterFound = true;
    }
    ss << ",bv";
    if (options.counterType().compare("rad") == 0) {
      learnNetwork<RadCounter>(n, m, std::move(reader), options, comm);
      counterFound = true;
    }
    ss << ",rad";
    if (!counterFound) {
      throw std::runtime_error("Requested counter not found. Supported counter types are: {" + ss.str() + "}");
    }

  }
  catch (const std::runtime_error& e) {
    std::cerr << "Encountered runtime error during execution:" << std::endl;
    std::cerr << e.what() << std::endl;
    std::cerr << "Aborting." << std::endl;
    return 1;
  }

  return 0;
}
