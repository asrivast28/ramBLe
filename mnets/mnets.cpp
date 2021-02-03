/**
 * @file mnets.cpp
 * @brief The implementation of the main function for mnets,
 *        and other functions that drive the program execution.
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
#include "RawData.hpp"
#include "Genomica.hpp"
#include "LemonTree.hpp"
#include "ProgramOptions.hpp"
#include "UintSet.hpp"

#include "mxx/env.hpp"
#include "utils/Logging.hpp"
#include "utils/Timer.hpp"

#include <boost/asio/ip/host_name.hpp>

#include <iostream>
#include <memory>
#include <vector>


/**
 * @brief Gets a pointer to the object of the required module network learning algorithm.
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
std::unique_ptr<ModuleNetworkLearning<Data, Var, Set>>
getAlgorithm(
  const std::string& algoName,
  const mxx::comm& comm,
  const Data& data
)
{
  std::stringstream ss;
  if (algoName.compare("lemontree") == 0) {
    return std::make_unique<LemonTree<Data, Var, Set>>(comm, data);
  }
  ss << "lemontree";
  if (algoName.compare("genomica") == 0) {
    return std::make_unique<Genomica<Data, Var, Set>>(comm, data);
  }
  ss << ",genomica";
  throw std::runtime_error("Requested algorithm not found. Supported algorithms are: {" + ss.str() + "}");
  return std::unique_ptr<ModuleNetworkLearning<Data, Var, Set>>();
}

/**
 * @brief Learns the module network with the given parameters
 *        and writes it to the given file.
 *
 * @tparam Var Type of the variables (expected to be an integral type).
 * @param options Program options provider.
 */
template <typename Var, typename Size, typename Data>
void
learnNetwork(
  const ProgramOptions& options,
  const mxx::comm& comm,
  const Data& data
)
{
  auto algo = getAlgorithm<Var, UintSet<Var, Size>>(options.algoName(), comm, data);
  comm.barrier();
  TIMER_DECLARE(tNetwork);
  algo->learnNetwork((comm.size() > 1) || options.forceParallel(), options.algoConfigs(), options.outputDir());
  comm.barrier();
  if (comm.is_first()) {
    TIMER_ELAPSED("Time taken in learning the network: ", tNetwork);
  }
}

/**
 * @brief Learns the module network with the given parameters
 *        and writes it to the given file.
 *
 * @tparam FileType Type of the file to be read.
 * @param options Program options provider.
 * @param reader File data reader.
 */
template <template <typename> class Reader, typename DataType>
void
learnNetwork(
  const ProgramOptions& options,
  const mxx::comm& comm,
  std::unique_ptr<Reader<DataType>>&& reader
)
{
  auto n = options.numVars();
  auto m = options.numObs();
  auto s = std::max(n, m);
  if ((s - 1) <= UintSet<uint8_t, std::integral_constant<int, (maxSize<uint8_t>() >> 2)>>::capacity()) {
    RawData<DataType, uint8_t> data(reader->data(), reader->varNames(), static_cast<uint8_t>(n), static_cast<uint8_t>(m));
    learnNetwork<uint8_t, std::integral_constant<int, (maxSize<uint8_t>() >> 2)>>(options, comm, data);
  }
  else if ((s - 1) <= UintSet<uint8_t, std::integral_constant<int, (maxSize<uint8_t>() >> 1)>>::capacity()) {
    RawData<DataType, uint8_t> data(reader->data(), reader->varNames(), static_cast<uint8_t>(n), static_cast<uint8_t>(m));
    learnNetwork<uint8_t, std::integral_constant<int, (maxSize<uint8_t>() >> 1)>>(options, comm, data);
  }
  else if ((s - 1) <= UintSet<uint8_t>::capacity()) {
    RawData<DataType, uint8_t> data(reader->data(), reader->varNames(), static_cast<uint8_t>(n), static_cast<uint8_t>(m));
    learnNetwork<uint8_t, std::integral_constant<int, maxSize<uint8_t>()>>(options, comm, data);
  }
  else if ((s - 1) <= UintSet<uint16_t, std::integral_constant<int, (maxSize<uint16_t>() >> 7)>>::capacity()) {
    RawData<DataType, uint16_t> data(reader->data(), reader->varNames(), static_cast<uint16_t>(n), static_cast<uint16_t>(m));
    learnNetwork<uint16_t, std::integral_constant<int, (maxSize<uint16_t>() >> 7)>>(options, comm, data);
  }
  else if ((s - 1) <= UintSet<uint16_t, std::integral_constant<int, (maxSize<uint16_t>() >> 6)>>::capacity()) {
    RawData<DataType, uint16_t> data(reader->data(), reader->varNames(), static_cast<uint16_t>(n), static_cast<uint16_t>(m));
    learnNetwork<uint16_t, std::integral_constant<int, (maxSize<uint16_t>() >> 6)>>(options, comm, data);
  }
  else if ((s - 1) <= UintSet<uint16_t, std::integral_constant<int, (maxSize<uint16_t>() >> 5)>>::capacity()) {
    RawData<DataType, uint16_t> data(reader->data(), reader->varNames(), static_cast<uint16_t>(n), static_cast<uint16_t>(m));
    learnNetwork<uint16_t, std::integral_constant<int, (maxSize<uint16_t>() >> 5)>>(options, comm, data);
  }
  else if ((s - 1) <= UintSet<uint16_t, std::integral_constant<int, (maxSize<uint16_t>() >> 4)>>::capacity()) {
    RawData<DataType, uint16_t> data(reader->data(), reader->varNames(), static_cast<uint16_t>(n), static_cast<uint16_t>(m));
    learnNetwork<uint16_t, std::integral_constant<int, (maxSize<uint16_t>() >> 4)>>(options, comm, data);
  }
  else if ((s - 1) <= UintSet<uint16_t, std::integral_constant<int, (maxSize<uint16_t>() >> 3)>>::capacity()) {
    RawData<DataType, uint16_t> data(reader->data(), reader->varNames(), static_cast<uint16_t>(n), static_cast<uint16_t>(m));
    learnNetwork<uint16_t, std::integral_constant<int, (maxSize<uint16_t>() >> 3)>>(options, comm, data);
  }
  else if ((s - 1) <= UintSet<uint16_t, std::integral_constant<int, (maxSize<uint16_t>() >> 2)>>::capacity()) {
    RawData<DataType, uint16_t> data(reader->data(), reader->varNames(), static_cast<uint16_t>(n), static_cast<uint16_t>(m));
    learnNetwork<uint16_t, std::integral_constant<int, (maxSize<uint16_t>() >> 2)>>(options, comm, data);
  }
  else if ((s - 1) <= UintSet<uint16_t, std::integral_constant<int, (maxSize<uint16_t>() >> 1)>>::capacity()) {
    RawData<DataType, uint16_t> data(reader->data(), reader->varNames(), static_cast<uint16_t>(n), static_cast<uint16_t>(m));
    learnNetwork<uint16_t, std::integral_constant<int, (maxSize<uint16_t>() >> 1)>>(options, comm, data);
  }
  else if ((s - 1) <= UintSet<uint16_t, std::integral_constant<int, maxSize<uint16_t>()>>::capacity()) {
    RawData<DataType, uint16_t> data(reader->data(), reader->varNames(), static_cast<uint16_t>(n), static_cast<uint16_t>(m));
    learnNetwork<uint16_t, std::integral_constant<int, maxSize<uint16_t>()>>(options, comm, data);
  }
  else {
    throw std::runtime_error("The given number of variables and observations is not supported.");
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
    std::unique_ptr<DataReader<double>> reader;
    constexpr auto varMajor = true;
    if (options.colObs()) {
      reader.reset(new ColumnObservationReader<double>(options.dataFile(), n, m, options.separator(),
                                                       options.varNames(), options.obsIndices(), varMajor, options.parallelRead()));
    }
    else {
      reader.reset(new RowObservationReader<double>(options.dataFile(), n, m, options.separator(),
                                                    options.varNames(), options.obsIndices(), varMajor, options.parallelRead()));
    }
    comm.barrier();
    if (comm.is_first()) {
      TIMER_ELAPSED("Time taken in reading the file: ", tRead);
    }
    learnNetwork(options, comm, std::move(reader));
  }
  catch (const std::runtime_error& e) {
    std::cerr << "Encountered runtime error during execution:" << std::endl;
    std::cerr << e.what() << std::endl;
    std::cerr << "Aborting." << std::endl;
    return 1;
  }

  return 0;
}
