/**
 * @file ProgramOptions.hpp
 * @brief Declaration of functionality for parsing command line options.
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
#ifndef PROGRAMOPTIONS_HPP_
#define PROGRAMOPTIONS_HPP_

#include <string>

#include <boost/program_options.hpp>
#include <boost/property_tree/ptree.hpp>

namespace po = boost::program_options;
namespace pt = boost::property_tree;

/**
 * Utility class for parsing command line arguments.
 */
class ProgramOptions {
public:
  ProgramOptions();

  void
  parse(int, char**);

  uint32_t
  numVars() const;

  uint32_t
  numObs() const;

  const std::string&
  dataFile() const;

  bool
  parallelRead() const;

  bool
  colObs() const;

  bool
  varNames() const;

  bool
  obsIndices() const;

  char
  separator() const;

  const std::string&
  algoName() const;

  bool
  learnNetwork() const;

  const std::string&
  outputDir() const;

  const std::string&
  counterType() const;

  bool
  forceParallel() const;

  bool
  hostNames() const;

  bool
  warmupMPI() const;

  const std::string&
  logLevel() const;

  const std::string&
  logFile() const;

  const pt::ptree&
  algoConfigs() const;

  ~ProgramOptions();

private:
  po::options_description m_desc;
  pt::ptree m_algoConfigs;
  std::string m_logLevel;
  std::string m_logFile;
  std::string m_dataFile;
  std::string m_algoName;
  std::string m_outputDir;
  std::string m_counterType;
  std::string m_configFile;
  uint32_t m_numVars;
  uint32_t m_numObs;
  char m_separator;
  bool m_parallelRead;
  bool m_colObs;
  bool m_varNames;
  bool m_obsIndices;
  bool m_learnNetwork;
  bool m_directEdges;
  bool m_forceParallel;
  bool m_hostNames;
  bool m_warmupMPI;
}; // class ProgramOptions

#endif // PROGRAMOPTIONS_HPP_
