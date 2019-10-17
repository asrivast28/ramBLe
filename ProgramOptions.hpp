/**
 * @file ProgramOptions.hpp
 * @brief Declaration of functionality for parsing command line options.
 */
#ifndef PROGRAMOPTIONS_HPP_
#define PROGRAMOPTIONS_HPP_

#include <string>

#include <boost/program_options.hpp>

namespace po = boost::program_options;

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
  fileName() const;

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

  const std::string&
  targetVar() const;

  bool
  discoverMB() const;

  bool
  learnNetwork() const;

  const std::string&
  outputFile() const;

  bool
  directEdges() const;

  const std::string&
  logLevel() const;

  bool
  wallTime() const;

  ~ProgramOptions();

private:
  po::options_description m_desc;
  std::string m_logLevel;
  std::string m_fileName;
  std::string m_algoName;
  std::string m_targetVar;
  std::string m_outputFile;
  uint32_t m_numVars;
  uint32_t m_numObs;
  char m_separator;
  bool m_colObs;
  bool m_varNames;
  bool m_obsIndices;
  bool m_discoverMB;
  bool m_learnNetwork;
  bool m_directEdges;
  bool m_wallTime;
}; // class ProgramOptions

#endif // PROGRAMOPTIONS_HPP_
