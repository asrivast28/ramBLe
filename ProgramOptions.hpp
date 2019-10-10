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

  const std::string&
  logLevel() const;

  const std::string&
  fileName() const;

  bool
  varNames() const;

  char
  separator() const;

  uint32_t
  numVars() const;

  uint32_t
  numRows() const;

  const std::string&
  algoName() const;

  const std::string&
  targetVar() const;

  bool
  discoverMB() const;

  const std::string&
  outputFile() const;

  bool
  directEdges() const;

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
  uint32_t m_numRows;
  char m_separator;
  bool m_varNames;
  bool m_discoverMB;
  bool m_directEdges;
  bool m_wallTime;
}; // class ProgramOptions

#endif // PROGRAMOPTIONS_HPP_
