/**
 * @file ProgramOptions.cpp
 * @brief Implementation of functionality for parsing command line options.
 */
#include "ProgramOptions.hpp"

#include <boost/filesystem.hpp>


namespace fs = boost::filesystem;

ProgramOptions::ProgramOptions(
) : m_desc("Discover Markov Blankets"),
    m_fileName(),
    m_algoName(),
    m_targetVar(),
    m_numVars(),
    m_numRows()
{
  m_desc.add_options()
    ("help,h", "Print this message.")
    ("log,l", po::value<std::string>(&m_logLevel)->default_value("error"), "Level of logging (when logging is enabled).")
    ("file,f", po::value<std::string>(&m_fileName), "Name of the file from which dataset is to be read.")
    ("nvars,n", po::value<uint32_t>(&m_numVars), "Number of variables in the dataset.")
    ("nrows,m", po::value<uint32_t>(&m_numRows), "Number of rows (observations) in the dataset.")
    ("algorithm,a", po::value<std::string>(&m_algoName)->default_value("gs"), "Name of the algorithm to be used.")
    ("target,t", po::value<std::string>(&m_targetVar), "Name of the target variable.")
    ;
}

void
ProgramOptions::parse(
  int argc,
  char** argv
)
{
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, m_desc), vm);
  po::notify(vm);

  if ((argc == 1) || (vm.count("help") > 0)) {
    std::stringstream ss;
    ss << m_desc;
    throw po::error(ss.str());
  }
  if (!fs::exists(fs::path(m_fileName))) {
    throw po::error("Couldn't find the data file.");
  }
}

const std::string&
ProgramOptions::logLevel(
) const
{
  return m_logLevel;
}

const std::string&
ProgramOptions::fileName(
) const
{
  return m_fileName;
}

uint32_t
ProgramOptions::numVars(
) const
{
  return m_numVars;
}

uint32_t
ProgramOptions::numRows(
) const
{
  return m_numRows;
}

const std::string&
ProgramOptions::algoName(
) const
{
  return m_algoName;
}

const std::string&
ProgramOptions::targetVar(
) const
{
  return m_targetVar;
}

ProgramOptions::~ProgramOptions(
)
{
}
