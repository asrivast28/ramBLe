/**
 * @file ProgramOptions.cpp
 * @brief Implementation of functionality for parsing command line options.
 */
#include "ProgramOptions.hpp"

#include <boost/filesystem.hpp>


namespace fs = boost::filesystem;

ProgramOptions::ProgramOptions(
) : m_desc("Constraint-based Causal Discovery"),
    m_fileName(),
    m_algoName(),
    m_targetVar(),
    m_outputFile(),
    m_numVars(),
    m_numRows(),
    m_separator(),
    m_varNames(),
    m_discoverMB(),
    m_directEdges(),
    m_wallTime()
{
  m_desc.add_options()
    ("help,h", "Print this message.")
    ("log,l", po::value<std::string>(&m_logLevel)->default_value("error"), "Level of logging (when logging is enabled).")
    ("file,f", po::value<std::string>(&m_fileName), "Name of the file from which dataset is to be read.")
    ("separator,s", po::value<char>(&m_separator)->default_value(','), "Delimiting character in the file.")
    ("varnames,v", po::bool_switch(&m_varNames)->default_value(false), "Read variable names from the first row of the file.")
    ("nvars,n", po::value<uint32_t>(&m_numVars), "Number of variables in the dataset.")
    ("nrows,m", po::value<uint32_t>(&m_numRows), "Number of rows (observations) in the dataset.")
    ("algorithm,a", po::value<std::string>(&m_algoName)->default_value("gs"), "Name of the algorithm to be used.")
    ("target,t", po::value<std::string>(&m_targetVar), "Name of the target variable.")
    ("blanket,b", po::bool_switch(&m_discoverMB)->default_value(false), "Find MB instead of PC for the target var.")
    ("output,o", po::value<std::string>(&m_outputFile), "Name of the file to which the learned network should be written.")
    ("directed,d", po::bool_switch(&m_directEdges)->default_value(false), "Orient the edges in the learned network.")
    ("walltime,w", po::bool_switch(&m_wallTime)->default_value(false), "Time the top level operations.")
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
  if ((vm.count("target") == 0) && (vm.count("output") == 0)) {
    throw po::error("At least one of --target or --output should be specified.");
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

bool
ProgramOptions::varNames(
) const
{
  return m_varNames;
}

char
ProgramOptions::separator(
) const
{
  return m_separator;
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

bool
ProgramOptions::discoverMB(
) const
{
  return m_discoverMB;
}

const std::string&
ProgramOptions::outputFile(
) const
{
  return m_outputFile;
}

bool
ProgramOptions::directEdges(
) const
{
  return m_directEdges;
}

bool
ProgramOptions::wallTime(
) const
{
  return m_wallTime;
}

ProgramOptions::~ProgramOptions(
)
{
}
