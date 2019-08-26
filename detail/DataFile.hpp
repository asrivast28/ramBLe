/**
 * @file DataFile.hpp
 * @brief Implementation of the functions for reading files.
 */
#ifndef DETAIL_DATAFILE_HPP_
#define DETAIL_DATAFILE_HPP_

#include <fstream>


template <typename DataType>
/**
 * @brief Constructor that reads the files and creates the object.
 *
 * @param fileName The name of the file to be read.
 * @param numCols The number of columns (variables) in the file.
 * @param numRows The number of rows (observations) in the file.
 * @param sep The character used for delimiting data points.
 * @param header If the file contains a header row.
 * @param columnMajor If the data should be stored in column-major format.
 */
SeparatedFile<DataType>::SeparatedFile(
  const std::string& fileName,
  const uint32_t numCols,
  const uint32_t numRows,
  const char sep,
  const bool header,
  const bool columnMajor
) : m_data(numCols*numRows),
    m_header(numCols)
{
  std::ifstream dataFile(fileName);
  std::string line;
  if (header) {
    // Consume the header line before reading the data
    std::getline(dataFile, line);
    std::stringstream ss(line);
    std::string item;
    auto i = 0u;
    while (std::getline(ss, item, sep)) {
      // Remove any quotes from the header variables
      item.erase(std::remove(item.begin(), item.end(), '"'), item.end());
      m_header[i++] = item;
    }
  }
  else {
    // Create default headers [0, ..., numCols-1]
    for (auto i = 0u; i < numCols; ++i) {
      m_header[i] = std::to_string(i);
    }
  }
  auto j = 0u;
  while (std::getline(dataFile, line)) {
    std::stringstream ss(line);
    std::string item;
    auto i = 0u;
    while (std::getline(ss, item, sep)) {
      std::istringstream is(item);
      if (columnMajor) {
        // Store the data in column major format
        is >> m_data[i*numRows + j];
      }
      else {
        // Store the data in row major format
        is >> m_data[j*numCols + i];
      }
      ++i;
      //std::cout << *data << sep;
    }
    ++j;
    //std::cout << std::endl;
  }
}

template <typename DataType>
/**
 * @brief Returns the data stored in the file.
 */
const std::vector<DataType>&
SeparatedFile<DataType>::data(
) const
{
  return m_data;
}

template <typename DataType>
/**
 * @brief Returns the header data of the file.
 */
const std::vector<std::string>&
SeparatedFile<DataType>::header(
) const
{
  return m_header;
}

#endif // DETAIL_DATAFILE_HPP_
