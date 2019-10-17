/**
 * @file DataFile.hpp
 * @brief Implementation of the functions for reading files.
 */
#ifndef DETAIL_DATAFILE_HPP_
#define DETAIL_DATAFILE_HPP_

#include <fstream>


template <typename DataType>
/**
 * @brief Constructor that allocates space for storing the read data.
 *
 * @param n The number of variables in the file.
 * @param m The number of observations in the file.
 */
DataFile<DataType>::DataFile(
  const uint32_t n,
  const uint32_t m
) : m_data(n*m),
    m_varNames(n)
{
}

template <typename DataType>
/**
 * @brief Returns the data stored in the file.
 */
const std::vector<DataType>&
DataFile<DataType>::data(
) const
{
  return m_data;
}

template <typename DataType>
/**
 * @brief Returns the variable names corresponding to the data.
 */
const std::vector<std::string>&
DataFile<DataType>::varNames(
) const
{
  return m_varNames;
}

template <typename DataType>
/**
 * @brief Constructor that reads data from the file.
 *
 * @param fileName The name of the file to be read.
 * @param numCols The number of columns (variables) in the file.
 * @param numRows The number of rows (observations) in the file.
 * @param sep The character used for delimiting data points.
 * @param varNames If the file contains variable names in the top row.
 * @param columnMajor If the data should be stored in column-major format.
 */
RowObservationFile<DataType>::RowObservationFile(
  const std::string& fileName,
  const uint32_t numCols,
  const uint32_t numRows,
  const char sep,
  const bool varNames,
  const bool columnMajor
) : DataFile<DataType>(numCols, numRows)
{
  std::ifstream dataFile(fileName);
  std::string line;
  if (varNames) {
    // Consume the variable names before reading the data
    std::getline(dataFile, line);
    std::stringstream ss(line);
    std::string name;
    auto i = 0u;
    while (std::getline(ss, name, sep)) {
      // Remove any quotes from the variable names
      name.erase(std::remove(name.begin(), name.end(), '"'), name.end());
      this->m_varNames[i++] = name;
    }
  }
  else {
    // Create default variable names [0, ..., numCols-1]
    for (auto i = 0u; i < numCols; ++i) {
      this->m_varNames[i] = std::to_string(i);
    }
  }
  auto t = 0u;
  auto j = 0u;
  while (std::getline(dataFile, line)) {
    std::stringstream ss(line);
    std::string item;
    auto i = 0u;
    while (std::getline(ss, item, sep)) {
      std::istringstream is(item);
      if (columnMajor) {
        // Store the data in column major format
        is >> this->m_data[i*numRows + j];
        ++t;
      }
      else {
        // Store the data in row major format
        is >> this->m_data[j*numCols + i];
        ++t;
      }
      ++i;
      //std::cout << *data << sep;
    }
    ++j;
    //std::cout << std::endl;
  }
  if (t != (numRows  * numCols)) {
    throw std::runtime_error("Read file did not match the expected dimensions.");
  }
}

template <typename DataType>
/**
 * @brief Constructor that reads data from the file.
 *
 * @param fileName The name of the file to be read.
 * @param numRows The number of rows (variables) in the file.
 * @param numCols The number of columns (observations) in the file.
 * @param sep The character used for delimiting data points.
 * @param varNames If the file contains variable names in the first column.
 * @param columnMajor If the data should be stored in column-major format.
 */
ColumnObservationFile<DataType>::ColumnObservationFile(
  const std::string& fileName,
  const uint32_t numRows,
  const uint32_t numCols,
  const char sep,
  const bool varNames,
  const bool columnMajor
) : DataFile<DataType>(numRows, numCols)
{
  std::ifstream dataFile(fileName);
  std::string line;
  auto t = 0u;
  auto i = 0u;
  while (std::getline(dataFile, line)) {
    std::stringstream ss(line);
    std::string item;
    auto j = 0u;
    while (std::getline(ss, item, sep)) {
      std::istringstream is(item);
      if (varNames) {
        std::string name;
        is >> name;
        // Remove any quotes from the variable names
        name.erase(std::remove(name.begin(), name.end(), '"'), name.end());
        this->m_varNames[i] = name;
      }
      if (columnMajor) {
        // Store the data in column major format
        is >> this->m_data[i*numCols + j];
        ++t;
      }
      else {
        // Store the data in row major format
        is >> this->m_data[j*numRows + i];
        ++t;
      }
      ++j;
      //std::cout << *data << sep;
    }
    ++i;
    //std::cout << std::endl;
  }
  if (t != (numRows  * numCols)) {
    throw std::runtime_error("Read file did not match the expected dimensions.");
  }
  if (!varNames) {
    // Create default variable names [0, ..., numCols-1]
    for (auto i = 0u; i < numRows; ++i) {
      this->m_varNames[i] = std::to_string(i);
    }
  }
}

#endif // DETAIL_DATAFILE_HPP_
