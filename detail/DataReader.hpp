/**
 * @file DataReader.hpp
 * @brief Implementation of the functions for reading files.
 */
#ifndef DETAIL_DATAREADER_HPP_
#define DETAIL_DATAREADER_HPP_

#include "mxx/collective.hpp"

#include <algorithm>
#include <fstream>


template <typename DataType>
/**
 * @brief Constructor that allocates space for storing the read data.
 *
 * @param n The number of variables in the file.
 * @param m The number of observations in the file.
 */
DataReader<DataType>::DataReader(
  const uint32_t n,
  const uint32_t m
) : m_data(n*m),
    m_varNames(n)
{
}

template <typename DataType>
/**
 * @brief Extract the data from the given stream and store it at the given index.
 *
 * @param index The index at which the data is to be stored.
 * @param is The stream from which data is to be extracted.
 */
void
DataReader<DataType>::data(
  const size_t index,
  std::istringstream& is
)
{
  if (std::is_same<DataType, uint8_t>::value) {
    uint32_t temp;
    is >> temp;
    m_data[index] = static_cast<uint8_t>(temp);
  }
  else {
    is >> m_data[index];
  }
}

template <typename DataType>
/**
 * @brief Returns the data read from the file.
 */
const std::vector<DataType>&
DataReader<DataType>::data(
) const
{
  return m_data;
}

template <typename DataType>
/**
 * @brief Returns the variable names corresponding to the data.
 */
const std::vector<std::string>&
DataReader<DataType>::varNames(
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
 * @param obsIndices If the file contains observation indices in the first column.
 * @param columnMajor If the data should be stored in column-major format.
 */
RowObservationReader<DataType>::RowObservationReader(
  const std::string& fileName,
  const uint32_t numCols,
  const uint32_t numRows,
  const char sep,
  const bool varNames,
  const bool obsIndices,
  const bool columnMajor
) : DataReader<DataType>(numCols, numRows)
{
  mxx::comm comm;
  // Read data in processor 0
  if (comm.is_first()) {
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
      // Create default variable names [V1, ..., V<numCols>]
      for (auto i = 0u; i < numCols; ++i) {
        this->m_varNames[i] = "V" + std::to_string(i + 1);
      }
    }
    auto t = 0u;
    auto j = 0u;
    while (std::getline(dataFile, line)) {
      std::stringstream ss(line);
      std::string item;
      if (obsIndices) {
        // First get the observation index
        std::getline(ss, item, sep);
      }
      auto i = 0u;
      while (std::getline(ss, item, sep)) {
        std::istringstream is(item);
        if (columnMajor) {
          // Store the data in column major format
          this->data(i*numRows + j, is);
          ++t;
        }
        else {
          // Store the data in row major format
          this->data(j*numCols + i, is);
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
  // Broadcast the read data
  mxx::bcast(this->m_data, 0, comm);
  // Can not broadcast a vector of strings
  // Broadcast each string separately
  for (auto& name: this->m_varNames) {
    mxx::bcast(name, 0, comm);
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
 * @param obsIndices If the file contains observation indices in the first row.
 * @param columnMajor If the data should be stored in column-major format.
 */
ColumnObservationReader<DataType>::ColumnObservationReader(
  const std::string& fileName,
  const uint32_t numRows,
  const uint32_t numCols,
  const char sep,
  const bool varNames,
  const bool obsIndices,
  const bool columnMajor
) : DataReader<DataType>(numRows, numCols)
{
  mxx::comm comm;
  // Read data from file in processor 0
  if (comm.is_first()) {
    std::ifstream dataFile(fileName);
    std::string line;
    auto t = 0u;
    auto i = 0u;
    if (obsIndices) {
      std::getline(dataFile, line);
    }
    while (std::getline(dataFile, line)) {
      std::stringstream ss(line);
      if (varNames) {
        std::string name;
        std::getline(ss, name, sep);
        // Remove any quotes from the variable names
        name.erase(std::remove(name.begin(), name.end(), '"'), name.end());
        this->m_varNames[i] = name;
      }
      auto j = 0u;
      std::string item;
      while (std::getline(ss, item, sep)) {
        std::istringstream is(item);
        if (columnMajor) {
          // Store the data in column major format
          this->data(i*numCols + j, is);
          ++t;
        }
        else {
          // Store the data in row major format
          this->data(j*numRows + i, is);
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
      // Create default variable names [V1, ..., V<numCols>]
      for (auto i = 0u; i < numRows; ++i) {
        this->m_varNames[i] = "V" + std::to_string(i + 1);
      }
    }
  }
  // Broadcast the read data
  mxx::bcast(this->m_data, 0, comm);
  // Can not broadcast a vector of strings
  // Broadcast each string separately
  for (auto& name: this->m_varNames) {
    mxx::bcast(name, 0, comm);
  }
}

#endif // DETAIL_DATAREADER_HPP_
