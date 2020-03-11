/**
 * @file DataReader.hpp
 * @brief Implementation of the functions for reading files.
 */
#ifndef DETAIL_DATAREADER_HPP_
#define DETAIL_DATAREADER_HPP_

#include "mxx/collective.hpp"

#include <boost/filesystem.hpp>

#include <algorithm>
#include <fstream>


template <typename DataType>
/**
 * @brief Constructor that allocates space for storing the read data.
 *
 * @param n The number of variables in the file.
 * @param varMajor If the data for the variables is stored contiguously.
 */
DataReader<DataType>::DataReader(
  const uint32_t n,
  const bool varNames,
  const bool varMajor
) : m_data(),
    m_varNames(n),
    m_varMajor(varMajor)
{
  if (!varNames) {
    // Create default variable names [V1, ..., V<numCols>] on every processor
    for (auto i = 0u; i < n; ++i) {
      this->m_varNames[i] = "V" + std::to_string(i + 1);
    }
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
 * @brief Returns if the data for variables is stored contiguously.
 */
bool
DataReader<DataType>::varMajor(
) const
{
  return m_varMajor;
}

template <typename DataType>
/**
 * @brief Store the given name for the variable with the given index.
 *
 * @param idx The index of the variable.
 * @param name Name of the variable with the given index.
 */
void
DataReader<DataType>::varName(
  const uint32_t idx,
  std::string& name
)
{
  // Remove any quotes from the variable names
  name.erase(std::remove(name.begin(), name.end(), '"'), name.end());
  this->m_varNames[idx] = name;
}

template <typename DataType>
/**
 * @brief Extract the data from the given stream.
 *
 * @param is The stream from which data is to be extracted.
 */
DataType
DataReader<DataType>::get(
  const std::string& item
)
{
  std::istringstream is(item);
  if (std::is_same<DataType, uint8_t>::value) {
    uint32_t temp;
    is >> temp;
    return static_cast<uint8_t>(temp);
  }
  else {
    DataType temp;
    is >> temp;
    return temp;
  }
}

template <typename DataType>
/**
 * @brief Read all the data items from the given line.
 *
 * @param ss The stream corresponding to the line.
 * @param sep The character used for delimiting data points.
 * @param buffer Storage for the data points.
 *
 * @return The number of data points read from the line.
 */
uint32_t
DataReader<DataType>::readLine(
  std::stringstream& ss,
  const char sep,
  DataType* const buffer
)
{
  auto idx = 0u;
  std::string item;
  while (std::getline(ss, item, sep)) {
    buffer[idx++] = get(item);
  }
  return idx;
}

template <typename DataType>
/**
 * @brief Gathers data read on each process to get the
 *        complete dataset on every process.
 *
 * @param buffer Storage for data read by this process.
 * @param numRows Total number of rows in the file.
 * @param rowBlocks Block distributed count of rows.
 * @param numCols Total number of columns in the file.
 * @param comm mxx wrapper for the MPI communicator.
 */
void
DataReader<DataType>::gatherData(
  const std::vector<DataType>& buffer,
  const uint32_t numRows,
  const uint32_t rowBlocks,
  const uint32_t numCols,
  const mxx::comm& comm
)
{
  this->m_data.resize(numRows * numCols);
  if (numRows % comm.size() == 0) {
    // All the processors read the same number of rows
    mxx::allgather(&buffer[0], buffer.size(), &this->m_data[0], comm);
  }
  else {
    // Different processors read different number of rows
    // But, we can figure the sizes out using local computations
    std::vector<size_t> recvSizes(comm.size());
    for (int p = 0; p < comm.size(); ++p) {
      auto pOffset = std::min(p * rowBlocks, numRows);
      recvSizes[p] = std::min(rowBlocks, numRows - pOffset) * numCols;
    }
    mxx::allgatherv(&buffer[0], recvSizes[comm.rank()], &this->m_data[0], recvSizes, comm);
  }
}

template <typename DataType>
/**
 * @brief Transpose the given dataset and store it.
 *
 * @param buffer Storage for the original dataset.
 * @param a Major dimension for the transposed dataset.
 * @param b Minor dimension for the transposed dataset.
 */
void
DataReader<DataType>::transposeData(
  const std::vector<DataType>&& buffer,
  const uint32_t a,
  const uint32_t b
)
{
  this->m_data.resize(a * b);
  for (auto x = 0u; x < a; ++x) {
    for (auto y = 0u; y < b; ++y) {
      this->m_data[x * b + y] = buffer[y * a + x];
    }
  }
}

template <typename DataType>
/**
 * @brief Constructor that reads data from the file.
 *
 * @param fileName The name of the file to be read.
 * @param numCols Total number of columns (variables) in the file.
 * @param numRows Total number of rows (observations) in the file.
 * @param sep The character used for delimiting data points.
 * @param varNames If the file contains variable names in the top row.
 * @param obsIndices If the file contains observation indices in the first column.
 * @param varMajor If the data for variables should be stored contiguously.
 * @param parallelRead If the data should be read in parallel.
 */
RowObservationReader<DataType>::RowObservationReader(
  const std::string& fileName,
  const uint32_t numCols,
  const uint32_t numRows,
  const char sep,
  const bool varNames,
  const bool obsIndices,
  const bool varMajor,
  const bool parallelRead
) : DataReader<DataType>(numCols, varNames, varMajor)
{
  mxx::comm comm;
  decltype(this->m_data) buffer;
  auto rowBlocks = numRows;
  auto myRows = numRows;
  auto myOffset = 0u;
  if (parallelRead || comm.is_first()) {
    std::ifstream dataFile(boost::filesystem::canonical(fileName).string());
    std::string line;
    if (varNames) {
      // Consume the variable names before reading the data
      std::getline(dataFile, line);
      std::stringstream ss(line);
      std::string name;
      auto i = 0u;
      while (std::getline(ss, name, sep)) {
        this->varName(i++, name);
      }
    }
    if (parallelRead) {
      rowBlocks = (numRows / comm.size()) + (numRows % comm.size() ? 1 : 0);
      myOffset = std::min(comm.rank() * rowBlocks, numRows);
      for (auto i = 0u; i < myOffset; ++i) {
        dataFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
      }
      myRows = std::min(rowBlocks, numRows - myOffset);
    }

    buffer.resize(myRows * numCols);
    auto t = 0u;
    for (auto j = myOffset; j < myOffset + myRows; ++j) {
      if (!std::getline(dataFile, line)) {
        throw std::runtime_error("Unexpected EOF encountered.");
      }
      std::stringstream ss(line);
      std::string item;
      if (obsIndices) {
        // First get the observation index
        std::getline(ss, item, sep);
      }
      t += this->readLine(ss, sep, &buffer[t]);
    }
    if (parallelRead) {
      t = mxx::reduce(t, 0, comm);
    }
    if (comm.is_first() && (t != (numRows  * numCols))) {
      throw std::runtime_error("Read file did not match the expected dimensions.");
    }
  }
  if (parallelRead) {
    this->gatherData(buffer, numRows, rowBlocks, numCols, comm);
  }
  else {
    if (comm.is_first()) {
      // We are not moving the buffer because we may need it later
      std::swap(this->m_data, buffer);
    }
    else {
      this->m_data.resize(numRows * numCols);
    }
    // Broadcast the read data
    mxx::bcast(&this->m_data[0], numRows * numCols, 0, comm);
    if (varNames) {
      // Can not broadcast a vector of strings
      // Broadcast each string separately
      for (auto& name : this->m_varNames) {
        mxx::bcast(name, 0, comm);
      }
    }
  }
  if (varMajor) {
    // We can transpose now that every processor has all the elements
    std::swap(this->m_data, buffer);
    this->transposeData(std::move(buffer), numCols, numRows);
  }
}

template <typename DataType>
/**
 * @brief Constructor that reads data from the file.
 *
 * @param fileName The name of the file to be read.
 * @param numRows Total number of rows (variables) in the file.
 * @param numCols Total number of columns (observations) in the file.
 * @param sep The character used for delimiting data points.
 * @param varNames If the file contains variable names in the first column.
 * @param obsIndices If the file contains observation indices in the first row.
 * @param varMajor If the data for variables should be stored contiguously.
 * @param parallelRead If the data should be read in parallel.
 */
ColumnObservationReader<DataType>::ColumnObservationReader(
  const std::string& fileName,
  const uint32_t numRows,
  const uint32_t numCols,
  const char sep,
  const bool varNames,
  const bool obsIndices,
  const bool varMajor,
  const bool parallelRead
) : DataReader<DataType>(numRows, varNames, varMajor)
{
  mxx::comm comm;
  decltype(this->m_data) buffer;
  auto rowBlocks = numRows;
  auto myRows = numRows;
  auto myOffset = 0u;
  if (parallelRead || comm.is_first()) {
    std::ifstream dataFile(boost::filesystem::canonical(fileName).string());
    std::string line;
    if (obsIndices) {
      // Discard the observation indices
      std::getline(dataFile, line);
    }
    if (parallelRead) {
      rowBlocks = (numRows / comm.size()) + (numRows % comm.size() ? 1 : 0);
      myOffset = std::min(comm.rank() * rowBlocks, numRows);
      for (auto i = 0u; i < myOffset; ++i) {
        dataFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
      }
      myRows = std::min(rowBlocks, numRows - myOffset);
    }

    buffer.resize(myRows * numCols);
    auto t = 0u;
    for (auto i = myOffset; i < myOffset + myRows; ++i) {
      if (!std::getline(dataFile, line)) {
        throw std::runtime_error("Unexpected EOF encountered.");
      }
      std::stringstream ss(line);
      if (varNames) {
        std::string name;
        std::getline(ss, name, sep);
        this->varName(i, name);
      }
      t += this->readLine(ss, sep, &buffer[t]);
    }
    if (parallelRead) {
      t = mxx::reduce(t, 0, comm);
    }
    if (comm.is_first() && (t != (numRows  * numCols))) {
      throw std::runtime_error("Read file did not match the expected dimensions.");
    }
  }
  if (parallelRead) {
    this->gatherData(buffer, numRows, rowBlocks, numCols, comm);
  }
  else {
    if (comm.is_first()) {
      // We are not moving the buffer because we may need it later
      std::swap(this->m_data, buffer);
    }
    else {
      this->m_data.resize(numRows * numCols);
    }
    // Broadcast the read data from process 0
    mxx::bcast(&this->m_data[0], numRows * numCols, 0, comm);
  }
  if (!varMajor) {
    // We can transpose now that every processor has all the elements
    std::swap(this->m_data, buffer);
    this->transposeData(std::move(buffer), numRows, numCols);
  }
  if (varNames) {
    auto i = 0u;
    for (int p = 0; p < comm.size(); ++p) {
      auto pOffset = std::min(p * rowBlocks, numRows);
      auto pRows = std::min(rowBlocks, numRows - pOffset);
      for (auto r = 0u; r < pRows; ++r, ++i) {
        mxx::bcast(this->m_varNames[i], p, comm);
      }
    }
  }
}

#endif // DETAIL_DATAREADER_HPP_
