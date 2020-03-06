/**
 * @file DataReader.hpp
 * @brief Declaration of the functions for reading files.
 */
#ifndef DATAREADER_HPP_
#define DATAREADER_HPP_

#include <sstream>
#include <vector>

/**
 * @brief Base for all the classes that read data from a file.
 *
 * @tparam DataType The type of the data to be read.
 */
template <typename DataType>
class DataReader {
public:
  DataReader(const uint32_t, const bool, const bool);

  const std::vector<DataType>&
  data() const;

  const std::vector<std::string>&
  varNames() const;

  bool
  varMajor() const;

protected:
  static
  DataType
  get(const std::string&);

  static
  uint32_t
  readLine(std::stringstream&, const char, DataType* const);

protected:
  void
  varName(const uint32_t, std::string&);

  void
  gatherData(const std::vector<DataType>&, const uint32_t, const uint32_t, const uint32_t, const mxx::comm&);

  void
  transposeData(const std::vector<DataType>&& buffer, const uint32_t, const uint32_t);

protected:
  std::vector<DataType> m_data;
  std::vector<std::string> m_varNames;
  const bool m_varMajor;
}; // class DataReader

/**
 * @brief Class that reads a file with observations arranged in rows.
 *
 * @tparam DataType The type of the data to be read.
 */
template <typename DataType>
class RowObservationReader : public DataReader<DataType> {
public:
  RowObservationReader(const std::string&, const uint32_t, const uint32_t, const char = '\t', const bool = false, const bool = false, const bool = false, const bool = false);
}; // class RowObservationReader

/**
 * @brief Class that reads a file with observations arranged in columns.
 *
 * @tparam DataType The type of the data to be read.
 */
template <typename DataType>
class ColumnObservationReader : public DataReader<DataType> {
public:
  ColumnObservationReader(const std::string&, const uint32_t, const uint32_t, const char = '\t', const bool = false, const bool = false, const bool = false, const bool = false);
}; // class ColumnObservationReader

#include "detail/DataReader.hpp"

#endif // DATAREADER_HPP_
