/**
 * @file DataFile.hpp
 * @brief Declaration of the functions for reading files.
 */
#ifndef DATAFILE_HPP_
#define DATAFILE_HPP_

#include <vector>

/**
 * @brief Base for all the classes that read data from a file.
 *
 * @tparam DataType The type of the data to be read.
 */
template <typename DataType>
class DataFile {
public:
  DataFile(const uint32_t, const uint32_t);

  const std::vector<DataType>&
  data() const;

  const std::vector<std::string>&
  varNames() const;

protected:
  std::vector<DataType> m_data;
  std::vector<std::string> m_varNames;
}; // class DataFile

/**
 * @brief Class that reads a file with observations arranged in rows.
 *
 * @tparam DataType The type of the data to be read.
 */
template <typename DataType>
class RowObservationFile : public DataFile<DataType> {
public:
  RowObservationFile(const std::string&, const uint32_t, const uint32_t, const char = '\t', const bool = false, const bool = false, const bool = false);
}; // class RowObservationFile

/**
 * @brief Class that reads a file with observations arranged in columns.
 *
 * @tparam DataType The type of the data to be read.
 */
template <typename DataType>
class ColumnObservationFile : public DataFile<DataType> {
public:
  ColumnObservationFile(const std::string&, const uint32_t, const uint32_t, const char = '\t', const bool = false, const bool = false, const bool = false);
}; // class ColumnObservationFile

#include "detail/DataFile.hpp"

#endif // DATAFILE_HPP_
