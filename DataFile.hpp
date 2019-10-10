/**
 * @file DataFile.hpp
 * @brief Declaration of the functions for reading files.
 */
#ifndef DATAFILE_HPP_
#define DATAFILE_HPP_

#include <vector>


/**
 * @brief Class that provides functionality for reading a separator delimited file.
 *
 * @tparam DataType The type of the data to be read.
 */
template <typename DataType>
class SeparatedFile {
public:
  SeparatedFile(const std::string&, const uint32_t, const uint32_t, const char = '\t', const bool = false, const bool = false);

  const std::vector<DataType>&
  data() const;

  const std::vector<std::string>&
  varNames() const;

private:
  std::vector<DataType> m_data;
  std::vector<std::string> m_varNames;
};

#include "detail/DataFile.hpp"

#endif // DATAFILE_HPP_
