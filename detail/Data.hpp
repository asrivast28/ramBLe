/**
 * @file Data.hpp
 * @brief Implementation of the the functions used for querying the data.
 */
#ifndef DETAIL_DATA_HPP_
#define DETAIL_DATA_HPP_

#include "../SetUtils.hpp"
#include "../UintSet.hpp"

#include "utils/Logging.hpp"

#include <boost/math/distributions/chi_squared.hpp>


/**
 * @brief Class that provides an iterator over all combinations of
 *        the states that different variables can take.
 *
 * @tparam DataType Type of the variable states, expected to be an integral type.
 *
 * This class provides a (supposedly) lightweight way of iterating over
 * all the combinations of states of different variables. Variable i
 * is supposed to take all the states between 0 and bounds[i] - 1.
 *
 * It further allows including a single user-provided state for another
 * variable x, whose index is provided by idxX.
 */
template <typename DataType>
class StateIterator {
public:
  StateIterator(const std::vector<DataType>& bounds)
    : m_bounds(bounds),
      m_state(bounds.size(), 0),
      m_valid(true)
  {
  }

  void
  next()
  {
    auto it = 0u;
    m_valid = false;
    while (it < m_bounds.size()) {
      ++m_state[it];
      if (m_state[it] == m_bounds[it]) {
        m_state[it] = 0;
        ++it;
      }
      else {
        m_valid = true;
        break;
      }
    }
  }

  bool
  valid()
  {
    return m_valid;
  }

  const std::vector<DataType>&
  state() const
  {
    return m_state;
  }

private:
  const std::vector<DataType> m_bounds;
  std::vector<DataType> m_state;
  bool m_valid;
}; // class StateIterator

template <typename CounterType, typename VarType>
/**
 * @brief Default constructor.
 */
Data<CounterType, VarType>::Data(
) : m_counter(),
    m_varNames(),
    m_threshold()
{
}

template <typename CounterType, typename VarType>
/**
 * @brief Constructs the object for querying the given dataset.
 *
 * @param counter Object that executes counting queries.
 * @param varNames Names of all the variables.
 * @param threshold Target nominal type I error rate.
 */
Data<CounterType, VarType>::Data(
  const CounterType& counter,
  const std::vector<std::string>& varNames,
  const double threshold
) : m_counter(counter),
    m_varNames(varNames),
    m_threshold(threshold)
{
  LOG_MESSAGE_IF(numVars() != varNames.size(), error, "Number of variables (%d) != Number of variable names", counter.n(), varNames.size());
}

template <typename CounterType, typename VarType>
/**
 * @return Number of variables in the given dataset.
 */
VarType
Data<CounterType, VarType>::numVars(
) const
{
  return static_cast<VarType>(m_counter.n());
}

template <typename CounterType, typename VarType>
/**
 * @return Number of rows (observations) in the given dataset.
 */
uint32_t
Data<CounterType, VarType>::numRows(
) const
{
  return static_cast<uint32_t>(m_counter.m());
}

template <typename CounterType, typename VarType>
/**
 * @brief Returns the name of a variable.
 *
 * @param x The index of the query variable.
 *
 * @return The name of the query variable.
 */
const std::string&
Data<CounterType, VarType>::varName(
  const VarType x
) const
{
  LOG_MESSAGE_IF(x >= m_varNames.size(), error, "Variable index %d out of range.", static_cast<uint32_t>(x));
  return m_varNames[x];
}

template <typename CounterType, typename VarType>
/**
 * @brief Returns the names of all the variables in the given set.
 *
 * @tparam SetType The type of container for the variable indices.
 * @param vars The indices of all the query variable.
 *
 * @return The name of all the query variables.
 */
template <typename SetType>
std::vector<std::string>
Data<CounterType, VarType>::varNames(
  const SetType& vars
) const
{
  std::vector<std::string> names(vars.size());
  auto i = 0u;
  for (const auto var: vars) {
    LOG_MESSAGE_IF(var >= m_varNames.size(), error, "Variable index %d out of range.", static_cast<uint32_t>(var));
    names[i++] = m_varNames[var];
  }
  return names;
}

template <typename CounterType, typename VarType>
/**
 * @brief Returns the index of a variable.
 *
 * @param name The name of the query variable.
 *
 * @return The index of the query variable.
 */
VarType
Data<CounterType, VarType>::varIndex(
  const std::string& name
) const
{
  VarType x = 0u;
  for (const auto& var: m_varNames) {
    if (var.compare(name) == 0) {
      break;
    }
    ++x;
  }
  LOG_MESSAGE_IF(x == numVars(), error, "Variable with name %s not found.", name);
  return x;
}

template <typename CounterType, typename VarType>
/**
 * @brief Returns indices of multiple variables.
 *
 * @tparam SetType The type of container to be used for returning indices.
 * @param names The list of names of the query variables.
 *
 * @return The indices of all the queried variables.
 */
template <typename SetType>
SetType
Data<CounterType, VarType>::varIndices(
  const std::vector<std::string>& names
) const
{
  auto indices = set_init(SetType(), numVars());
  for (const auto& name: names) {
    indices.insert(indices.end(), this->varIndex(name));
  }
  return indices;
}

template <typename CounterType, typename VarType>
/**
 * @brief Computes the p-value for the variables, given the conditioning set,
 *        and the corresponding degree of freedom.
 *
 * @tparam SetType The type of the container used for indices of the given variables.
 * @param x The index of the first variable.
 * @param y The index of the second variable.
 * @param given The indices of the variables to be conditioned on.
 *
 * @return Pair of degree of freedom and the computed G^2 value.
 */
template <typename SetType>
std::pair<uint32_t, double>
Data<CounterType, VarType>::gSquare(
  const VarType x,
  const VarType y,
  const SetType& given
) const
{
  using data_type = typename CounterType::data_type;

  auto r_x = m_counter.r(x);
  auto r_y = m_counter.r(y);

  uint32_t df = (r_x - 1) * (r_y - 1);
  LOG_MESSAGE(trace, "r_x = %d, r_y = %d", r_x, r_y);
  double gSquare = 0.0;

  std::vector<data_type> r(given.size());
  std::vector<int> pa(given.size());
  auto k = 0u;
  for (auto xk = given.begin(); xk != given.end(); ++xk, ++k) {
    pa[k] = *xk;
    r[k] = m_counter.r(*xk);
    df *= r[k];
  }

  for (auto c = StateIterator<data_type>(r); c.valid(); c.next()) {
    auto base = m_counter.common(pa, c.state());
    auto sk = m_counter.count(base);
    if (sk == 0) {
      continue;
    }
    for (data_type a = 0; a < r_x; ++a) {
      auto count_x = m_counter.common(base, static_cast<int>(x), a);
      auto sik = m_counter.count(count_x);
      if (sik == 0) {
        continue;
      }
      for (data_type b = 0; b < r_y; ++b) {
        auto sjk = m_counter.count(m_counter.common(base, static_cast<int>(y), b));
        if (sjk == 0) {
          continue;
        }
        auto sijk = m_counter.count(m_counter.common(count_x, static_cast<int>(y), b));
        if (sijk == 0) {
          continue;
        }
        LOG_MESSAGE(trace, "a = %d, b = %d", static_cast<int>(a), static_cast<int>(b));
        LOG_MESSAGE(trace, "sk = %d, sik = %d, sjk = %d, sijk = %d", sk, sik, sjk, sijk);
        if (sijk * sk == sik * sjk) {
          LOG_MESSAGE(trace, "component = 0.0");
          continue;
        }
        auto component = sijk * (log(sijk) + log(sk) - log(sik) - log(sjk));
        gSquare += component;
        LOG_MESSAGE(trace, "component = %g", component);
      }
    }
  }
  gSquare *= 2.0;
  LOG_MESSAGE(debug, "df = %d, G-square = %g", df, gSquare);
  return std::make_pair(df, gSquare);
}

template <typename CounterType, typename VarType>
/**
 * @brief Computes the p-value for the variables, given the conditioning set.
 *
 * @tparam SetType The type of the container used for indices of the given variables.
 * @param x The index of the first variable.
 * @param y The index of the second variable.
 * @param given The indices of the variables to be conditioned on.
 */
template <typename SetType>
double
Data<CounterType, VarType>::pValue(
  const VarType x,
  const VarType y,
  const SetType& given
) const
{
  auto ret = this->gSquare(x, y, given);
  if (std::fpclassify(ret.second) == FP_ZERO) {
    return 1.0;
  }
  boost::math::chi_squared dist(ret.first);
  double pValue = 1.0 - boost::math::cdf(dist, ret.second);
  LOG_MESSAGE(debug, "p-value = %g", pValue);
  return pValue;
}

template <typename CounterType, typename VarType>
template <typename SetType>
/**
 * @brief Computes the association score for the given variables, given the conditioning set.
 *
 * @tparam SetType The type of the container used for indices of the given variables.
 * @param x The index of the first variable.
 * @param y The index of the second variable.
 * @param given The indices of the variables to be conditioned on.
 *
 * @return A score in the range [0.0, 1.0] which quantifies
 *         the strength of association between the given variables.
 */
double
Data<CounterType, VarType>::assocScore(
  const VarType x,
  const VarType y,
  const SetType& given
) const
{
  return (1.0 - this->pValue(x, y, given));
}

template <typename CounterType, typename VarType>
template <typename SetType>
/**
 * @brief Checks if the variables are independent, given the conditioning set.
 *
 * @tparam SetType The type of the container used for indices of the given variables.
 * @param x The index of the first variable.
 * @param y The index of the second variable.
 * @param given The indices of the variables to be conditioned on.
 */
bool
Data<CounterType, VarType>::isIndependent(
  const VarType x,
  const VarType y,
  const SetType& given
) const
{
  return std::isgreater(this->pValue(x, y, given), m_threshold);
}

template <typename CounterType, typename VarType>
/**
 * @brief Checks for independence, given the association score.
 *
 * @param assocScore The association score, expected to be in the range [0.0, 1.0].
 */
bool
Data<CounterType, VarType>::isIndependent(
  const double assocScore
) const
{
  return std::isgreater(1.0 - assocScore, m_threshold);
}

template <typename CounterType, typename VarType>
/**
 * @brief Finds the minimum strength of association between the given variables,
 *        conditioned on any subset of the given conditioning set.
 *
 * @tparam SetType The type of the container used for indices of the given variables.
 * @param x The index of the first variable.
 * @param y The index of the second variable.
 * @param given The indices of the variables to be conditioned on.
 * @param maxSize The maximum size of the subset to be tested.
 *
 * @return The computed minimum association score.
 */
template <typename SetType>
double
Data<CounterType, VarType>::minAssocScore(
  const VarType x,
  const VarType y,
  const SetType& given,
  const uint32_t maxSize
) const
{
  auto subsetSize = std::min(static_cast<uint32_t>(given.size()), maxSize);
  double minScore = std::numeric_limits<double>::max();
  for (auto i = 0u; (i <= subsetSize) && std::isgreater(minScore, m_threshold); ++i) {
    SubsetIterator<SetType, VarType> sit(given, i);
    do {
      double thisScore = this->assocScore(x, y, sit.get());
      minScore = std::min(thisScore, minScore);
      sit.next();
    } while (sit.valid() && std::isgreater(minScore, m_threshold));
  }
  LOG_MESSAGE(debug, "minAssocScore = %g", minScore);
  return minScore;
}

template <typename CounterType, typename VarType>
/**
 * @brief Finds the minimum strength of association between the given variables,
 *        conditioned on any subset of the given conditioning set.
 *
 * @tparam SetType The type of the container used for indices of the given variables.
 * @param x The index of the first variable.
 * @param y The index of the second variable.
 * @param given The indices of the variables to be conditioned on.
 * @param seed The indices of the variables to be included in every subset.
 * @param maxSize The maximum size of the subset to be tested.
 *
 * @return The computed minimum association score.
 */
template <typename SetType>
double
Data<CounterType, VarType>::minAssocScore(
  const VarType x,
  const VarType y,
  const SetType& given,
  const SetType& seed,
  const uint32_t maxSize
) const
{
  auto subsetSize = std::min(static_cast<uint32_t>(given.size()), maxSize);
  auto minScore = std::numeric_limits<double>::max();
  for (auto i = 0u; (i <= subsetSize) && std::isgreater(minScore, m_threshold); ++i) {
    SubsetIterator<SetType, VarType> sit(given, i);
    do {
      auto subset = sit.get();
      auto condition = SetType(subset.begin(), subset.end());
      // Always include the seed set in the conditioning set
      condition = set_union(condition, seed);
      auto thisScore = this->assocScore(x, y, condition);
      minScore = std::min(thisScore, minScore);
      sit.next();
    } while (sit.valid() && std::isgreater(minScore, m_threshold));
  }
  LOG_MESSAGE(debug, "minAssocScore = %g", minScore);
  return minScore;
}

template <typename CounterType, typename VarType>
/**
 * @brief Finds the subset of the given conditioning set that minimizes
 *        the strength of association between the given variables.
 *
 * @tparam SetType The type of the container used for indices of the given variables.
 * @param x The index of the first variable.
 * @param y The index of the second variable.
 * @param given The indices of the variables to be conditioned on.
 * @param maxSize The maximum size of the subset to be tested.
 *
 * @return A pair with the computed minimum score and the corresponding subset.
 */
template <typename SetType>
std::pair<double, SetType>
Data<CounterType, VarType>::minAssocScoreSubset(
  const VarType x,
  const VarType y,
  const SetType& given,
  const uint32_t maxSize
) const
{
  auto subsetSize = std::min(static_cast<uint32_t>(given.size()), maxSize);
  auto minScore = std::numeric_limits<double>::max();
  auto z = set_init(SetType(), numVars());
  for (auto i = 0u; (i <= subsetSize) && std::isgreater(minScore, m_threshold); ++i) {
    SubsetIterator<SetType, VarType> sit(given, i);
    do {
      auto thisScore = this->assocScore(x, y, sit.get());
      if (std::isless(thisScore, minScore)) {
        minScore = thisScore;
        auto subset = sit.get();
        z = set_init(SetType(subset.begin(), subset.end()), numVars());
      }
      sit.next();
    } while (sit.valid() && std::isgreater(minScore, m_threshold));
  }
  LOG_MESSAGE(debug, "minAssocScore = %g", minScore);
  return std::make_pair(minScore, z);
}

template <typename CounterType, typename VarType>
/**
 * @brief Checks if the given variables are independent, given any
 *        subset of the given conditioning subset.
 *
 * @tparam SetType The type of the container used for indices of the given variables.
 * @param x The index of the first variable.
 * @param y The index of the second variable.
 * @param given The indices of the variables to be conditioned on.
 * @param maxSize The maximum size of the subset to be tested.
 */
template <typename SetType>
bool
Data<CounterType, VarType>::isIndependentAnySubset(
  const VarType x,
  const VarType y,
  const SetType& given,
  const uint32_t maxSize
) const
{
  auto minScore = this->minAssocScore(x, y, given, maxSize);
  return this->isIndependent(minScore);
}

template <typename CounterType, typename VarType>
/**
 * @brief Checks if the given variables are independent, given any
 *        subset of the given conditioning subset.
 *
 * @tparam SetType The type of the container used for indices of the given variables.
 * @param x The index of the first variable.
 * @param y The index of the second variable.
 * @param given The indices of the variables to be conditioned on.
 * @param seed The indices of the variables to be included in every subset.
 * @param maxSize The maximum size of the subset to be tested.
 */
template <typename SetType>
bool
Data<CounterType, VarType>::isIndependentAnySubset(
  const VarType x,
  const VarType y,
  const SetType& given,
  const SetType& seed,
  const uint32_t maxSize
) const
{
  auto minScore = this->minAssocScore(x, y, given, seed, maxSize);
  return this->isIndependent(minScore);
}

#endif // DETAIL_DATA_HPP_
