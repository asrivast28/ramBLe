/**
 * @file DiscreteData.hpp
 * @brief Implementation of the the functions used for querying the data.
 */
#ifndef DETAIL_DISCRETEDATA_HPP_
#define DETAIL_DISCRETEDATA_HPP_

#include "../SetUtils.hpp"
#include "../UintSet.hpp"

#include "utils/Logging.hpp"

#include <boost/math/distributions/chi_squared.hpp>


template <typename Counter, typename Var>
/**
 * @brief Default constructor.
 */
DiscreteData<Counter, Var>::DiscreteData(
) : m_counter(),
    m_varNames(),
    m_threshold()
{
  TIMER_RESET(m_timer);
}

template <typename Counter, typename Var>
/**
 * @brief Constructs the object for querying the given dataset.
 *
 * @param counter Object that executes counting queries.
 * @param varNames Names of all the variables.
 * @param threshold Target nominal type I error rate.
 */
DiscreteData<Counter, Var>::DiscreteData(
  const Counter& counter,
  const std::vector<std::string>& varNames,
  const double threshold
) : m_counter(counter),
    m_varNames(varNames),
    m_threshold(threshold)
{
  LOG_MESSAGE_IF(numVars() != varNames.size(), error, "Number of variables (%d) != Number of variable names", counter.n(), varNames.size());
  TIMER_RESET(m_timer);
}

template <typename Counter, typename Var>
/**
 * @return Number of variables in the given dataset.
 */
Var
DiscreteData<Counter, Var>::numVars(
) const
{
  return static_cast<Var>(m_counter.n());
}

template <typename Counter, typename Var>
/**
 * @return Number of rows (observations) in the given dataset.
 */
uint32_t
DiscreteData<Counter, Var>::numRows(
) const
{
  return static_cast<uint32_t>(m_counter.m());
}

template <typename Counter, typename Var>
/**
 * @brief Returns the name of a variable.
 *
 * @param x The index of the query variable.
 *
 * @return The name of the query variable.
 */
const std::string&
DiscreteData<Counter, Var>::varName(
  const Var x
) const
{
  LOG_MESSAGE_IF(x >= m_varNames.size(), error, "Variable index %d out of range.", static_cast<uint32_t>(x));
  return m_varNames[x];
}

template <typename Counter, typename Var>
/**
 * @brief Returns the names of all the variables in the given set.
 *
 * @tparam Set The type of container for the variable indices.
 * @param vars The indices of all the query variable.
 *
 * @return The name of all the query variables.
 */
template <typename Set>
std::vector<std::string>
DiscreteData<Counter, Var>::varNames(
  const Set& vars
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

template <typename Counter, typename Var>
/**
 * @brief Returns the names of all the variables.
 */
const std::vector<std::string>&
DiscreteData<Counter, Var>::varNames(
) const
{
  return m_varNames;
}

template <typename Counter, typename Var>
/**
 * @brief Returns the index of a variable.
 *
 * @param name The name of the query variable.
 *
 * @return The index of the query variable.
 */
Var
DiscreteData<Counter, Var>::varIndex(
  const std::string& name
) const
{
  Var x = 0u;
  for (const auto& var: m_varNames) {
    if (var.compare(name) == 0) {
      break;
    }
    ++x;
  }
  LOG_MESSAGE_IF(x == numVars(), error, "Variable with name %s not found.", name);
  return x;
}

template <typename Counter, typename Var>
/**
 * @brief Returns indices of multiple variables.
 *
 * @tparam Set The type of container to be used for returning indices.
 * @param names The list of names of the query variables.
 *
 * @return The indices of all the queried variables.
 */
template <typename Set>
Set
DiscreteData<Counter, Var>::varIndices(
  const std::vector<std::string>& names
) const
{
  auto indices = set_init(Set(), numVars());
  for (const auto& name: names) {
    indices.insert(indices.end(), this->varIndex(name));
  }
  return indices;
}

template <typename Counter, typename Var>
/**
 * @brief Computes the p-value for the variables, given the conditioning set.
 *
 * @tparam Set The type of the container used for indices of the given variables.
 * @param x The index of the first variable.
 * @param y The index of the second variable.
 * @param given The indices of the variables to be conditioned on.
 */
template <typename Set>
double
DiscreteData<Counter, Var>::pValue(
  const Var x,
  const Var y,
  const Set& given
) const
{
  TIMER_START(m_timer);
  auto ret = computeGSquare(m_counter, x, y, given);
  TIMER_PAUSE(m_timer);
  if (std::fpclassify(ret.second) == FP_ZERO) {
    return 1.0;
  }
  boost::math::chi_squared dist(ret.first);
  double pValue = 1.0 - boost::math::cdf(dist, ret.second);
  LOG_MESSAGE(debug, "p-value = %g", pValue);
  return pValue;
}

template <typename Counter, typename Var>
template <typename Set>
/**
 * @brief Computes the association score for the given variables, given the conditioning set.
 *
 * @tparam Set The type of the container used for indices of the given variables.
 * @param x The index of the first variable.
 * @param y The index of the second variable.
 * @param given The indices of the variables to be conditioned on.
 *
 * @return A score in the range [0.0, 1.0] which quantifies
 *         the strength of association between the given variables.
 */
double
DiscreteData<Counter, Var>::assocScore(
  const Var x,
  const Var y,
  const Set& given
) const
{
  return (1.0 - this->pValue(x, y, given));
}

template <typename Counter, typename Var>
template <typename Set>
/**
 * @brief Checks if the variables are independent, given the conditioning set.
 *
 * @tparam Set The type of the container used for indices of the given variables.
 * @param x The index of the first variable.
 * @param y The index of the second variable.
 * @param given The indices of the variables to be conditioned on.
 */
bool
DiscreteData<Counter, Var>::isIndependent(
  const Var x,
  const Var y,
  const Set& given
) const
{
  return std::isgreater(this->pValue(x, y, given), m_threshold);
}

template <typename Counter, typename Var>
/**
 * @brief Checks for independence, given the association score.
 *
 * @param assocScore The association score, expected to be in the range [0.0, 1.0].
 */
bool
DiscreteData<Counter, Var>::isIndependent(
  const double assocScore
) const
{
  return std::isgreater(1.0 - assocScore, m_threshold);
}

template <typename Counter, typename Var>
/**
 * @brief Finds the minimum strength of association between the given variables,
 *        conditioned on any subset of the given conditioning set.
 *
 * @tparam Set The type of the container used for indices of the given variables.
 * @param x The index of the first variable.
 * @param y The index of the second variable.
 * @param given The indices of the variables to be conditioned on.
 * @param maxSize The maximum size of the subset to be tested.
 *
 * @return The computed minimum association score.
 */
template <typename Set>
double
DiscreteData<Counter, Var>::minAssocScore(
  const Var x,
  const Var y,
  const Set& given,
  const Var maxSize
) const
{
  auto subsetSize = std::min(static_cast<Var>(given.size()), maxSize);
  double minScore = std::numeric_limits<double>::max();
  for (auto i = 0u; (i <= subsetSize) && std::isgreater(minScore, m_threshold); ++i) {
    SubsetIterator<Set, Var> sit(given, i);
    do {
      double thisScore = this->assocScore(x, y, sit.get());
      minScore = std::min(thisScore, minScore);
      sit.next();
    } while (sit.valid() && std::isgreater(minScore, m_threshold));
  }
  LOG_MESSAGE(debug, "minAssocScore = %g", minScore);
  return minScore;
}

template <typename Counter, typename Var>
/**
 * @brief Finds the minimum strength of association between the given variables,
 *        conditioned on any subset of the given conditioning set.
 *
 * @tparam Set The type of the container used for indices of the given variables.
 * @param x The index of the first variable.
 * @param y The index of the second variable.
 * @param given The indices of the variables to be conditioned on.
 * @param seed The indices of the variables to be included in every subset.
 * @param maxSize The maximum size of the subset to be tested.
 *
 * @return The computed minimum association score.
 */
template <typename Set>
double
DiscreteData<Counter, Var>::minAssocScore(
  const Var x,
  const Var y,
  const Set& given,
  const Set& seed,
  const Var maxSize
) const
{
  auto subsetSize = std::min(static_cast<Var>(given.size()), maxSize);
  auto minScore = std::numeric_limits<double>::max();
  for (auto i = 0u; (i <= subsetSize) && std::isgreater(minScore, m_threshold); ++i) {
    SubsetIterator<Set, Var> sit(given, i);
    do {
      auto subset = sit.get();
      auto condition = Set(subset.begin(), subset.end());
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

template <typename Counter, typename Var>
/**
 * @brief Finds the subset of the given conditioning set that minimizes
 *        the strength of association between the given variables.
 *
 * @tparam Set The type of the container used for indices of the given variables.
 * @param x The index of the first variable.
 * @param y The index of the second variable.
 * @param given The indices of the variables to be conditioned on.
 * @param maxSize The maximum size of the subset to be tested.
 *
 * @return A pair with the computed minimum score and the corresponding subset.
 */
template <typename Set>
std::pair<double, Set>
DiscreteData<Counter, Var>::minAssocScoreSubset(
  const Var x,
  const Var y,
  const Set& given,
  const Var maxSize
) const
{
  auto subsetSize = std::min(static_cast<Var>(given.size()), maxSize);
  auto minScore = std::numeric_limits<double>::max();
  auto z = set_init(Set(), numVars());
  for (auto i = 0u; (i <= subsetSize) && std::isgreater(minScore, m_threshold); ++i) {
    SubsetIterator<Set, Var> sit(given, i);
    do {
      auto thisScore = this->assocScore(x, y, sit.get());
      if (std::isless(thisScore, minScore)) {
        minScore = thisScore;
        auto subset = sit.get();
        z = set_init(Set(subset.begin(), subset.end()), numVars());
      }
      sit.next();
    } while (sit.valid() && std::isgreater(minScore, m_threshold));
  }
  LOG_MESSAGE(debug, "minAssocScore = %g", minScore);
  return std::make_pair(minScore, z);
}

template <typename Counter, typename Var>
/**
 * @brief Checks if the given variables are independent, given any
 *        subset of the given conditioning subset.
 *
 * @tparam Set The type of the container used for indices of the given variables.
 * @param x The index of the first variable.
 * @param y The index of the second variable.
 * @param given The indices of the variables to be conditioned on.
 * @param maxSize The maximum size of the subset to be tested.
 */
template <typename Set>
bool
DiscreteData<Counter, Var>::isIndependentAnySubset(
  const Var x,
  const Var y,
  const Set& given,
  const Var maxSize
) const
{
  auto minScore = this->minAssocScore(x, y, given, maxSize);
  return this->isIndependent(minScore);
}

template <typename Counter, typename Var>
/**
 * @brief Checks if the given variables are independent, given any
 *        subset of the given conditioning subset.
 *
 * @tparam Set The type of the container used for indices of the given variables.
 * @param x The index of the first variable.
 * @param y The index of the second variable.
 * @param given The indices of the variables to be conditioned on.
 * @param seed The indices of the variables to be included in every subset.
 * @param maxSize The maximum size of the subset to be tested.
 */
template <typename Set>
bool
DiscreteData<Counter, Var>::isIndependentAnySubset(
  const Var x,
  const Var y,
  const Set& given,
  const Set& seed,
  const Var maxSize
) const
{
  auto minScore = this->minAssocScore(x, y, given, seed, maxSize);
  return this->isIndependent(minScore);
}

template <typename Counter, typename Var>
/**
 * @brief Default destructor.
 */
DiscreteData<Counter, Var>::~DiscreteData(
)
{
  TIMER_ELAPSED_NONZERO("Time taken in G-square computations: ", m_timer);
}

#endif // DETAIL_DISCRETEDATA_HPP_
