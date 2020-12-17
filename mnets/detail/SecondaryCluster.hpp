/**
 * @file SecondaryCluster.hpp
 * @brief Implementation of functionality for storing secondary clusters.
 * @author Ankit Srivastava <asrivast@gatech.edu>
 *
 * Copyright 2020 Georgia Institute of Technology
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
#ifndef DETAIL_SECONDARYCLUSTER_HPP_
#define DETAIL_SECONDARYCLUSTER_HPP_

#include "Cluster.hpp"


/**
 * @brief Computes the log-likelihood from the given statistics.
 *
 * @param count Number of data points.
 * @param sum Sum of the data points.
 * @param sum2 Sum of squares of the data points.
 */
double
computeLogLikelihood(
  const uint32_t count,
  const double sum,
  const double sum2
)
{
  // Fixed parameters
  static constexpr double lambda0 = 0.1;
  static constexpr double alpha0 = 0.1;
  static constexpr double beta0 = 0.1;
  static constexpr double mu = 0.0;
  static constexpr double log2pi = log(2 * acos(-1));
  // Function parameter independent computations
  // lgamma is not a constexpr
  static const auto fixed = 0.5 * log(lambda0) + alpha0 * log(beta0) - lgamma(alpha0);
  // Dependent parameters
  const auto lambda1 = lambda0 + count;
  const auto alpha1 = alpha0 + 0.5 * count;
  const auto beta1 = beta0 + 0.5 * (sum2 - pow(sum, 2) / count) +
                     lambda0 * pow(sum - mu * count, 2) / (2 * lambda1 * count);
  // Log-likelihood computation
  auto logLikelihood = -0.5 * count * log2pi + fixed + lgamma(alpha1) -
                       alpha1 * log(beta1) - 0.5 * log(lambda1);
  return std::isnan(logLikelihood) ? 0.0 : logLikelihood;
}

/**
 * @brief Class that provides functionality for storing
 *        primary clusters and computing their score.
 *
 * @tparam Data Type of the data provider.
 * @tparam Var Type of variables stored in the cluster.
 * @tparam Set Type of container used to store the clusters.
 */
template <typename Data, typename Var, typename Set>
class SecondaryCluster : public Cluster<Data, Var, Set> {
public:
  SecondaryCluster(const Data&, const Var);

  SecondaryCluster(const SecondaryCluster&);

  SecondaryCluster(const SecondaryCluster&, const SecondaryCluster&);

  ~SecondaryCluster();

  double
  score(const Cluster<Data, Var, Set>&);

  double
  scoreMerge(const Cluster<Data, Var, Set>&, const SecondaryCluster&, const bool = false);

  double
  scoreInsertPrimary(const Cluster<Data, Var, Set>&, const Var, const bool = false);

  double
  scoreInsertPrimary(const Cluster<Data, Var, Set>&, const Set&, const bool = false);

  double
  scoreInsertSecondary(const Cluster<Data, Var, Set>&, const Var, const bool = false);

  double
  scoreErasePrimary(const Cluster<Data, Var, Set>&, const Var, const bool = false);

  double
  scoreEraseSecondary(const Cluster<Data, Var, Set>&, const Var, const bool = false);

private:
  void
  scoreCache(const Cluster<Data, Var, Set>&);

private:
  Set m_primary;
  double m_score;
  double m_sum;
  double m_sum2;
  uint32_t m_count;
};

template <typename Data, typename Var, typename Set>
/**
 * @brief Constructs an empty secondary cluster.
 *
 * @param data The data provider.
 * @param numSecondary Number of secondary variables.
 */
SecondaryCluster<Data, Var, Set>::SecondaryCluster(
  const Data& data,
  const Var numSecondary
) : Cluster<Data, Var, Set>(data, nullptr, numSecondary),
    m_primary(),
    m_score(std::nan("")),
    m_sum(),
    m_sum2(),
    m_count()
{
}

template <typename Data, typename Var, typename Set>
/**
 * @brief Copy constructor.
 *
 * @param other The secondary cluster to be copied.
 */
SecondaryCluster<Data, Var, Set>::SecondaryCluster(
  const SecondaryCluster<Data, Var, Set>& other
) : Cluster<Data, Var, Set>(other),
    m_primary(other.m_primary),
    m_score(other.m_score),
    m_sum(other.m_sum),
    m_sum2(other.m_sum2),
    m_count(other.m_count)
{
}

template <typename Data, typename Var, typename Set>
/**
 * @brief Merge constructor creates a new secondary cluster
 *        by merging the two given clusters.
 *
 * @param first The first secondary cluster to be merged.
 * @param second The second secondary cluster to be merged.
 */
SecondaryCluster<Data, Var, Set>::SecondaryCluster(
  const SecondaryCluster<Data, Var, Set>& first,
  const SecondaryCluster<Data, Var, Set>& second
) : Cluster<Data, Var, Set>(first, second),
    m_primary(first.m_primary),
    m_score(std::nan("")),
    m_sum(first.m_sum + second.m_sum),
    m_sum2(first.m_sum2 + second.m_sum2),
    m_count(first.m_count + second.m_count)
{
  LOG_MESSAGE_IF(first.m_primary != second.m_primary, error,
                 "Merging secondary clusters with different primary clusters");
  m_score = computeLogLikelihood(m_count, m_sum, m_sum2);
}

template <typename Data, typename Var, typename Set>
/**
 * @brief Default destructor.
 */
SecondaryCluster<Data, Var, Set>::~SecondaryCluster(
)
{
}

template <typename Data, typename Var, typename Set>
/**
 * @brief Caches the score of this cluster corresponding
 *        to the given elements in the primary cluster.
 *
 * @param primary The primary cluster to be used for score computations.
 */
void
SecondaryCluster<Data, Var, Set>::scoreCache(
  const Cluster<Data, Var, Set>& primary
)
{
  // We need to compute the score if it has never been computed
  // before or if the elements in the primary cluster have changed
  if (std::isnan(m_score) || (m_primary != primary.elements())) {
    if (!m_primary.empty()) {
      // Reset counts for previous data
      m_sum = 0.0;
      m_sum2 = 0.0;
      m_count = 0u;
    }
    // Store a snapshot of the elements in the primary cluster
    m_primary = primary.elements();
    for (const auto p : m_primary) {
      for (const auto s : this->m_elements) {
        auto d = this->m_data(p, s);
        if (!std::isnan(d)) {
          m_sum += d;
          m_sum2 += d * d;
          ++m_count;
        }
      }
    }
    m_score = computeLogLikelihood(m_count, m_sum, m_sum2);
  }
}

template <typename Data, typename Var, typename Set>
/**
 * @brief Returns the score of this secondary cluster corresponding
 *        to the given primary cluster elements.
 *
 * @param primary The primary cluster to be used for score computations.
 */
double
SecondaryCluster<Data, Var, Set>::score(
  const Cluster<Data, Var, Set>& primary
)
{
  this->scoreCache(primary);
  return m_score;
}

template <typename Data, typename Var, typename Set>
/**
 * @brief Computes the score of this cluster when another secondary cluster
 *        is merged with it, optionally updating the cached score.
 *
 * @param primary The primary cluster to be used for score computations.
 * @param other The secondary cluster to be merged.
 * @param cache If the cached score should be updated.
 *
 * @return The changed score of this cluster after merging the clusters.
 */
double
SecondaryCluster<Data, Var, Set>::scoreMerge(
  const Cluster<Data, Var, Set>& primary,
  const SecondaryCluster<Data, Var, Set>& other,
  const bool cache
)
{
  this->scoreCache(primary);
  LOG_MESSAGE_IF(m_primary != other.m_primary, error,
                 "Merging secondary clusters with different primary clusters");
  auto score = computeLogLikelihood(m_count + other.m_count,
                                    m_sum + other.m_sum,
                                    m_sum2 + other.m_sum2);
  if (cache) {
    m_sum += other.m_sum;
    m_sum2 += other.m_sum2;
    m_count += other.m_count;
    m_score = score;
  }
  return score;
}

template <typename Data, typename Var, typename Set>
/**
 * @brief Computes the score of this secondary cluster when a primary variable
 *        is inserted, optionally updating the cached score.
 *
 * @param primary The primary cluster to be used for score computations.
 * @param given The index of the primary variable to be inserted.
 * @param cache If the cached score should be updated.
 *
 * @return The changed score of this cluster after inserting the variable.
 */
double
SecondaryCluster<Data, Var, Set>::scoreInsertPrimary(
  const Cluster<Data, Var, Set>& primary,
  const Var given,
  const bool cache
)
{
  this->scoreCache(primary);
  auto sum = 0.0;
  auto sum2 = 0.0;
  auto count = 0;
  for (const auto s : this->m_elements) {
    auto d = this->m_data(given, s);
    if (!std::isnan(d)) {
      sum += d;
      sum2 += d * d;
      ++count;
    }
  }
  auto score = computeLogLikelihood(m_count + count, m_sum + sum, m_sum2 + sum2);
  if (cache) {
    m_primary.insert(given);
    m_sum += sum;
    m_sum2 += sum2;
    m_count += count;
    m_score = score;
  }
  return score;
}

template <typename Data, typename Var, typename Set>
/**
 * @brief Computes the score of this secondary cluster when multiple primary variables
 *        are inserted, optionally updating the cached score.
 *
 * @param primary The primary cluster to be used for score computations.
 * @param newElements A container with the elements to be inserted in the primary cluster.
 * @param cache If the cached score should be updated.
 *
 * @return The changed score of this cluster after inserting the variables.
 */
double
SecondaryCluster<Data, Var, Set>::scoreInsertPrimary(
  const Cluster<Data, Var, Set>& primary,
  const Set& newElements,
  const bool cache
)
{
  this->scoreCache(primary);
  auto sum = 0.0;
  auto sum2 = 0.0;
  auto count = 0;
  // We assume that the sets are mutually exclusive
  for (const auto p : newElements) {
    for (const auto s : this->m_elements) {
      auto d = this->m_data(p, s);
      if (!std::isnan(d)) {
        sum += d;
        sum2 += d * d;
        ++count;
      }
    }
  }
  auto score = computeLogLikelihood(m_count + count, m_sum + sum, m_sum2 + sum2);
  if (cache) {
    m_primary = set_union(primary.elements(), newElements);
    m_sum += sum;
    m_sum2 += sum2;
    m_count += count;
    m_score = score;
  }
  return score;
}

template <typename Data, typename Var, typename Set>
/**
 * @brief Computes the score of this primary cluster when a secondary variable
 *        is inserted, optionally updating the cached score.
 *
 * @param primary The primary cluster to be used for score computations.
 * @param given The index of the secondary variable to be inserted.
 * @param cache If the cached score should be updated.
 *
 * @return The changed score of this cluster after inserting the variable.
 */
double
SecondaryCluster<Data, Var, Set>::scoreInsertSecondary(
  const Cluster<Data, Var, Set>& primary,
  const Var given,
  const bool cache
)
{
  this->scoreCache(primary);
  auto sum = 0.0;
  auto sum2 = 0.0;
  auto count = 0;
  for (const auto p : m_primary) {
    auto d = this->m_data(p, given);
    if (!std::isnan(d)) {
      sum += d;
      sum2 += d * d;
      ++count;
    }
  }
  auto score = computeLogLikelihood(m_count + count, m_sum + sum, m_sum2 + sum2);
  if (cache) {
    m_sum += sum;
    m_sum2 += sum2;
    m_count += count;
    m_score = score;
  }
  return score;
}

template <typename Data, typename Var, typename Set>
/**
 * @brief Computes the score of this cluster when a primary variable
 *        is erased, optionally updating the cached score.
 *
 * @param primary The primary cluster to be used for score computations.
 * @param given The index of the primary variable to be erased.
 * @param cache If the cached score should be updated.
 *
 * @return The changed score of this cluster after erasing the variable.
 */
double
SecondaryCluster<Data, Var, Set>::scoreErasePrimary(
  const Cluster<Data, Var, Set>& primary,
  const Var given,
  const bool cache
)
{
  this->scoreCache(primary);
  auto sum = 0.0;
  auto sum2 = 0.0;
  auto count = 0;
  for (const auto s : this->m_elements) {
    auto d = this->m_data(given, s);
    if (!std::isnan(d)) {
      sum += d;
      sum2 += d * d;
      ++count;
    }
  }
  auto score = computeLogLikelihood(m_count - count, m_sum - sum, m_sum2 - sum2);
  if (cache) {
    m_primary.erase(given);
    m_sum -= sum;
    m_sum2 -= sum2;
    m_count -= count;
    m_score = score;
  }
  return score;
}

template <typename Data, typename Var, typename Set>
/**
 * @brief Computes the score of this cluster when a secondary variable
 *        is erased, optionally updating the cached score.
 *
 * @param primary The primary cluster to be used for score computations.
 * @param given The index of the secondary variable to be erased.
 * @param cache If the cached score should be updated.
 *
 * @return The changed score of this cluster after erasing the variable.
 */
double
SecondaryCluster<Data, Var, Set>::scoreEraseSecondary(
  const Cluster<Data, Var, Set>& primary,
  const Var given,
  const bool cache
)
{
  this->scoreCache(primary);
  auto sum = 0.0;
  auto sum2 = 0.0;
  auto count = 0;
  for (const auto p : m_primary) {
    auto d = this->m_data(p, given);
    if (!std::isnan(d)) {
      sum += d;
      sum2 += d * d;
      ++count;
    }
  }
  auto score = computeLogLikelihood(m_count - count, m_sum - sum, m_sum2 - sum2);
  if (cache) {
    m_sum -= sum;
    m_sum2 -= sum2;
    m_count -= count;
    m_score = score;
  }
  return score;
}

#endif // DETAIL_SECONDARYCLUSTER_HPP_
