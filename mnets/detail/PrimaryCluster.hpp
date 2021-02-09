/**
 * @file PrimaryCluster.hpp
 * @brief Implementation of functionality for storing primary clusters.
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
#ifndef DETAIL_PRIMARYCLUSTER_HPP_
#define DETAIL_PRIMARYCLUSTER_HPP_

#include "Cluster.hpp"
#include "SecondaryCluster.hpp"


/**
 * @brief Class that provides functionality for storing
 *        primary clusters and computing their score.
 *
 * @tparam Data Type of the data provider.
 * @tparam Var Type of variables stored in the cluster.
 * @tparam Set Type of container used to store the clusters.
 */
template <typename Data, typename Var, typename Set>
class PrimaryCluster : public Cluster<Data, Var, Set> {
public:
  PrimaryCluster(const Data&, std::mt19937* const, const Var, const Var);

  PrimaryCluster(const Data&, std::mt19937* const, const Set&, const Var);

  PrimaryCluster(const PrimaryCluster&);

  PrimaryCluster(const PrimaryCluster&, const PrimaryCluster&);

  ~PrimaryCluster();

  double
  score();

  void
  scoreClear();

  double
  scoreMerge(const PrimaryCluster&, const bool = false);

  double
  scoreInsertPrimary(const Var, const bool = false);

  double
  scoreErasePrimary(const Var, const bool = false);

  void
  clearSecondary();

  void
  singleSecondary();

  void
  randomSecondary(const Var);

  void
  clusterSecondary(const uint32_t = 1);

  const std::list<SecondaryCluster<Data, Var, Set>>&
  secondaryClusters() const;

private:
  void
  removeEmptyClusters();

  void
  reassignSecondary(const Var);

  bool
  mergeCluster(typename std::list<SecondaryCluster<Data, Var, Set>>::iterator&);

private:
  std::list<SecondaryCluster<Data, Var, Set>> m_cluster;
  std::vector<typename std::list<SecondaryCluster<Data, Var, Set>>::iterator> m_membership;
  double m_score;
  const Var m_numSecondary;
}; // class PrimaryCluster

template <typename Data, typename Var, typename Set>
/**
 * @brief Constructs an empty primary cluster.
 *
 * @param data The data provider.
 * @param generator The random number generator.
 * @param numPrimary Number of primary variables.
 * @param numSecondary Number of secondary variables.
 */
PrimaryCluster<Data, Var, Set>::PrimaryCluster(
  const Data& data,
  std::mt19937* const generator,
  const Var numPrimary,
  const Var numSecondary
) : Cluster<Data, Var, Set>(data, generator, numPrimary),
    m_cluster(),
    m_membership(numSecondary, m_cluster.end()),
    m_score(std::nan("")),
    m_numSecondary(numSecondary)
{
}

template <typename Data, typename Var, typename Set>
/**
 * @brief Constructs a primary cluster from the given elements.
 *
 * @param data The data provider.
 * @param generator The random number generator.
 * @param primaryElements The primary variables in the cluster.
 * @param numSecondary Number of secondary variables.
 */
PrimaryCluster<Data, Var, Set>::PrimaryCluster(
  const Data& data,
  std::mt19937* const generator,
  const Set& primaryElements,
  const Var numSecondary
) : Cluster<Data, Var, Set>(data, generator, primaryElements),
    m_cluster(),
    m_membership(numSecondary, m_cluster.end()),
    m_score(std::nan("")),
    m_numSecondary(numSecondary)
{
}

template <typename Data, typename Var, typename Set>
/**
 * @brief Copy constructor.
 *
 * @param other The primary cluster to be copied.
 */
PrimaryCluster<Data, Var, Set>::PrimaryCluster(
  const PrimaryCluster<Data, Var, Set>& other
) : Cluster<Data, Var, Set>(other),
    m_cluster(),
    m_membership(other.m_numSecondary, m_cluster.end()),
    m_score(other.m_score),
    m_numSecondary(other.m_numSecondary)
{
  // Create a copy of the secondary clusters
  for (auto cit = other.m_cluster.begin(); cit != other.m_cluster.end(); ++cit) {
    m_cluster.emplace_back(*cit);
    for (const auto e : cit->elements()) {
      m_membership[e] = std::prev(m_cluster.end());
    }
  }
}

template <typename Data, typename Var, typename Set>
/**
 * @brief Merge constructor creates a new primary cluster
 *        by merging the two given clusters.
 *
 * @param first The first primary cluster to be merged.
 * @param second The second primary cluster to be merged.
 */
PrimaryCluster<Data, Var, Set>::PrimaryCluster(
  const PrimaryCluster<Data, Var, Set>& first,
  const PrimaryCluster<Data, Var, Set>& second
) : Cluster<Data, Var, Set>(first, second),
    m_cluster(),
    m_membership(first.m_numSecondary, m_cluster.end()),
    m_score(std::nan("")),
    m_numSecondary(first.m_numSecondary)
{
  // Just one secondary cluster with all the elements is created
  this->singleSecondary();
}

template <typename Data, typename Var, typename Set>
/**
 * @brief Default destructor.
 */
PrimaryCluster<Data, Var, Set>::~PrimaryCluster(
)
{
}

template <typename Data, typename Var, typename Set>
/**
 * @brief Computes the score of this cluster, if not cached,
 *        and returns it.
 */
double
PrimaryCluster<Data, Var, Set>::score(
)
{
  if (std::isnan(m_score)) {
    m_score = 0.0;
    // Iterate over all the secondary clusters and
    // compute the score for each one separately
    for (auto& secondary : m_cluster) {
      m_score += secondary.score(*this);
    }
  }
  return m_score;
}

template <typename Data, typename Var, typename Set>
/**
 * @brief Clears the cached score for this cluster.
 */
void
PrimaryCluster<Data, Var, Set>::scoreClear(
)
{
  m_score = std::nan("");
}

template <typename Data, typename Var, typename Set>
/**
 * @brief Computes the score of this primary cluster when a primary variable
 *        is inserted, optionally updating the cached score.
 *
 * @param given The index of the primary variable to be inserted.
 * @param cache If the cached score should be updated.
 *
 * @return The changed score of this cluster after inserting the variable.
 */
double
PrimaryCluster<Data, Var, Set>::scoreInsertPrimary(
  const Var given,
  const bool cache
)
{
  auto score = 0.0;
  for (auto& secondary : m_cluster) {
    score += secondary.scoreInsertPrimary(*this, given, cache);
  }
  if (cache) {
    m_score = score;
  }
  return score;
}

template <typename Data, typename Var, typename Set>
/**
 * @brief Computes the score of this cluster when a primary variable
 *        is erased, optionally updating the cached score.
 *
 * @param given The index of the primary variable to be erased.
 * @param cache If the cached score should be updated.
 *
 * @return The changed score of this cluster after erasing the variable.
 */
double
PrimaryCluster<Data, Var, Set>::scoreErasePrimary(
  const Var given,
  const bool cache
)
{
  auto score = 0.0;
  for (auto& secondary : m_cluster) {
    score += secondary.scoreErasePrimary(*this, given, cache);
  }
  if (cache) {
    m_score = score;
  }
  return score;
}

template <typename Data, typename Var, typename Set>
/**
 * @brief Computes the score of this cluster when another primary cluster
 *        is merged with it, optionally updating the cached score.
 *
 * @param other The primary cluster to be merged.
 * @param cache If the cached score should be updated.
 *
 * @return The changed score of this cluster after merging the clusters.
 */
double
PrimaryCluster<Data, Var, Set>::scoreMerge(
  const PrimaryCluster<Data, Var, Set>& other,
  const bool cache
)
{
  auto score = 0.0;
  // Iterate over all the secondary clusters and
  // compute the score for each one separately
  for (auto& secondary : m_cluster) {
    score += secondary.scoreInsertPrimary(*this, other.elements(), cache);
  }
  if (cache) {
    m_score = score;
  }
  return score;
}

template <typename Data, typename Var, typename Set>
/**
 * @brief Clears all the secondary clusters for this primary cluster.
 *        Also clears the cached score for this cluster.
 */
void
PrimaryCluster<Data, Var, Set>::clearSecondary(
)
{
  LOG_MESSAGE(info, "Clearing all secondary clusters");
  for (auto& m : m_membership) {
    m = m_cluster.end();
  }
  m_cluster.clear();
  this->scoreClear();
}

template <typename Data, typename Var, typename Set>
/**
 * @brief Creates a single secondary cluster for this primary cluster.
 *        Also clears the cached score for this cluster.
 */
void
PrimaryCluster<Data, Var, Set>::singleSecondary(
)
{
  LOG_MESSAGE(info, "Assigning all secondary variables to the same cluster");
  m_cluster.clear();
  m_cluster.emplace_back(this->m_data, m_numSecondary);
  auto single = m_cluster.begin();
  for (Var e = 0u; e < m_numSecondary; ++e) {
    single->insert(e);
    m_membership[e] = single;
  }
  this->scoreClear();
}

template <typename Data, typename Var, typename Set>
/**
 * @brief Randomly initializes the secondary clusters for this primary cluster.
 *        Also clears the cached score for this cluster.
 */
void
PrimaryCluster<Data, Var, Set>::randomSecondary(
  const Var numClusters
)
{
  LOG_MESSAGE(info, "Randomly assigning secondary variables to %u clusters", static_cast<uint32_t>(numClusters));
  std::vector<SecondaryCluster<Data, Var, Set>> cluster(numClusters, SecondaryCluster<Data, Var, Set>(this->m_data, m_numSecondary));
  std::uniform_int_distribution<Var> clusterDistrib(0, numClusters - 1);
  for (Var e = 0u; e < m_numSecondary; ++e) {
    auto c = clusterDistrib(*(this->m_generator));
    cluster[c].insert(e);
  }
  m_cluster = std::list<SecondaryCluster<Data, Var, Set>>(cluster.begin(), cluster.end());
  this->removeEmptyClusters();
  for (auto cit = m_cluster.begin(); cit != m_cluster.end(); ++cit) {
    for (const auto e : cit->elements()) {
      m_membership[e] = cit;
    }
  }
  this->scoreClear();
  LOG_MESSAGE(info, "Assigned secondary variables to %u clusters", m_cluster.size());
}

template <typename Data, typename Var, typename Set>
/**
 * @brief Performs a Gibbs clustering step for the secondary variables
 *        corresponding to this primary cluster.
 *
 * @param numReps Number of times clustering of secondary variables
 *                should be repeated.
 */
void
PrimaryCluster<Data, Var, Set>::clusterSecondary(
  const uint32_t numReps
)
{
  std::uniform_int_distribution<Var> varDistrib(0, m_numSecondary - 1);
  for (auto r = 0u; r < numReps; ++r) {
    // Reassign a random secondary variable for n iterations
    LOG_MESSAGE(info, "Reassigning secondary variables");
    for (auto i = 0u; i < m_numSecondary; ++i) {
      auto v = varDistrib(*(this->m_generator));
      this->reassignSecondary(v);
    }
    LOG_MESSAGE(info, "Done reassigning secondary variables");
    LOG_MESSAGE(info, "Merging secondary clusters (number of clusters = %u)", m_cluster.size());
    for (auto cit = m_cluster.begin(); (cit != m_cluster.end()) && (m_cluster.size() > 1); ) {
      if (this->mergeCluster(cit)) {
        cit = m_cluster.erase(cit);
      }
      else {
        ++cit;
      }
    }
    LOG_MESSAGE(info, "Done merging secondary clusters (number of clusters = %u)", m_cluster.size());
  }
  this->scoreClear();
}

template <typename Data, typename Var, typename Set>
/**
 * @brief Returns the secondary clusters corresponding to this primary cluster.
 */
const std::list<SecondaryCluster<Data, Var, Set>>&
PrimaryCluster<Data, Var, Set>::secondaryClusters(
) const
{
  return m_cluster;
}

template <typename Data, typename Var, typename Set>
/**
 * @brief Removes all the empty secondary clusters from this primary cluster.
 */
void
PrimaryCluster<Data, Var, Set>::removeEmptyClusters(
)
{
  auto emptyCluster = [] (const SecondaryCluster<Data, Var, Set>& cluster)
                         { return cluster.empty(); };
  m_cluster.remove_if(emptyCluster);
}

template <typename Data, typename Var, typename Set>
/**
 * @brief Moves the given secondary variable to a different
 *        secondary cluster in this primary cluster.
 *
 * @param given The index of the secondary variable to be moved.
 */
void
PrimaryCluster<Data, Var, Set>::reassignSecondary(
  const Var given
)
{
  LOG_MESSAGE(debug, "Reassigning secondary variable %u", static_cast<uint32_t>(given));
  // Remove the given var from the old cluster
  auto oldCluster = m_membership[given];
  m_membership[given] = m_cluster.end();
  // Create a new cluster with only the given var
  SecondaryCluster<Data, Var, Set> newCluster(this->m_data, m_numSecondary);
  newCluster.insert(given);
  if (oldCluster->size() > 1) {
    // Remove the element and update the score of the cluster
    oldCluster->scoreEraseSecondary(*this, given, true);
    oldCluster->erase(given);
  }
  else {
    // Remove the cluster if this variable was its only element
    LOG_MESSAGE(debug, "Removing the old cluster of the variable");
    m_cluster.erase(oldCluster);
  }
  // Compute the weight of the var being in its separate cluster
  // as well as all the existing clusters
  std::vector<double> weight(m_cluster.size() + 1);
  weight[0] = 1.0;
  auto wit = weight.begin() + 1;
  auto maxDiff = weight[0];
  auto singleScore = newCluster.score(*this);
  // Only compute score diffs for existing clusters
  for (auto cit = m_cluster.begin(); cit != m_cluster.end(); ++cit, ++wit) {
    auto thisDiff = cit->scoreInsertSecondary(*this, given) -
                    cit->score(*this) - singleScore;
    *wit = thisDiff;
    maxDiff = std::max(thisDiff, maxDiff);
  }
  for (auto& w : weight) {
    w = exp(w - maxDiff);
  }
  // Pick a cluster using the computed weights
  std::discrete_distribution<Var> clusterDistrib(weight.begin(), weight.end());
  auto c = clusterDistrib(*(this->m_generator));
  if (c == 0) {
    // The variable will stay in its own cluster
    LOG_MESSAGE(info, "Secondary variable %u assigned to a newly created cluster", static_cast<uint32_t>(given));
    m_cluster.push_back(std::move(newCluster));
    m_membership[given] = std::prev(m_cluster.end());
  }
  else {
    // Add the variable to the chosen cluster
    LOG_MESSAGE(info, "Secondary variable %u assigned to the existing cluster %u",
                      static_cast<uint32_t>(given), static_cast<uint32_t>(c - 1));
    auto chosen = std::next(m_cluster.begin(), c - 1);
    chosen->scoreInsertSecondary(*this, given, true);
    chosen->insert(given);
    m_membership[given] = chosen;
  }
}

template <typename Data, typename Var, typename Set>
/**
 * @brief Merges the given secondary cluster with another secondary cluster
 *        in this primary cluster.
 */
bool
PrimaryCluster<Data, Var, Set>::mergeCluster(
  typename std::list<SecondaryCluster<Data, Var, Set>>::iterator& given
)
{
  // Compute the weight of merging this cluster with
  // all the other clusters which are not empty
  auto givenScore = given->score(*this);
  std::vector<double> weight(m_cluster.size(), 0.0);
  auto wit = weight.begin();
  for (auto cit = m_cluster.begin(); cit != m_cluster.end(); ++cit, ++wit) {
    if (cit != given) {
      auto thisDiff = cit->scoreMerge(*this, *given) -
                      cit->score(*this) - givenScore;
      *wit = exp(thisDiff);
    }
    else {
      *wit = 1.0;
    }
  }
  // Choose a cluster using the computed weights
  std::discrete_distribution<Var> clusterDistrib(weight.begin(), weight.end());
  auto c = clusterDistrib(*(this->m_generator));
  auto chosen = std::next(m_cluster.begin(), c);
  if (chosen != given) {
    LOG_MESSAGE(info, "Merging given cluster with cluster %u", static_cast<uint32_t>(c));
    // Merge this cluster with the chosen cluster
    // and update the membership of all the moved elements
    for (const auto e : given->elements()) {
      m_membership[e] = chosen;
    }
    chosen->scoreMerge(*this, *given, true);
    chosen->merge(*given);
    return true;
  }
  else {
    LOG_MESSAGE(info, "Not merging given cluster");
    return false;
  }
}

#endif // DETAIL_PRIMARYCLUSTER_HPP_
