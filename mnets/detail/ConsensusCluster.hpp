/**
 * @file ConsensusCluster.hpp
 * @brief Implementation of consensus clustering using Armadillo.
 * @author Sriram P. Chockalingam <srirampc@gatech.edu>
 * @version 1.0
 * @date 2021-01-05
 *
 * Copyright 2021 Georgia Institute of Technology
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
#ifndef DETAIL_CONSENSUSCLUSTER_HPP_
#define DETAIL_CONSENSUSCLUSTER_HPP_

#include "utils/Logging.hpp"

#include <armadillo>


template <typename MatType>
MatType
load_csv_mat(
  const char *file_name
)
{
  // load matrix from csv file
  MatType X;
  X.load(file_name, arma::csv_ascii);
  return X;
}

template <bool Norm, typename RandomIt, typename MatType>
double
clusterScore(const MatType&, const RandomIt, const RandomIt);

template <bool Norm, typename RandomIt>
double
clusterScore(
  const arma::mat& A,
  const RandomIt begin,
  const RandomIt end
)
{
  if (begin == end) {
    return 0.0;
  }
  arma::colvec v = arma::zeros<arma::colvec>(A.n_rows);
  for (auto it = begin; it != end; ++it) {
    v(*it) = 1.0;
  }
  auto w = A * v;
  auto score = arma::dot(v, w);
  auto denom = static_cast<double>(Norm ? pow(arma::norm(v), 2) :
                                          std::distance(begin, end));
  return score / denom;
}

template <bool Norm, typename RandomIt>
double
clusterScore(
  const arma::sp_mat& A,
  const RandomIt begin,
  const RandomIt end
)
{
  if (begin == end) {
    return 0.0;
  }
  arma::sp_mat v(A.n_rows, 1);
  for (auto it = begin; it != end; ++it) {
    v(*it, 0) = 1.0;
  }
  auto w = A * v;
  auto score = arma::dot(v, w);
  auto denom = static_cast<double>(Norm ? pow(arma::norm(v), 2) :
                                          std::distance(begin, end));
  return score / denom;
}

// Identify the cluster with best score by finding the top k
// vertices with the highest values of eigen vector that can maximize the
// cluster score
template <typename Var, typename MatType>
std::pair<std::vector<Var>, double>
bestCluster(
  const MatType& A,
  const arma::sp_mat&& v,
  const double tolerance
)
{
  LOG_MESSAGE_IF(v.n_cols != 1, error, "Perron vector should have one column");
  LOG_MESSAGE_IF(v.n_rows != A.n_rows, error, "Mismatch between the given row vector and matrix rows");
  LOG_MESSAGE_IF(v.n_rows != A.n_cols, error, "Mismatch between the given row vector and matrix columns");
  // Store non-zero values in sorted descending order
  std::set<std::pair<double, Var>, std::greater<std::pair<double, Var>>> vals;
  for (auto vit = v.begin_col(0); vit != v.end_col(0); ++vit) {
    // We only want to retain the elements which are greater than tolerance
    if (std::isgreater(*vit, tolerance)) {
      vals.insert(std::make_pair(*vit, vit.row()));
    }
  }
  // Now, copy all the eligible indices in the sorted order
  std::vector<Var> sortedIdx(vals.size());
  std::transform(vals.begin(), vals.end(), sortedIdx.begin(),
                 [] (const std::pair<Var, double>& vl)
                    { return vl.second; });

  auto maxScore = 0.0;
  auto numElements = 0u;
  // Find the top k vertices in the sorted order that maximzes the
  // cluster score
  const auto first = sortedIdx.cbegin();
  auto last = first + 1;
  for (auto k = 0u; k < sortedIdx.size(); ++k, ++last) {
    auto thisScore = clusterScore<true>(A, first, last);
    if (std::isgreaterequal(thisScore, maxScore)) {
      maxScore = thisScore;
      numElements = k + 1;
    }
  }
  LOG_MESSAGE_IF(numElements > 0, info, "Score: %g; Cluster Size: %u",
                                        clusterScore<false>(A, first, std::next(first, numElements)), numElements);
  sortedIdx.resize(numElements);
  return std::make_pair(sortedIdx, maxScore);
}

// Compute dominant eigen vector using the power method
template <typename MatType>
arma::sp_mat
perronVector(
  const MatType& A,
  const double tolerance,
  const uint32_t maxSteps
)
{
  // Initial vector : unit vector with 1.0 at vertex w. the maximum weight
  MatType rx = arma::sum(A, 1);
  arma::sp_mat v(rx.n_rows, 1);
  auto i = rx.index_max();
  v(i, 0) = 1.0;

  auto mu = 1.0;
  auto diff = 1.0;
  auto step = maxSteps;
  for (step = 0u; (step < maxSteps) && (diff > tolerance); ++step) {
    // Matrix - vector multiplication
    v = A * v;
    auto muNext = arma::norm(v);
    // normalize
    v = arma::normalise(v); // , p = 2);
    // Convergence parameter
    diff = std::abs(1.0 - muNext / mu);
    // Update mu
    mu = muNext;
  }
  LOG_MESSAGE(info, "Number of Steps / Max Steps: %u / %u", step, maxSteps);
  LOG_MESSAGE_IF(step == maxSteps, info, "Maximum number of steps reached with error = %g", diff);
  return v;
}


template <typename Var, typename MatType>
MatType
getSubmatrix(const MatType&, const std::vector<Var>&);

template <typename Var>
arma::mat
getSubmatrix(
  const arma::mat& Ain,
  const std::vector<Var>& indices
)
{
  // TODO : find a better way ?
  arma::mat Aout(indices.size(), indices.size());
  for (auto ix = 0u; ix < Aout.n_rows; ++ix) {
    for (auto jx = 0u; jx < Aout.n_cols; ++jx) {
      Aout(ix, jx) = Ain(indices[ix], indices[jx]);
    }
  }
  return Aout;
}

template <typename Var>
arma::sp_mat
getSubmatrix(
  const arma::sp_mat& Ain,
  const std::vector<Var>& indices
)
{
  // TODO : find a better way ?
  const auto nvertices = Ain.n_rows;
  arma::uvec map_idx(nvertices);
  map_idx.fill(nvertices);
  for (uint32_t ix = 0; ix < indices.size(); ++ix) {
    map_idx[indices[ix]] = ix;
  }
  uint32_t submat_size = 0;
  for (auto it = Ain.begin(); it != Ain.end(); ++it) {
    if (map_idx[it.row()] < nvertices && map_idx[it.col()] < nvertices) {
      ++submat_size;
    }
  }
  arma::umat locations(2, submat_size, arma::fill::zeros);
  arma::vec values(submat_size, arma::fill::zeros);
  uint32_t idx = 0;
  for (auto it = Ain.begin(); it != Ain.end(); ++it) {
    if ((map_idx[it.row()] < nvertices) && (map_idx[it.col()] < nvertices)) {
      locations(0, idx) = map_idx[it.row()];
      locations(1, idx) = map_idx[it.col()];
      values(idx) = *it; ++idx;
    }
  }
  return arma::sp_mat(locations, values, indices.size(), indices.size());
}

// Lemon tree cluster tightening algorithm
template <typename Var, typename MatType>
std::multimap<Var, Var>
perronCluster(
  MatType&& A,
  const double tolerance,
  const uint32_t maxSteps,
  const uint32_t minClustSize,
  const double minClustScore
)
{
  // set diagonal to ones as done by Lemon Tree
  A.diag().ones();
  // (to keep track of the original vertex id when we construct sub matrix)
  std::vector<Var> vertexMapping(A.n_rows);
  for (Var i = 0; i < A.n_rows; ++i) {
    vertexMapping[i] = i;
  }
  std::multimap<Var, Var> vertexClusters;
  Var clusterId = 0;
  auto numRemaining = A.n_rows;
  while (numRemaining >= minClustSize) {
    // Compute PF vector and best cluster
    auto v = perronVector(A, tolerance,  maxSteps);
    auto best = bestCluster<Var>(A, std::move(v), tolerance);
    auto& clusterElements = best.first;
    numRemaining = A.n_rows - clusterElements.size();
    // Identify remaining rows/cols
    std::vector<Var> currElements(A.n_rows);
    for (Var i = 0; i < A.n_rows; ++i) {
      currElements[i] = i;
    }
    std::sort(clusterElements.begin(), clusterElements.end());
    std::vector<Var> remainingElements(numRemaining);
    // XXX: We can possibly do better here
    std::set_difference(currElements.begin(), currElements.end(),
                        clusterElements.begin(), clusterElements.end(),
                        remainingElements.begin());
    A = getSubmatrix(A, remainingElements);
    if ((clusterElements.size() >= minClustSize) &&
        (best.second >= minClustScore)) {
      // Update the cluster ids of the clusters
      for (const auto ce : clusterElements) {
        auto vid = vertexMapping[ce];
        vertexClusters.emplace(clusterId, vid);
      }
      ++clusterId;
    }
    // Resize matrix and update mapping
    vertexMapping.resize(numRemaining);
    for (auto i = 0u; i < remainingElements.size(); ++i) {
      auto vid = vertexMapping[remainingElements[i]];
      vertexMapping[i] = vid;
    }
  }

  // XXX: Assign unique ids to the rest of the vertices ?
  return vertexClusters;
}

template <typename MatType, typename Set, typename Var>
void
fillCoclusteringWeights(
  MatType& C,
  const std::list<std::list<Set>>&& sampledClusters,
  const Var n,
  const double minWeight
)
{
  std::list<std::unordered_map<Var, uint32_t>> varClusterMaps;
  for (const auto& varClusters : sampledClusters) {
    std::unordered_map<Var, uint32_t> thisMap;
    auto c = 0u;
    for (auto cit = varClusters.begin(); cit != varClusters.end(); ++cit, ++c) {
      for (const auto var : *cit) {
        thisMap[var] = c;
      }
    }
    varClusterMaps.push_back(thisMap);
  }
  for (auto u = 0u; u < n; ++u) {
    for (auto v = u + 1; v < n; ++v) {
      auto cooccurrence = 0u;
      for (const auto& varCluster : varClusterMaps) {
        if (varCluster.at(u) == varCluster.at(v)) {
          ++cooccurrence;
        }
      }
      auto weight = static_cast<double>(cooccurrence) / sampledClusters.size();
      if (std::isgreater(weight, minWeight)) {
        C(u, v) = weight;
        C(v, u) = weight;
      }
    }
  }
}

template <typename Var, typename Set>
std::multimap<Var, Var>
consensusCluster(
  const std::list<std::list<Set>>&& sampledClusters,
  const Var numVars,
  const double minWeight,
  const double tolerance,
  const uint32_t maxSteps,
  const uint32_t minClustSize,
  const double minClustScore,
  const bool sparse = true
)
{
  if (sparse) {
    arma::sp_mat C(numVars, numVars);
    fillCoclusteringWeights(C, std::move(sampledClusters), numVars, minWeight);
    return perronCluster<Var>(std::move(C), tolerance, maxSteps, minClustSize, minClustScore);
  }
  else {
    arma::mat C(numVars, numVars, arma::fill::zeros);
    fillCoclusteringWeights(C, std::move(sampledClusters), numVars, minWeight);
    return perronCluster<Var>(std::move(C), tolerance, maxSteps, minClustSize, minClustScore);
  }
}

#endif // DETAIL_CONSENSUSCLUSTER_HPP_
