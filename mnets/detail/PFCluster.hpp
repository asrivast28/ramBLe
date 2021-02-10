/**
 * @file PFCluster.hpp
 * @brief Implementation of tighten clusters using Armadillo.
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
#ifndef DETAIL_PFCLUSTER_HPP_
#define DETAIL_PFCLUSTER_HPP_

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

// Lemon tree has too different ways to score a cluster
// 1. Using average
// 2. Using norm defn.
template <typename MatType, typename RandomIt>
double
clusterScore1(
  const MatType& A,
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
  return score / std::distance(begin, end);
}

template <typename MatType, typename RandomIt>
double
clusterScore2(
  const MatType& A,
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
  return score / pow(arma::norm(v), 2);
}

// Identify the cluster with best score by finding the top k
// vertices with the highest values of eigen vector that can maximize the
// cluster score
template <typename Var, typename MatType>
std::pair<std::vector<Var>, double>
bestCluster(
  const MatType& A,
  const arma::vec& v,
  const double tolerance
)
{
  auto nvals = v.n_elem;
  LOG_MESSAGE_IF(nvals != A.n_rows, error, "Mismatch between the given row vector and matrix rows");
  LOG_MESSAGE_IF(nvals != A.n_cols, error, "Mismatch between the given row vector and matrix columns");
  std::vector<std::pair<double, Var>> vals(nvals);
  for (Var i = 0; i < nvals; ++i) {
    vals[i] = std::make_pair(v[i], i);
  }
  // Sort the eigen vector in descending order
  std::sort(vals.begin(), vals.end(), std::greater<std::pair<double, Var>>());
  // We only want to retain the elements which are greater than tolerance
  // Find iterator to the last element which is greater than tolerance
  auto it = std::partition_point(vals.begin(), vals.end(),
                                 [&tolerance] (const std::pair<double, Var>& v)
                                              { return std::isgreaterequal(v.first, tolerance); });
  // Now, copy all the eligible indices in sorted order
  std::vector<Var> sortedIdx(std::distance(vals.begin(), it));
  std::transform(vals.begin(), it, sortedIdx.begin(),
                 [] (const std::pair<Var, double>& v)
                    { return v.second; });

  auto maxScore = 0.0;
  auto numElements = 0u;
  // Find the top k vertices in the sorted order that maximzes the
  // cluster score
  const auto first = sortedIdx.cbegin();
  auto last = first + 1;
  for (auto k = 0u; k < sortedIdx.size(); ++k, ++last) {
    auto thisScore = clusterScore2(A, first, last);
    if (std::isgreater(thisScore, maxScore)) {
      maxScore = thisScore;
      numElements = k + 1;
    }
  }
  LOG_MESSAGE_IF(numElements > 0, info, "Cluster Size: %u; Score: %g",
                                        numElements, clusterScore1(A, first, std::next(first, numElements)));
  sortedIdx.resize(numElements);
  return std::make_pair(sortedIdx, maxScore);
}

// Compute dominant eigen vector using the power method
template <typename MatType>
arma::vec
perronVector(
  const MatType& A,
  const double tolerance,
  const uint32_t maxSteps
)
{
  // Initial vector : unit vector with 1.0 at vertex w. the maximum weight
  MatType rx = arma::sum(A, 1);
  arma::colvec v = arma::zeros<arma::colvec>(rx.n_rows);
  auto i = rx.index_max();
  v(i) = 1.0;

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
getSubmatrix(const MatType& Ain, const std::vector<Var>&);

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
    auto best = bestCluster<Var>(A, v, tolerance);
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
      // Resize matrix and update mapping
      vertexMapping.resize(numRemaining);
      for (auto i = 0u; i < remainingElements.size(); ++i) {
        auto vid = vertexMapping[remainingElements[i]];
        vertexMapping[i] = vid;
      }
      ++clusterId;
    }
  }

  // XXX: Assign unique ids to the rest of the vertices ?
  return vertexClusters;
}

#if 0
std::vector<VertClusterId>
sparse_perron(
  const std::vector<double>& coMatrix,
  const uint32_t numVars,
  const double tolerance,
  const uint32_t maxSteps,
  const uint32_t minClustSize,
  const double minClustScore
)
{
  arma::sp_mat C = load_csv_mat<arma::sp_mat>(p.mat_file.c_str());
  return perronCluster<Var>(std::move(C), tolerance, maxSteps, minClustSize, minClustScore);
}
#endif

template <typename Var>
std::multimap<Var, Var>
densePerron(
  const std::vector<double>&& coMatrix,
  const Var numVars,
  const double tolerance,
  const uint32_t maxSteps,
  const uint32_t minClustSize,
  const double minClustScore
)
{
  arma::mat C(&coMatrix[0], numVars, numVars);
  return perronCluster<Var>(std::move(C), tolerance, maxSteps, minClustSize, minClustScore);
}

#endif // DETAIL_PFCLUSTER_HPP_
