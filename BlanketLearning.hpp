/**
 * @file BlanketLearning.hpp
 * @brief Declaration of the classes for blanket learning algorithms.
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
#ifndef BLANKETLEARNING_HPP_
#define BLANKETLEARNING_HPP_

#include "LocalLearning.hpp"

#include "utils/Timer.hpp"


/**
 * @brief Abstract base class for causal discovery by first learning MB sets.
 *
 * @tparam Data Type of the object which is used for querying the data.
 * @tparam Var Type of variable indices (expected to be an integer type).
 * @tparam Set Type of set container.
 */
template <typename Data, typename Var, typename Set>
class BlanketLearning : public LocalLearning<Data, Var, Set> {
public:
  BlanketLearning(const mxx::comm&, const Data&, const Var);

  virtual
  ~BlanketLearning();

protected:
  Set
  shrinkMB(const Var, Set&) const;

  virtual
  std::pair<Var, double>
  pickBestCandidate(const Var, const Set&, const Set&) const;

  virtual
  void
  updatePValues(std::vector<std::tuple<Var, Var, double>>&, const std::unordered_map<Var, Set>&) const;

  std::set<std::tuple<Var, Var, double>>
  growAll(const std::vector<std::tuple<Var, Var, double>>&, std::unordered_map<Var, Set>&) const;

  std::set<std::pair<Var, Var>>
  shrinkAll(std::unordered_map<Var, Set>&) const;

  virtual
  void
  growShrink(std::vector<std::tuple<Var, Var, double>>&&, std::unordered_map<Var, Set>&, std::set<std::pair<Var, Var>>&, const double) const;

private:
  Set
  getCandidatePC(const Var, Set&&) const override;

  bool
  evaluateCandidatePC(const Var, const Var, const Set&, const Set&) const;

  BayesianNetwork<Var>
  getSkeleton_parallel(const double) const override;

  std::pair<bool, double>
  checkCollider(const Var, const Var, const Var) const override;

protected:
  TIMER_DECLARE(m_tGrow, mutable);
  TIMER_DECLARE(m_tShrink, mutable);
  TIMER_DECLARE(m_tDist, mutable);
  TIMER_DECLARE(m_tSymmetry, mutable);
  TIMER_DECLARE(m_tSync, mutable);
  TIMER_DECLARE(m_tBlankets, mutable);
  TIMER_DECLARE(m_tNeighbors, mutable);
}; // class BlanketLearning

/**
 * @brief Class that implements Grow-Shrink strategy for MB discovery,
 *        as described by Margaritis & Thrun.
 *
 * @tparam Data Type of the object which is used for querying the data.
 * @tparam Var Type of variable indices (expected to be an integer type).
 * @tparam Set Type of set container.
 */
template <typename Data, typename Var, typename Set>
class GS : public BlanketLearning<Data, Var, Set> {
public:
  GS(const mxx::comm&, const Data&, const Var = std::numeric_limits<Var>::max());

private:
  Set
  getCandidateMB(const Var, Set&&) const override;

  void
  updatePValues(std::vector<std::tuple<Var, Var, double>>&, const std::unordered_map<Var, Set>&) const override;

  std::pair<Var, double>
  pickBestCandidate(const Var, const Set&, const Set&) const override;
}; // class GS

/**
 * @brief Class that implements Incremental Association strategy for MB discovery,
 *        as described by Tsamardinos et al.
 *
 * @tparam Data Type of the object which is used for querying the data.
 * @tparam Var Type of variable indices (expected to be an integer type).
 * @tparam Set Type of set container.
 */
template <typename Data, typename Var, typename Set>
class IAMB : public BlanketLearning<Data, Var, Set> {
public:
  IAMB(const mxx::comm&, const Data&, const Var = std::numeric_limits<Var>::max());

private:
  Set
  getCandidateMB(const Var, Set&&) const override;
}; // class IAMB

/**
 * @brief Class that implements Interleaved Incremental Association strategy
 *        for MB discovery, as described by Tsamardinos et al.
 *
 * @tparam Data Type of the object which is used for querying the data.
 * @tparam Var Type of variable indices (expected to be an integer type).
 * @tparam Set Type of set container.
 */
template <typename Data, typename Var, typename Set>
class InterIAMB : public BlanketLearning<Data, Var, Set> {
public:
  InterIAMB(const mxx::comm&, const Data&, const Var = std::numeric_limits<Var>::max());

private:
  Set
  getCandidateMB(const Var, Set&&) const override;

  void
  growShrink(std::vector<std::tuple<Var, Var, double>>&&, std::unordered_map<Var, Set>&, std::set<std::pair<Var, Var>>&, const double) const override;
}; // class InterIAMB

#include "detail/BlanketLearning.hpp"

#endif // BLANKETLEARNING_HPP_
