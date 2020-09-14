/**
 * @file DirectLearning.hpp
 * @brief Declaration of the classes for direct learning algorithms.
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
#ifndef DIRECTLEARNING_HPP_
#define DIRECTLEARNING_HPP_

#include "ConstraintBasedLearning.hpp"


/**
 * @brief Abstract base class for causal discovery by directly learning PC sets.
 *
 * @tparam Data Type of the object which is used for querying the data.
 * @tparam Var Type of variable indices (expected to be an integer type).
 * @tparam Set Type of set container.
 */
template <typename Data, typename Var, typename Set>
class DirectLearning : public ConstraintBasedLearning<Data, Var, Set> {
public:
  DirectLearning(const mxx::comm&, const Data&, const Var);

  Set
  removeFalsePC(const Var, Set&) const;

  Set
  getCandidateMB(const Var, Set&&) const override;

  virtual
  ~DirectLearning();

protected:
  void
  updateMaxPValues(const Var, std::vector<std::pair<double, Var>>&, const Set&, const Set&) const;

  void
  updateMyPValues(std::vector<std::tuple<Var, Var, double>>&, const std::unordered_map<Var, Set>&, const std::unordered_map<Var, Set>&) const;

  template <typename Compare>
  std::set<std::tuple<Var, Var, double>>
  forwardPhase(const std::vector<std::tuple<Var, Var, double>>&, const Compare&, const bool, std::unordered_map<Var, Set>&) const;

  std::set<std::pair<Var, Var>>
  backwardPhase(std::unordered_map<Var, Set>&) const;

  virtual
  void
  forwardBackward(std::vector<std::tuple<Var, Var, double>>&&, std::unordered_map<Var, Set>&, std::set<std::pair<Var, Var>>&, const double) const { };

private:
  BayesianNetwork<Var>
  getSkeleton_parallel(const double, std::unordered_map<Var, Set>&, std::unordered_map<Var, Set>&) const override;

protected:
  TIMER_DECLARE(m_tForward, mutable);
  TIMER_DECLARE(m_tBackward, mutable);
  TIMER_DECLARE(m_tDist, mutable);
  TIMER_DECLARE(m_tSymmetry, mutable);
  TIMER_DECLARE(m_tSync, mutable);
  TIMER_DECLARE(m_tNeighbors, mutable);
}; // class DirectLearning


/**
 * @brief Class that implements Max-Min algorithm for PC discovery,
 *        as described by Tsamardinos et al.
 *
 * @tparam Data Type of the object which is used for querying the data.
 * @tparam Var Type of variable indices (expected to be an integer type).
 * @tparam Set Type of set container.
 */
template <typename Data, typename Var, typename Set>
class MMPC : public DirectLearning<Data, Var, Set> {
public:
  MMPC(const mxx::comm&, const Data&, const Var = std::numeric_limits<Var>::max());

  Set
  getCandidatePC(const Var, Set&&) const override;

private:
  void
  forwardBackward(std::vector<std::tuple<Var, Var, double>>&&, std::unordered_map<Var, Set>&, std::set<std::pair<Var, Var>>&, const double) const override;
}; // class MMPC

/**
 * @brief Class that implements HITON algorithm for PC discovery,
 *        as described by Aliferis et al.
 *
 * @tparam Data Type of the object which is used for querying the data.
 * @tparam Var Type of variable indices (expected to be an integer type).
 * @tparam Set Type of set container.
 */
template <typename Data, typename Var, typename Set>
class HITON : public DirectLearning<Data, Var, Set> {
public:
  HITON(const mxx::comm&, const Data&, const Var = std::numeric_limits<Var>::max());

private:
  Set
  getCandidatePC(const Var, Set&&) const override;
}; // class HITON

/**
 * @brief Class that implements Semi-interleaved HITON algorithm for PC discovery,
 *        as described by Aliferis et al.
 *
 * @tparam Data Type of the object which is used for querying the data.
 * @tparam Var Type of variable indices (expected to be an integer type).
 * @tparam Set Type of set container.
 */
template <typename Data, typename Var, typename Set>
class SemiInterleavedHITON : public DirectLearning<Data, Var, Set> {
public:
  SemiInterleavedHITON(const mxx::comm&, const Data&, const Var = std::numeric_limits<Var>::max());

private:
  Set
  getCandidatePC(const Var, Set&&) const override;

  void
  forwardBackward(std::vector<std::tuple<Var, Var, double>>&&, std::unordered_map<Var, Set>&, std::set<std::pair<Var, Var>>&, const double) const override;
}; // class SemiInterleavedHITON

/**
 * @brief Class that implements GetPC algorithm for PC discovery,
 *        as described by Pena et al.
 *
 * @tparam Data Type of the object which is used for querying the data.
 * @tparam Var Type of variable indices (expected to be an integer type).
 * @tparam Set Type of set container.
 */
template <typename Data, typename Var, typename Set>
class GetPC : public DirectLearning<Data, Var, Set> {
public:
  GetPC(const mxx::comm&, const Data&, const Var = std::numeric_limits<Var>::max());

private:
  Set
  getCandidatePC(const Var, Set&&) const override;
}; // class GetPC

#include "detail/DirectLearning.hpp"

#endif // DIRECTLEARNING_HPP_
