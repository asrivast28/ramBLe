/**
 * @file DirectLearning.hpp
 * @brief Declaration of the classes for direct learning algorithms.
 */
#ifndef DIRECTLEARNING_HPP_
#define DIRECTLEARNING_HPP_

#include "ConstraintBasedLearning.hpp"


/**
 * @brief Abstract base class for causal discovery by first learning PC sets.
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
  ~DirectLearning() { }
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
