/**
 * @file GSquare.hpp
 * @brief Definition of the functions used for computing G^2.
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
#ifndef GSQUARE_HPP_
#define GSQUARE_HPP_

#include "common/CTCounter.hpp"
#include "common/ext/BVCounter.hpp"
#include "common/ext/RadCounter.hpp"
#include "utils/Logging.hpp"


/**
 * @brief Class that provides an iterator over all combinations of
 *        the states that different variables can take.
 *
 * @tparam State Type of the variable states, expected to be an integral type.
 *
 * This class provides a (supposedly) lightweight way of iterating over
 * all the combinations of states of different variables. Variable i
 * is supposed to take all the states between 0 and bounds[i] - 1.
 *
 * It further allows including a single user-provided state for another
 * variable x, whose index is provided by idxX.
 */
template <typename State>
class StateIterator {
public:
  StateIterator(const std::vector<State>& bounds)
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

  const std::vector<State>&
  state() const
  {
    return m_state;
  }

private:
  const std::vector<State> m_bounds;
  std::vector<State> m_state;
  bool m_valid;
}; // class StateIterator


/**
 * @brief Class that implements the G^2 computations
 *        using the different counters.
 */
class GSquare {
public:
  /**
   * @brief Default constructor.
   */
  GSquare()
    : m_zz(nullptr),
      m_cc(nullptr),
      m_cx(nullptr),
      m_cy(nullptr),
      m_cz(nullptr),
      m_zzSize(),
      m_ccSize(),
      m_cxSize(),
      m_cySize(),
      m_czSize()
  {
  }

  /**
   * @brief Default destructor that deallocates all the buffers.
   */
  ~GSquare()
  {
    free(m_zz);
    free(m_cc);
    free(m_cx);
    free(m_cy);
    free(m_cz);
  }

  /**
   * @brief Computes the G^2 statistic for the variables, given the conditioning set,
   *        and the corresponding degree of freedom using radix counters.
   *
   * @tparam CounterType Type of the counter to be used.
   * @tparam Var Type of the variables (expected to be an integral type).
   * @tparam Set The type of the container used for indices of the given variables.
   * @tparam CounterArgs Type of the arguments list for the counter type.
   * @param x The index of the first variable.
   * @param y The index of the second variable.
   * @param given The indices of the variables to be conditioned on.
   *
   * @return Pair of degree of freedom and the computed G^2 value.
   */
  template <template <typename...> class CounterType, typename Var, typename Set, typename... CounterArgs>
  typename std::enable_if<
    std::is_same<CounterType<CounterArgs...>, RadCounter<CounterArgs...>>::value,
    std::pair<uint32_t, double>
  >::type
  compute(
    const CounterType<CounterArgs...>& counter,
    const Var x,
    const Var y,
    const Set& given
  ) const
  {
    return computeSABNAtk(counter, x, y, given);
  }

  /**
   * @brief Computes the G^2 statistic for the variables, given the conditioning set,
   *        and the corresponding degree of freedom using bitvector counters.
   *
   * @tparam CounterType Type of the counter to be used.
   * @tparam Var Type of the variables (expected to be an integral type).
   * @tparam Set The type of the container used for indices of the given variables.
   * @tparam CounterArgs Type of the arguments list for the counter type.
   * @param x The index of the first variable.
   * @param y The index of the second variable.
   * @param given The indices of the variables to be conditioned on.
   *
   * @return Pair of degree of freedom and the computed G^2 value.
   */
  template <template <typename...> class CounterType, typename Var, typename Set, typename... CounterArgs>
  typename std::enable_if<
    std::is_same<CounterType<CounterArgs...>, BVCounter<CounterArgs...>>::value,
    std::pair<uint32_t, double>
  >::type
  compute(
    const CounterType<CounterArgs...>& counter,
    const Var x,
    const Var y,
    const Set& given
  ) const
  {
    return computeSABNAtk(counter, x, y, given);
  }

  /**
   * @brief Computes the G^2 statistic for the variables, given the conditioning set,
   *        and the corresponding degree of freedom using contingency tables.
   *
   * @tparam CounterType Type of the counter to be used.
   * @tparam Var Type of the variables (expected to be an integral type).
   * @tparam Set The type of the container used for indices of the given variables.
   * @tparam CounterArgs Type of the arguments list for the counter type.
   * @param x The index of the first variable.
   * @param y The index of the second variable.
   * @param given The indices of the variables to be conditioned on.
   *
   * @return Pair of degree of freedom and the computed G^2 value.
   */
  template <template <typename...> class CounterType, typename Var, typename Set, typename... CounterArgs>
  typename std::enable_if<
    std::is_same<CounterType<CounterArgs...>, CTCounter<CounterArgs...>>::value,
    std::pair<uint32_t, double>
  >::type
  compute(
    const CounterType<CounterArgs...>& counter,
    const Var x,
    const Var y,
    const Set& given
  ) const
  {
    if (given.size() > 0) {
      LOG_MESSAGE(trace, "Computing conditional G-square");
      return conditional(counter, x, y, given);
    }
    else {
      LOG_MESSAGE(trace, "Computing marginal G-square");
      return marginal(counter, x, y);
    }
  }

private:
  /**
   * @brief Computes the G^2 statistic using the SABNAtk counters.
   */
  template <template <typename...> class CounterType, typename Var, typename Set, typename... CounterArgs>
  std::pair<uint32_t, double>
  computeSABNAtk(
    const CounterType<CounterArgs...>& counter,
    const Var x,
    const Var y,
    const Set& given
  ) const
  {
    using data_type = typename CounterType<CounterArgs...>::data_type;

    auto r_x = counter.r(x);
    auto r_y = counter.r(y);

    uint32_t df = (r_x - 1) * (r_y - 1);
    LOG_MESSAGE(trace, "r_x = %d, r_y = %d", r_x, r_y);

    std::vector<data_type> r(given.size());
    std::vector<int> pa(given.size());
    auto k = 0u;
    for (auto xk = given.begin(); xk != given.end(); ++xk, ++k) {
      pa[k] = *xk;
      r[k] = counter.r(*xk);
      df *= r[k];
    }

    double gSquare = 0.0;
    for (auto c = StateIterator<data_type>(r); c.valid(); c.next()) {
      auto base = counter.common(pa, c.state());
      auto sk = counter.count(base);
      if (sk == 0) {
        continue;
      }
      for (data_type a = 0; a < r_x; ++a) {
        auto count_x = counter.common(base, static_cast<int>(x), a);
        auto sik = counter.count(count_x);
        if (sik == 0) {
          continue;
        }
        for (data_type b = 0; b < r_y; ++b) {
          auto sjk = counter.count(counter.common(base, static_cast<int>(y), b));
          auto s = counter.count(counter.common(count_x, static_cast<int>(y), b));
          if (s * sjk != 0) {
            LOG_MESSAGE(trace, "a = %d, b = %d", static_cast<int>(a), static_cast<int>(b));
            LOG_MESSAGE(trace, "sk = %d, sik = %d, sjk = %d, s = %d", sk, sik, sjk, s);
            if (s * sk != sik * sjk) {
              auto component = s * (log(s) + log(sk) - log(sik) - log(sjk));
              gSquare += component;
              LOG_MESSAGE(trace, "component = %g", component);
            }
            else {
              LOG_MESSAGE(trace, "component = 0.0");
            }
          }
        }
      }
    }
    gSquare *= 2.0;
    LOG_MESSAGE(debug, "df = %d, G-square = %g", df, gSquare);
    return std::make_pair(df, gSquare);
  }

  /**
   * @brief Computes the configurations of the conditioning set.
   */
  template <typename Set>
  uint32_t
  indexGiven(
    const CTCounter<>& counter,
    const Set& given,
    uint32_t* const indices
  ) const
  {
    const auto& data = counter.data();
    const auto nobs = counter.m();
    auto xk = given.begin();
    std::copy(&data[*xk * nobs], &data[(*xk + 1) * nobs], indices);
    uint32_t cumulative = counter.r(*xk);
    ++xk;
    for (; xk != given.end(); ++xk) {
      const auto* datak = &data[*xk * nobs];
      for (auto i = 0u; i < nobs; ++i) {
        indices[i] += (datak[i] * cumulative);
      }
      cumulative *= counter.r(*xk);
    }
    return cumulative;
  }

  /**
   * @brief Function for efficiently allocating buffers.
   *
   * @param buffer Pointer to the previously allocated buffer.
   * @param buffSize Size of the previously allocated buffer.
   * @param reqdSize Required size of the buffer.
   * @param setZero If the required buffer needs to be zeroed out.
   *
   * @return Pointer to the allocated buffer of the required size.
   */
  uint32_t*
  allocate(
    uint32_t* buffer,
    uint32_t& buffSize,
    const uint32_t reqdSize,
    const bool setZero = false
  ) const
  {
    if (buffSize == 0) {
      LOG_MESSAGE(debug, "Allocating a buffer of size %u", reqdSize);
      buffSize = reqdSize;
      buffer = static_cast<uint32_t*>(malloc(sizeof(uint32_t) * buffSize));
      if (buffer == nullptr) {
        throw std::runtime_error("Failed to allocate a buffer of size " + std::to_string(buffSize));
      }
      LOG_MESSAGE(debug, "Allocated the buffer");
    }
    else if (buffSize < reqdSize) {
      LOG_MESSAGE(debug, "Reallocating a buffer of size %u to %u", buffSize, reqdSize);
      buffSize = reqdSize;
      if (buffer == nullptr) {
        throw std::runtime_error("Reallocating a nullptr");
      }
      buffer = static_cast<uint32_t*>(realloc(buffer, sizeof(uint32_t) * buffSize));
      if (buffer == nullptr) {
        throw std::runtime_error("Failed to reallocate a buffer of size " + std::to_string(buffSize));
      }
      LOG_MESSAGE(debug, "Reallocated the buffer");
    }
    if (setZero) {
      memset(static_cast<void*>(buffer), 0, sizeof(uint32_t) * reqdSize);
    }
    return buffer;
  }

  /**
   * @brief Compute conditional G^2 statistic using contingency tables.
   */
  template <typename Var, typename Set>
  std::pair<uint32_t, double>
  conditional(
    const CTCounter<>& counter,
    const Var x,
    const Var y,
    const Set& given
  ) const
  {
    auto r_x = counter.r(x);
    auto r_y = counter.r(y);
    auto r_xy = r_x * r_y;
    const auto nobs = counter.m();

    uint32_t df = (r_x - 1) * (r_y - 1);
    LOG_MESSAGE(trace, "r_x = %d, r_y = %d", r_x, r_y);

    // Storage for the computed indices for all the observations
    m_zz = allocate(m_zz, m_zzSize, nobs);
    auto r_given = indexGiven(counter, given, m_zz);
    df *= r_given;

    // Storage for row counts corresponding to every possible configuration
    m_cc = allocate(m_cc, m_ccSize, r_given * r_xy, true);
    // Storage for marginal counts for the states of x and the given variables
    m_cx = allocate(m_cx, m_cxSize, r_given * r_x, true);
    // Storage for marginal counts for the states of y and the given variables
    m_cy = allocate(m_cy, m_cySize, r_given * r_y, true);
    // Storage for marginal counts for the states of the given variables
    m_cz = allocate(m_cz, m_czSize, r_given, true);

    auto xx = &counter.data()[x * nobs];
    auto yy = &counter.data()[y * nobs];
    for (auto k = 0u; k < nobs; ++k) {
      ++m_cc[m_zz[k] * r_xy + xx[k] * r_y + yy[k]];
    }

    for (auto c = 0u, idx = 0u; c < r_given; ++c) {
      for (auto a = 0u, i = c * r_x; a < r_x; ++a, ++i) {
        for (auto b = 0u, j = c * r_y; b < r_y; ++b, ++j) {
          m_cx[i] += m_cc[idx];
          m_cy[j] += m_cc[idx];
          m_cz[c] += m_cc[idx];
          ++idx;
        }
      }
    }

    double gSquare = 0.0;
    for (auto c = 0u, idx = 0u; c < r_given; ++c) {
      if (m_cz[c] == 0) {
        idx += r_xy;
        continue;
      }
      for (auto a = 0u, i = c * r_x; a < r_x; ++a, ++i) {
        auto first = static_cast<double>(m_cz[c]) / m_cx[i];
        for (auto b = 0u, j = c * r_y; b < r_y; ++b, ++j) {
          LOG_MESSAGE(trace, "a = %d, b = %d", static_cast<int>(a), static_cast<int>(b));
          LOG_MESSAGE(trace, "sk = %d, sik = %d, sjk = %d, s = %d", m_cz[c], m_cx[i], m_cy[j], m_cc[idx]);
          // XXX: This is one of the only two places where we multiply three observation counts
          // It can lead to overflow for 32-bit unsigned int
          // Use 64-bit unsigned int just for this for now
          if ((static_cast<uint64_t>(m_cc[idx]) * m_cx[i] * m_cy[j] != 0) && (m_cc[idx] * m_cz[c] != m_cx[i] * m_cy[j])) {
            auto component = m_cc[idx] * log((first * m_cc[idx]) / m_cy[j]);
            gSquare += component;
            LOG_MESSAGE(trace, "component = %g", component);
          }
          LOG_MESSAGE_IF((m_cc[idx] * m_cx[i] * m_cy[j] == 0) || (m_cc[idx] * m_cz[c] == m_cx[i] * m_cy[j]),
                         trace, "component = 0.0");
          ++idx;
        }
      }
    }
    gSquare *= 2.0;
    LOG_MESSAGE(debug, "df = %d, G-square = %g", df, gSquare);

    return std::make_pair(df, gSquare);
  }

  /**
   * @brief Compute marginal G^2 statistic using contingency tables.
   */
  template <typename Var>
  std::pair<uint32_t, double>
  marginal(
    const CTCounter<>& counter,
    const Var x,
    const Var y
  ) const
  {
    auto r_x = counter.r(x);
    auto r_y = counter.r(y);

    uint32_t df = (r_x - 1) * (r_y - 1);
    LOG_MESSAGE(trace, "r_x = %d, r_y = %d", r_x, r_y);

    // Storage for row counts corresponding to every possible configuration
    m_cc = allocate(m_cc, m_ccSize, r_x * r_y, true);
    // Storage for marginal counts for the states of x
    m_cx = allocate(m_cx, m_cxSize, r_x, true);
    // Storage for marginal counts for the states of y
    m_cy = allocate(m_cy, m_cySize, r_y, true);

    const auto nobs = counter.m();
    auto xx = &counter.data()[x * nobs];
    auto yy = &counter.data()[y * nobs];
    for (auto k = 0u; k < nobs; ++k) {
      ++m_cc[xx[k] * r_y + yy[k]];
    }

    for (auto a = 0u, idx = 0u; a < r_x; ++a) {
      for (auto b = 0u; b < r_y; ++b) {
        m_cx[a] += m_cc[idx];
        m_cy[b] += m_cc[idx];
        ++idx;
      }
    }

    double gSquare = 0.0;
    for (auto a = 0u, idx = 0u; a < r_x; ++a) {
      auto first = static_cast<double>(nobs) / m_cx[a];
      for (auto b = 0u; b < r_y; ++b) {
        LOG_MESSAGE(trace, "a = %d, b = %d", static_cast<int>(a), static_cast<int>(b));
        LOG_MESSAGE(trace, "si = %d, sj = %d, sij = %d", m_cx[a], m_cy[b], m_cc[idx]);
        // XXX: This is one of the only two places where we multiply three observation counts
        // It can lead to overflow for 32-bit unsigned int
        // Use 64-bit unsigned int just for this for now
        if ((static_cast<uint64_t>(m_cc[idx]) * m_cx[a] * m_cy[b] != 0) && (m_cc[idx] * nobs != m_cx[a] * m_cy[b])) {
          auto component = m_cc[idx] * log((first * m_cc[idx]) / m_cy[b]);
          gSquare += component;
          LOG_MESSAGE(trace, "component = %g", component);
        }
        else {
          LOG_MESSAGE(trace, "component = 0.0");
        }
        ++idx;
      }
    }
    gSquare *= 2.0;
    LOG_MESSAGE(debug, "df = %d, G-square = %g", df, gSquare);

    return std::make_pair(df, gSquare);
  }

private:
  mutable uint32_t* m_zz;
  mutable uint32_t* m_cc;
  mutable uint32_t* m_cx;
  mutable uint32_t* m_cy;
  mutable uint32_t* m_cz;
  mutable uint32_t m_zzSize;
  mutable uint32_t m_ccSize;
  mutable uint32_t m_cxSize;
  mutable uint32_t m_cySize;
  mutable uint32_t m_czSize;
}; // class GSquare

#endif // GSQUARE_HPP_
