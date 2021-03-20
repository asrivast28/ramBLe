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
#ifndef DETAIL_GSQUARE_HPP_
#define DETAIL_GSQUARE_HPP_

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
 * @brief Computes the G^2 statistic using the SABNAtk counters.
 */
template <template <typename...> class CounterType, typename Var, typename Set, typename... CounterArgs>
std::pair<uint32_t, double>
computeGSquareSABNAtk(
  const CounterType<CounterArgs...>& counter,
  const Var x,
  const Var y,
  const Set& given
)
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
computeGSquare(
  const CounterType<CounterArgs...>& counter,
  const Var x,
  const Var y,
  const Set& given
)
{
  return computeGSquareSABNAtk(counter, x, y, given);
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
computeGSquare(
  const CounterType<CounterArgs...>& counter,
  const Var x,
  const Var y,
  const Set& given
)
{
  return computeGSquareSABNAtk(counter, x, y, given);
}

/**
 * @brief Compute marginal G^2 statistic using contingency tables.
 */
template <typename Var>
std::pair<uint32_t, double>
marginalGSquare(
  const CTCounter<>& counter,
  const Var x,
  const Var y
)
{
  auto r_x = counter.r(x);
  auto r_y = counter.r(y);

  uint32_t df = (r_x - 1) * (r_y - 1);
  LOG_MESSAGE(trace, "r_x = %d, r_y = %d", r_x, r_y);

  // Storage for row counts corresponding to every possible configuration
  uint32_t* cc = new uint32_t[r_x * r_y]();
  // Storage for marginal counts for the states of x
  uint32_t* mcx = new uint32_t[r_x]();
  // Storage for marginal counts for the states of y
  uint32_t* mcy = new uint32_t[r_y]();

  const auto nobs = counter.m();
  auto xx = &counter.data()[x * nobs];
  auto yy = &counter.data()[y * nobs];
  for (auto k = 0u; k < nobs; ++k) {
    ++cc[xx[k] * r_y + yy[k]];
  }

  for (auto a = 0u, idx = 0u; a < r_x; ++a) {
    for (auto b = 0u; b < r_y; ++b) {
      mcx[a] += cc[idx];
      mcy[b] += cc[idx];
      ++idx;
    }
  }

  double gSquare = 0.0;
  for (auto a = 0u, idx = 0u; a < r_x; ++a) {
    auto first = static_cast<double>(nobs) / mcx[a];
    for (auto b = 0u; b < r_y; ++b) {
      LOG_MESSAGE(trace, "a = %d, b = %d", static_cast<int>(a), static_cast<int>(b));
      LOG_MESSAGE(trace, "si = %d, sj = %d, sij = %d", mcx[a], mcy[b], cc[idx]);
      // XXX: This is one of the only two places where we multiply three observation counts
      // It can lead to overflow for 32-bit unsigned int
      // Use 64-bit unsigned int just for this for now
      if ((static_cast<uint64_t>(cc[idx]) * mcx[a] * mcy[b] != 0) && (cc[idx] * nobs != mcx[a] * mcy[b])) {
        auto component = cc[idx] * log((first * cc[idx]) / mcy[b]);
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

  delete[] cc;
  delete[] mcx;
  delete[] mcy;

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
)
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
 * @brief Compute conditional G^2 statistic using contingency tables.
 */
template <typename Var, typename Set>
std::pair<uint32_t, double>
conditionalGSquare(
  const CTCounter<>& counter,
  const Var x,
  const Var y,
  const Set& given
)
{
  auto r_x = counter.r(x);
  auto r_y = counter.r(y);
  auto r_xy = r_x * r_y;

  uint32_t df = (r_x - 1) * (r_y - 1);
  LOG_MESSAGE(trace, "r_x = %d, r_y = %d", r_x, r_y);

  // Storage for the computed indices for all the observations
  uint32_t* zz = new uint32_t[counter.m()];
  auto r_given = indexGiven(counter, given, zz);
  df *= r_given;

  // Storage for row counts corresponding to every possible configuration
  uint32_t* cc = new uint32_t[r_given * r_xy]();
  // Storage for marginal counts for the states of x and the given variables
  uint32_t* mcx = new uint32_t[r_given * r_x]();
  // Storage for marginal counts for the states of y and the given variables
  uint32_t* mcy = new uint32_t[r_given * r_y]();
  // Storage for marginal counts for the states of the given variables
  uint32_t* mcz = new uint32_t[r_given]();

  auto xx = &counter.data()[x * counter.m()];
  auto yy = &counter.data()[y * counter.m()];
  for (auto k = 0u; k < counter.m(); ++k) {
    ++cc[zz[k] * r_xy + xx[k] * r_y + yy[k]];
  }

  for (auto c = 0u, idx = 0u; c < r_given; ++c) {
    for (auto a = 0u, i = c * r_x; a < r_x; ++a, ++i) {
      for (auto b = 0u, j = c * r_y; b < r_y; ++b, ++j) {
        mcx[i] += cc[idx];
        mcy[j] += cc[idx];
        mcz[c] += cc[idx];
        ++idx;
      }
    }
  }

  double gSquare = 0.0;
  for (auto c = 0u, idx = 0u; c < r_given; ++c) {
    if (mcz[c] == 0) {
      idx += r_xy;
      continue;
    }
    for (auto a = 0u, i = c * r_x; a < r_x; ++a, ++i) {
      auto first = static_cast<double>(mcz[c]) / mcx[i];
      for (auto b = 0u, j = c * r_y; b < r_y; ++b, ++j) {
        LOG_MESSAGE(trace, "a = %d, b = %d", static_cast<int>(a), static_cast<int>(b));
        LOG_MESSAGE(trace, "sk = %d, sik = %d, sjk = %d, s = %d", mcz[c], mcx[i], mcy[j], cc[idx]);
        // XXX: This is one of the only two places where we multiply three observation counts
        // It can lead to overflow for 32-bit unsigned int
        // Use 64-bit unsigned int just for this for now
        if ((static_cast<uint64_t>(cc[idx]) * mcx[i] * mcy[j] != 0) && (cc[idx] * mcz[c] != mcx[i] * mcy[j])) {
          auto component = cc[idx] * log((first * cc[idx]) / mcy[j]);
          gSquare += component;
          LOG_MESSAGE(trace, "component = %g", component);
        }
        LOG_MESSAGE_IF((cc[idx] * mcx[i] * mcy[j] == 0) || (cc[idx] * mcz[c] == mcx[i] * mcy[j]),
                       trace, "component = 0.0");
        ++idx;
      }
    }
  }
  gSquare *= 2.0;
  LOG_MESSAGE(debug, "df = %d, G-square = %g", df, gSquare);

  delete[] zz;
  delete[] cc;
  delete[] mcx;
  delete[] mcy;
  delete[] mcz;

  return std::make_pair(df, gSquare);
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
template <template <typename...> class CounterType, typename Var, typename Set, typename... CounterArgs> //, typename std::enable_if<sizeof...(CounterArgs) == 1, bool>::type = true>
typename std::enable_if<
  std::is_same<CounterType<CounterArgs...>, CTCounter<CounterArgs...>>::value,
  std::pair<uint32_t, double>
>::type
computeGSquare(
  const CounterType<CounterArgs...>& counter,
  const Var x,
  const Var y,
  const Set& given
)
{
  if (given.size() > 0) {
    LOG_MESSAGE(trace, "Computing conditional G-square");
    return conditionalGSquare(counter, x, y, given);
  }
  else {
    LOG_MESSAGE(trace, "Computing marginal G-square");
    return marginalGSquare(counter, x, y);
  }
}

#endif // DETAIL_GSQUARE_HPP_
