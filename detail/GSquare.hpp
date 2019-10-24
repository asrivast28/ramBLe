/**
 * @file GSquare.hpp
 * @brief Definition of the functions used for computing G^2.
 */
#ifndef DETAIL_GSQUARE_HPP_
#define DETAIL_GSQUARE_HPP_

#include "utils/Logging.hpp"
#include "BVCounter.hpp"
#include "RadCounter.hpp"


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
 * @brief Computes the G^2 statistic for the variables, given the conditioning set,
 *        and the corresponding degree of freedom.
 *
 * @tparam CounterType Type of the counter to be used.
 * @tparam Var Type of the variables (expected to be an integral type).
 * @tparam Set The type of the container used for indices of the given variables.
 * @param x The index of the first variable.
 * @param y The index of the second variable.
 * @param given The indices of the variables to be conditioned on.
 *
 * @return Pair of degree of freedom and the computed G^2 value.
 */
template <template <int, typename...> class CounterType, int N, typename Var, typename Set>
typename std::enable_if<
  std::is_same<CounterType<N>, BVCounter<N>>::value ||
  std::is_same<CounterType<N>, RadCounter<N>>::value,
  std::pair<uint32_t, double>
>::type
computeGSquare(
  const CounterType<N>& counter,
  const Var x,
  const Var y,
  const Set& given
)
{
  using data_type = typename CounterType<N>::data_type;

  auto r_x = counter.r(x);
  auto r_y = counter.r(y);

  uint32_t df = (r_x - 1) * (r_y - 1);
  LOG_MESSAGE(trace, "r_x = %d, r_y = %d", r_x, r_y);
  double gSquare = 0.0;

  std::vector<data_type> r(given.size());
  std::vector<int> pa(given.size());
  auto k = 0u;
  for (auto xk = given.begin(); xk != given.end(); ++xk, ++k) {
    pa[k] = *xk;
    r[k] = counter.r(*xk);
    df *= r[k];
  }

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
        if (sjk == 0) {
          continue;
        }
        auto sijk = counter.count(counter.common(count_x, static_cast<int>(y), b));
        if (sijk == 0) {
          continue;
        }
        LOG_MESSAGE(trace, "a = %d, b = %d", static_cast<int>(a), static_cast<int>(b));
        LOG_MESSAGE(trace, "sk = %d, sik = %d, sjk = %d, sijk = %d", sk, sik, sjk, sijk);
        if (sijk * sk == sik * sjk) {
          LOG_MESSAGE(trace, "component = 0.0");
          continue;
        }
        auto component = sijk * (log(sijk) + log(sk) - log(sik) - log(sjk));
        gSquare += component;
        LOG_MESSAGE(trace, "component = %g", component);
      }
    }
  }
  gSquare *= 2.0;
  LOG_MESSAGE(debug, "df = %d, G-square = %g", df, gSquare);
  return std::make_pair(df, gSquare);
}

#endif // DETAIL_GSQUARE_HPP_
