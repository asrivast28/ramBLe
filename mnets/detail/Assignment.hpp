/**
 * @file Assignment.hpp
 * @brief Implementation of parent and split assignment classes.
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
#ifndef DETAIL_ASSIGNMENTS_HPP_
#define DETAIL_ASSIGNMENTS_HPP_


template <typename Data, typename Var, typename Set>
class TreeNode;

template <typename Data, typename Var, typename Set>
class Assignment {
public:
  Assignment(
    const Data& data,
    const Var v,
    const TreeNode<Data, Var, Set>* const node
  ) : m_data(),
      m_sum(std::make_pair(0.0, 0.0)),
      m_missing(std::make_pair(0, 0)),
      m_varName(data.varName(v))
  {
    const auto& firstObs = node->children().first->observations();
    m_data.first = std::vector<double>(firstObs.size());
    auto d = m_data.first.begin();
    for (const auto o : firstObs) {
      if (!std::isnan(data(v, o))) {
        *d = data(v, o);
        m_sum.first += *d;
        ++d;
      }
    }
    m_missing.first = std::distance(d, m_data.first.end());
    m_data.first.resize(std::distance(m_data.first.begin(), d));
    const auto& secondObs = node->children().second->observations();
    m_data.second = std::vector<double>(secondObs.size());
    d = m_data.second.begin();
    for (const auto o : secondObs) {
      if (!std::isnan(data(v, o))) {
        *d = data(v, o);
        m_sum.second += *d;
        ++d;
      }
    }
    m_missing.second = std::distance(d, m_data.second.end());
    m_data.second.resize(std::distance(m_data.second.begin(), d));
  }

  const std::string&
  varName() const
  {
    return m_varName;
  }

  int
  sign(const double sv) const
  {
    int signFirst = -1, signSecond = -1;
    if (std::isless(m_sum.first, sv * m_data.first.size())) {
      signFirst = 1;
    }
    if (std::isgreater(m_sum.second, sv * m_data.second.size())) {
      signSecond = 1;
    }
    return (signFirst == signSecond) ? signFirst : 0;
  }

  double
  evaluate(
    const double sv,
    const int sign,
    const double beta
  ) const
  {
    auto sumFirst = 0.0;
    for (const auto x : m_data.first) {
      sumFirst += (x - sv) / (1 + exp(-1 * sign * beta * (x - sv)));
    }
    auto sumSecond = 0.0;
    for (const auto x : m_data.second) {
      sumSecond += (x - sv) / (1 + exp(sign * beta * (x - sv)));
    }
    return sumSecond - sumFirst;
  }

  double
  score(
    const double sv,
    const int sign,
    const double beta
  ) const
  {
    auto sumFirst = 0.0;
    for (const auto x : m_data.first) {
      sumFirst -= log(1 + exp(sign * beta * (x - sv)));
    }
    sumFirst -= (m_missing.first * log(2));
    auto sumSecond = 0.0;
    for (const auto x : m_data.second) {
      sumSecond -= log(1 + exp(-1 * sign * beta * (x - sv)));
    }
    sumSecond -= (m_missing.second * log(2));
    return sumFirst + sumSecond;
  }

private:
  std::pair<std::vector<double>, std::vector<double>> m_data;
  std::pair<double, double> m_sum;
  std::pair<uint32_t, uint32_t> m_missing;
  std::string m_varName;
}; // class Assignment<Data, Var, Set>

#endif // DETAIL_ASSIGNMENT_HPP_
