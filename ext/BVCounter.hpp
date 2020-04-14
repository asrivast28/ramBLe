/***
 *    $Id$
 **
 *    File: BVCounter.hpp
 *    Created: Oct 22, 2016
 *
 *    Authors: Matthew Eichhorn <maeichho@buffalo.edu>
 *             Blake Hurlburt <blakehur@buffalo.edu>
 *             Grant Iraci <grantira@buffalo.edu>
 *    Copyright (c) 2015-2017 SCoRe Group http://www.score-group.org/
 *    Distributed under the MIT License.
 *    See accompanying file LICENSE.
 */

#ifndef BV_COUNTER_HPP
#define BV_COUNTER_HPP

#include <algorithm>
#include <cmath>
#include <vector>

#include <bitvector.hpp>
#include <bit_util.hpp>


template <typename Size> class BVCounter {
public:
    typedef uint_type<Size::value> set_type;
    typedef uint8_t data_type;

    int n() const { return n_; }

    int m() const { return m_; }

    int r(int i) const { return data_[i].first.size(); }

    bool is_reorderable() { return true; }

    int count(const bitvector& res) const { return res.weight(); }

    template <typename data_type>
    bitvector common(const std::vector<int>& pa, const std::vector<data_type>& state_pa) const {
        vect res = vect::identity(m_);
        vect temp_res;
        for (auto i = 0u; i < state_pa.size(); ++i) {
            int id = r_idx_[pa[i]][state_pa[i]];
            intersect(res, data_[pa[i]].first[id], temp_res);
            if (temp_res.weight() == 0) return temp_res;
            res = temp_res;
        }
        return res;
    }

    template <typename data_type>
    bitvector common(const bitvector& prev_res, const int xi, const data_type state_xi) const {
        vect res;
        int r_id = r_idx_[xi][state_xi];
        intersect(data_[xi].first[r_id], prev_res, res);
        return res;
    }

    template <typename score_functor, typename data_type>
    void apply(const set_type& xi, const set_type& pa, const std::vector<data_type>& state_xi, const std::vector<data_type>& state_pa, std::vector<score_functor>& F) const {
        std::vector<int> xi_vect = as_vector(xi);

        int qpa = m_q__(pa);
        for (int i = 0; i < F.size(); ++i) F[i].init(r(xi_vect[i]), qpa);

        if (is_emptyset(pa)) {
            for (int i = 0; i < F.size(); ++i) {
                int r_id = r_idx_[xi_vect[i]][state_xi[i]];
                F[i](m_);
                F[i](data_[xi_vect[i]].first[r_id].weight(), m_);
                F[i].finalize(1);
            }
            return;
        }

        auto pa_vect = as_vector(pa);

        vect res = vect::identity(m_);
        vect temp_res;

        for (int i = 0; i < state_pa.size(); ++i) {
            int id = r_idx_[pa_vect[i]][state_pa[i]];
            intersect(res, data_[pa_vect[i]].first[id], temp_res);
            if (temp_res.weight() == 0) return;
            res = temp_res;
        }

        for (int i = 0; i < xi_vect.size(); ++i) {
            F[i](res.weight());
            int r_id = r_idx_[xi_vect[i]][state_xi[i]];
            intersect(data_[xi_vect[i]].first[r_id], res, temp_res);
            if (temp_res.weight() != 0) F[i](temp_res.weight(), res.weight());
        }
    } // apply (state specific queries)


    template <typename score_functor>
    void apply(const std::vector<int>& xi_vect, const set_type& pa, std::vector<score_functor>& F) const {
        int qpa = m_q__(pa);
        for (int i = 0; i < F.size(); ++i) F[i].init(r(xi_vect[i]), qpa);

        if (is_emptyset(pa)) {
            for (int i = 0; i < F.size(); ++i) {
                F[i](m_);
                for (const vect& v : data_[xi_vect[i]].first) {
                    F[i](v.weight(), m_);
                }
                F[i].finalize(1);
            }
            return;
        }

        using Iter = typename std::vector<vect>::const_iterator;

        std::vector<int> pa_sorted;

        for (int i : sorted_order_) {
            if (in_set(pa, i)) {
                pa_sorted.push_back(i);
            }
        }

        //the table stores pairs with
        //first = current bitvector at that level of the dfs traversal
        //second = iterator over bitvectors of next parent to intersect with
        std::vector<std::pair<vect, Iter>> table;
        table.reserve(pa_sorted.size());

        for (int i : pa_sorted) {
            table.push_back(std::make_pair(vect::identity(m_), std::begin(data_[i].first))); //start with identity
        }

        //[0] should be ?
        table.push_back(std::make_pair(vect::identity(m_), std::begin(data_[0].first)));
        table.push_back(std::make_pair(vect::identity(m_), std::begin(data_[0].first)));

        int layer = 0; //we start at the top of the tree
        int qi_obs = 0;

        vect container = vect::identity(m_); //bottom layer

        while (true) { //once we are finished, we will try to move one layer up from root

            while (layer >= 0 && table[layer].second == std::end(data_[pa_sorted[layer]].first)) {
                //reset the iterator for next traversal at this level
                table[layer].second = std::begin(data_[pa_sorted[layer]].first);
                --layer; //move up to continue traversal
            }

            if (layer < 0) break;

            const vect& root = table[layer].first; //fetch from the table

            intersect(root, *(table[layer].second), table[layer + 1].first);
            ++table[layer].second;

            int Nij = table[layer + 1].first.weight();

            if (Nij) { //if Nijk is zero we just move on
                if (layer == (pa_sorted.size() - 1)) { //we have explored all of the parents
                    for (int i = 0; i < F.size(); ++i) F[i](Nij);
                    ++qi_obs;
                    for (int i = 0; i < F.size(); ++i) {
                        for (const vect& v : data_[xi_vect[i]].first) {
                            intersect(table[layer + 1].first, v, container);
                            int Nijk = container.weight();
                            if (Nijk) {
                                F[i](Nijk, Nij);
                            }
                        }
                    }
                } else if (Nij == 1) {
                    //if we go down a level, we will find that this is the only call that would be made
                    for (int i = 0; i < F.size(); ++i) F[i](Nij);
                    ++qi_obs;
                    for (int i = 0; i < F.size(); ++i) {
                        F[i](1,1);
                    }
                } else { //there are more parents
                    ++layer; //move down a layer and keep going
                }
            }
        }

        for (int i = 0; i < F.size(); ++i) F[i].finalize(qi_obs);
    } // apply

    template <typename score_functor>
    void apply(const set_type& xi, const set_type& pa, std::vector<score_functor>& F) const {
        std::vector<int> xi_vect = as_vector(xi);
        apply(xi_vect, pa, F);
    } // apply

    template <typename score_functor>
    void apply(const std::vector<int>& xi_vect, const std::vector<int>& pa_vect, std::vector<score_functor>& F) const {
        set_type pa = as_set<set_type>(std::begin(pa_vect), std::end(pa_vect));
        apply(xi_vect, pa, F);
    } // apply

    template <typename score_functor>
    void apply(int xi, const set_type& pa, score_functor& F) const {
        std::vector<int> xi_vect{xi};
        std::vector<score_functor> F_vect{F};
        apply(xi_vect, pa, F_vect);
        F = F_vect[0];
    } // apply


    // reorder variables to improve expected query performance
    bool reorder(const std::vector<int>& order) {
        std::vector<std::vector<int>> temp_r_idx;
        std::vector<std::pair<std::vector<vect>, double>> vec;
        vec.reserve(n_);

        for (int i : order) {
            temp_r_idx.emplace_back(std::move(r_idx_[i]));
            vec.emplace_back(std::move(data_[i]));
        }

        r_idx_ = std::move(temp_r_idx);
        data_ = std::move(vec);

        return true;
    } // reorder

    // declaration of static method to create the counter
    template <typename Iter> static BVCounter<Size> create(int n, int m, Iter it);

private:
    typedef bitvector vect;

    class ent_order {
    public:
        ent_order(const std::vector<std::pair<std::vector<vect>, double>>& data) : data_(data) { }

        bool operator()(int lhs, int rhs) { return data_[lhs].second < data_[rhs].second; }

    private:
        const std::vector<std::pair<std::vector<vect>, double>>& data_;
    };

    // assume that q will not overflow, this will be checked by sabna calling code
    int m_q__(const set_type& pa) const {
        int q = 1;

        for (int i = 0; i < set_max_size<set_type>(); ++i) {
            if (in_set(pa, i)) {
                q *= r(i);
            }
        }

        return q;
    } // m_q__

    std::vector<std::vector<int>> r_idx_;
    std::vector<std::pair<std::vector<vect>, double>> data_;
    std::vector<int> sorted_order_;
    int n_ = -1;
    int m_ = -1;
}; // class BVCounter

template <typename Size>
template <typename Iter> BVCounter<Size> BVCounter<Size>::create(int n, int m, Iter it) {
    BVCounter<Size> p;

    int indices[256];
    int temp;
    int size;

    p.n_ = n;
    p.m_ = m;
    p.r_idx_.resize(n);

    for (int xi = 0; n > 0; --n, ++xi) {
        p.data_.push_back(std::make_pair(std::vector<typename BVCounter<Size>::vect>(), 0.0));

        size = 0;
        std::fill_n(indices, 256, -1);

        std::vector<std::pair<int, int>> temp_r;

        for (int j = 0; j < m; ++j) {
            temp = *it++;

            if (indices[temp] == -1) {
                temp_r.push_back({temp, size});
                indices[temp] = size++;
                p.data_.back().first.push_back(typename BVCounter<Size>::vect(m));
            }

            p.data_.back().first[indices[temp]].insert(j);
        }

        std::sort(temp_r.begin(), temp_r.end(), [] (const std::pair<int, int>& lhs, const std::pair<int, int>& rhs) { return lhs.first < rhs.first; });
        for (const auto x : temp_r) p.r_idx_[xi].push_back(x.second);

        double H = 0.0;
        for (auto& v : p.data_.back().first) {
            double px = (static_cast<double>(v.weight()) / m);
            H += px * (px == 0.0 ? 0 : std::log2(px));
        }

        p.data_.back().second = -H;
    }

    for (auto i = 0u; i < p.data_.size(); ++i) {
        p.sorted_order_.push_back(static_cast<int>(i));
    }

    std::sort(std::begin(p.sorted_order_), std::end(p.sorted_order_), typename BVCounter<Size>::ent_order(p.data_));

    return p;
} // BVCounter<Size>::create

#endif // BV_COUNTER_HPP
