/*********************************************************************
Author: Soonho Kong <soonhok@cs.cmu.edu>

dReal -- Copyright (C) 2013 - 2016, the dReal Team

dReal is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

dReal is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with dReal. If not, see <http://www.gnu.org/licenses/>.
*********************************************************************/

#pragma once

#include <iostream>
#include <unordered_set>
#include <vector>
#include <memory>
#include "constraint/constraint.h"
#include "ibex/ibex.h"
#include "opensmt/egraph/Enode.h"
#include "util/box.h"
#include "util/hash_combine.h"
#include "util/ibex_interval_hash.h"

namespace dreal {

// Generalized Box
class gbox {
private:
    std::vector<Enode *> m_vars;
    ibex::IntervalVector m_vec;
    std::vector<bool> m_gbits;

public:
    explicit gbox(box const & b);
    explicit gbox(std::vector<Enode *> const & vars);
    gbox(std::vector<Enode *> const & vars, ibex::IntervalVector const & vec);
    gbox(std::vector<Enode *> const & vars, ibex::IntervalVector const & vec, std::vector<bool> const & gbits);
    gbox(box const & b, std::unordered_set<std::shared_ptr<constraint>> const & used_ctrs, ibex::IntervalVector const & dom);
    unsigned size() const { return m_vec.size(); }
    bool operator==(gbox const & gb) const;
    bool operator!=(gbox const & gb) const { return !(*this == gb); }
    bool is_empty() const { return size() == 0 || m_vec.is_empty(); }
    bool is_subset(gbox const & gb) const;
    bool is_superset(gbox const & gb) const { return gb.is_subset(*this); }
    bool match(gbox const & ant, gbox const & con) const;
    bool is_strict_overlap(gbox const & gb) const {
        if (is_empty() || gb.is_empty()) {
            return false;
        }
        for (unsigned i = 0; i < size(); ++i) {
            if (!get_gbit(i) && !gb.get_gbit(i)) {
                // ibex::Interval::overlaps is checking strict_overlap
                if (!m_vec[i].overlaps(gb.m_vec[i])) {
                    return false;
                }
            }
            assert(m_vec[i].overlaps(gb.m_vec[i]));
        }
        return true;
    }
    bool generalized() const {
        return any_of(m_gbits.begin(), m_gbits.end(), [](bool const gbit) { return gbit; });
    }
    void generalize(gbox const & ant, gbox const & con);

    std::size_t hash() const;

    Enode * get_var(unsigned const i) const { return m_vars[i]; }
    std::vector<Enode *> const & get_vars() const { return m_vars; }
    std::vector<Enode *> get_vars() { return m_vars; }

    ibex::IntervalVector get_values() { return m_vec; }
    ibex::IntervalVector const & get_values() const { return m_vec; }
    ibex::Interval get_value(unsigned const i) { return m_vec[i]; }
    void set_value(unsigned const i, double const lb, double const ub);
    void set_value(unsigned const i, ibex::Interval const & iv);
    ibex::Interval const & get_value(unsigned const i) const { return m_vec[i]; }

    bool get_gbit(unsigned const i) const { return m_gbits[i]; }
    std::vector<bool> const & get_gbits() const { return m_gbits; }
    std::vector<bool> get_gbits() { return m_gbits; }
    void set_gbit(unsigned const i, bool const gbit) { m_gbits[i] = gbit; }

    std::unordered_set<gbox> set_minus(gbox const & gb) const;

    ibex::Interval const & operator[](unsigned const i) const { return m_vec[i]; }
    ibex::Interval const & operator[](unsigned const i) { return m_vec[i]; }
    friend std::ostream& operator<<(std::ostream& out, gbox const & gb);
};

std::unordered_set<gbox> set_minus(gbox gb1, gbox const & gb2);
gbox generalize(gbox gb, gbox const & ant, gbox const & con);
gbox intersection(gbox gb1, gbox const & gb2);
int find_merge_dim(gbox const & b1, gbox const & b2);
gbox merge(gbox b1, gbox const & b2, int const merge_dim);

std::ostream& operator<<(std::ostream& out, gbox const & gb);
}  // namespace dreal

namespace std {
template <>
struct hash<dreal::gbox> {
    std::size_t operator()(dreal::gbox const & gb) const { return gb.hash(); }
};
}  // namespace std
