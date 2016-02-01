/*********************************************************************
Author: Soonho Kong <soonhok@cs.cmu.edu>

dReal -- Copyright (C) 2013 - 2016, Soonho Kong, Sicun Gao, and Edmund Clarke

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

#include <functional>
#include <iostream>
#include <limits>
#include <map>
#include <random>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include "util/logging.h"
#include "util/box.h"
#include "util/stat.h"
#include "contractor/contractor.h"
#include "opensmt/smtsolvers/SMTConfig.h"
#include "picosat/picosat.h"

namespace dreal {
template <class T>
inline void hash_combine(std::size_t& seed, const T& v) {
    std::hash<T> hasher;
    seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}
}

namespace std {
template<>
struct hash<std::tuple<Enode*, double, bool>> {
    size_t operator () (const std::tuple<Enode*, double, bool> & v) const {
        std::size_t s = 23;
        dreal::hash_combine<Enode*>(s, std::get<0>(v));
        dreal::hash_combine<double>(s, std::get<1>(v));
        dreal::hash_combine<bool>(s, std::get<2>(v));
        return s;
    }
};
template<>
struct equal_to<std::tuple<Enode*, double, bool>> {
    bool operator() (const std::tuple<Enode*, double, bool> & v1, const std::tuple<Enode*, double, bool> & v2) const {
        return std::get<0>(v1) == std::get<0>(v2)
            && std::get<1>(v1) == std::get<1>(v2)
            && std::get<2>(v1) == std::get<2>(v2);
    }
};
}  // namespace std

namespace dreal {
class pred_abs {
private:
    int m_num_vars = 0;
    std::unordered_map<int, std::tuple<Enode*, double, bool>> m_con_map;
    std::unordered_map<std::tuple<Enode*, double, bool>, int> m_abs_map;

public:
    int add(Enode* v, double const bound, bool const le) {
        auto p = std::make_tuple(v, bound, le);
        auto const it = m_abs_map.find(p);
        if (it == m_abs_map.end()) {
            ++m_num_vars;
            m_con_map.emplace(m_num_vars, p);
            m_abs_map.emplace(p, m_num_vars);
            return m_num_vars;
        } else {
            return it->second;
        }
    }
    int get_num_vars() const {
        return m_num_vars;
    }
    int lookup(Enode * v, double const bound, bool const le) const {
        return m_abs_map.at(std::make_tuple(v, bound, le));
    }

    std::tuple<Enode*, double, bool> lookup(int const n) const {
        return m_con_map.at(n);
    }
    void debug_print() const {
        for (int i = 1; i <= m_num_vars; ++i) {
            auto t = lookup(i);
            DREAL_LOG_FATAL << "B" << i << "\t <---> \t"
                            << std::get<0>(t)
                            << (std::get<2>(t) ? " <= " : " >= ")
                            << std::get<1>(t);
        }
    }
};

class picosat_wrapper {
private:
    PicoSAT * m_psat;
    pred_abs m_store;
    std::vector<int> m_pmodel;

public:
    picosat_wrapper();
    ~picosat_wrapper();

    // Given a variable `v` and two constants `l` and `u` where l < u holds.
    void add_ordering(box const & b, Enode *v, double const l, double const u);
    void add_ordering(box const & b, std::unordered_set<Enode*> const & used_vars, Enode * v, double const l, double const u);
    void add_axiom(box const & b, Enode *v, double const c);
    void add_axiom(box const & b, std::unordered_set<Enode*> const & used_vars, Enode *v, double const c);

    // Add: v <= bound
    void add_le(Enode * v, double const bound);

    // Add: bound <= v
    void add_le(double const bound, Enode* v);

    // Add: v >= bound
    void add_ge(Enode * v, double const bound);

    // Add: bound >= v
    void add_ge(double const bound, Enode* v);

    // // Add: lb <= v <= ub
    // void add_intv(double const lb, Enode* v, double const ub);

    void add_box(box const & b);

    // Add blocking clause ¬B, but generalize it using used_vars
    void add_generalized_blocking_box(box const & b, std::unordered_set<Enode *> const & used_vars);

    // Add blocking clause B1 => B2, but generalize it using used_vars
    void add_generalized_blocking_box(box const & b1, box const & b2, std::unordered_set<Enode *> const & used_vars);

    // Add B => l1 \/ l2 \/ l3 \/ l4
    void add_imply(box const & b, std::vector<Enode *> const & used_vars, int const l1, int const l2 = 0, int const l3 = 0, int const l4 = 0);
    void add_imply(box const & b, std::unordered_set<Enode *> const & used_vars, int const l1, int const l2 = 0, int const l3 = 0, int const l4 = 0);

    // Add B => (B[v].lb <= v <= m) xor (m <= v <= B[v].ub)
    void add_branching(box const & b, Enode * v, double const m);

    void block_current_assignment();

    int check_sat();

    // Precondition: check_sat() == PICOSAT_SATISFIABLE
    // Reduce the given box b into a smaller box using SAT model
    box reduce_using_model(box b);

    void debug_print() const;
};
}  // namespace dreal
