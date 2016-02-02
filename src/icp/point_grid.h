/*********************************************************************
Author: Sicun Gao <sicung@cs.cmu.edu>
        Soonho Kong <soonhok@cs.cmu.edu>

dReal -- Copyright (C) 2013 - 2015, Soonho Kong, Sicun Gao, and Edmund Clarke

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
#include "util/hash_combine.h"

namespace std {
template<>
struct hash<std::pair<Enode*, double>> {
    size_t operator () (const std::pair<Enode*, double> & v) const {
        std::size_t s = 23;
        dreal::hash_combine<Enode*>(s, v.first);
        dreal::hash_combine<double>(s, v.second);
        return s;
    }
};
template<>
struct equal_to<std::pair<Enode*, double>> {
    bool operator() (const std::pair<Enode*, double> & v1, const std::pair<Enode*, double> & v2) const {
        return v1.first == v2.first && v1.second == v2.second;
    }
};
}  // namespace std


namespace dreal {

    class Grid { //Point Matrix
    private:
        //top_lit keeps the index of the last literal
        unsigned top_lit;
        //we need a map from a variable to the ordered points in its dimension
        std::unordered_map<Enode*, std::set<double>> point_rows;

        //encoding: lower bound constraints x>c use odd numbers. upper bounds use even numbers.
        std::unordered_map<Enode*, std::vector<int>>	lb_lits;
        std::unordered_map<Enode*, std::vector<int>>	ub_lits;

        /* a global mapping from (var, point) to lit index.
           lb and ub are automatically inferred from whether the index is odd or even
           by default it returns the lb constraint. Use +1 for the ub constraint. */
        std::unordered_map<std::pair<Enode*, double>, int> lb_lit_map;

        //clauses
        std::vector<std::vector<int>> linear_clauses;   //linear ordering on bounds
        std::vector<std::vector<int>> lu_clauses;       //lower and upper bounds should be consistent
        std::vector<int>              full_lb_clauses;
        std::vector<int>              full_ub_clauses;

        //new clauses that should be pushed
        std::vector<std::vector<int>> push_linear_clauses;
        std::vector<std::vector<int>> push_lu_clauses;

        //prepare for sat
        std::vector<int> current_formula;
        std::vector<int> push_formula;

    public:
        Grid(box const &);
        void add_box(box const &);
        void add_point(Enode *, double const);
	void add_initial_points(Enode *, double const, double const);
        std::set<double> get_point_row(Enode * v) { return point_rows[v]; }
        std::set<double> const & get_point_row(Enode * v) const { return point_rows.at(v); }
        inline int lookup_le(Enode * v, double const p) const {
            return lb_lit_map.at(std::make_pair(v, p));
        }
        inline int lookup_ge(Enode * v, double const p) const {
            return lookup_le(v, p) + 1;
        }

        // p <= v  ==  v >= p
        inline int lookup_le(double const p, Enode * v) const {
            return lookup_ge(v, p);
        }
        // p >= v  ==  v <= p
        inline int lookup_ge(double const p, Enode * v) const {
            return lookup_le(v, p) + 1;
        }

        inline void build_current_formula() {
            current_formula.clear();
            for (auto cl : linear_clauses)
                for (auto i : cl)
                    current_formula.push_back(i);
            for (auto cl : lu_clauses)
                for (auto i : cl)
                    current_formula.push_back(i);
            for (auto row : lb_lits) {
                for (auto i : row.second)
                    current_formula.push_back(i);
                current_formula.push_back(0);
            }
            for (auto row : ub_lits) {
                for (auto i : row.second)
                    current_formula.push_back(i);
                current_formula.push_back(0);
            }
        }
        inline void build_push_formula() {
            push_formula.clear();
            for (auto cl : push_linear_clauses)
                for (auto i : cl)
                    push_formula.push_back(i);
            for (auto cl : push_lu_clauses)
                for (auto i : cl)
                    push_formula.push_back(i);
            for (auto row : lb_lits) {
                for (auto i : row.second)
                    push_formula.push_back(i);
                push_formula.push_back(0);
            }
            for (auto row : ub_lits) {
                for (auto i : row.second)
                    push_formula.push_back(i);
                push_formula.push_back(0);
            }
        }
        inline void build_push_nobounds_formula() {
            push_formula.clear();
            for (auto cl : push_linear_clauses)
                for (auto i : cl)
                    push_formula.push_back(i);
            for (auto cl : push_lu_clauses)
                for (auto i : cl)
                    push_formula.push_back(i);
        }
        inline void build_push_bounds_only_formula() {
            push_formula.clear();
            for (auto row : lb_lits) {
                for (auto i : row.second)
                    push_formula.push_back(i);
                push_formula.push_back(0);
            }
            for (auto row : ub_lits) {
                for (auto i : row.second)
                    push_formula.push_back(i);
                push_formula.push_back(0);
            }
        }

        inline std::vector<int> const & get_current_formula() {
            build_current_formula();
            return current_formula;
        }
        inline std::vector<int> const & get_push_formula() {
            build_push_formula();
            return push_formula;
        }
        inline std::vector<int> const & get_push_nobounds_formula() {
            build_push_nobounds_formula();
            return push_formula;
        }
        inline std::vector<int> const & get_push_bounds_only_formula() {
            build_push_bounds_only_formula();
            return push_formula;
        }

        void debug_print_clause(std::vector<int> const & c) const;
        void debug_print() const;
    };
}
