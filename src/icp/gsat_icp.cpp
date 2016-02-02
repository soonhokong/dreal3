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

#include <functional>
#include <iostream>
#include <limits>
#include <map>
#include <random>
#include <unordered_set>
#include <vector>
#include <tuple>
#include "icp/gsat_icp.h"
#include "util/logging.h"
#include "util/stat.h"
#include "icp/point_grid.h"
#include "picosat/picosat.h"

using std::cerr;
using std::cout;
using std::endl;
using std::get;
using std::map;
using std::set;
using std::random_device;
using std::unordered_set;
using std::vector;
using std::numeric_limits;
using std::dynamic_pointer_cast;
using std::tuple;

namespace dreal {

gsat_icp::gsat_icp(box const & b) : m_initial_box(b), m_grid(b) {
    m_ps = picosat_init();
    picosat_save_original_clauses(m_ps);
}

gsat_icp::~gsat_icp() {
    picosat_reset(m_ps);
}

    // add (l <= v) /\ (v <= u)
void gsat_icp::add_interval(Enode * v, double const l, double const u) {
    picosat_add(m_ps, m_grid.lookup_le(l, v)); picosat_add(m_ps, 0);  // Add l <= v
    picosat_add(m_ps, m_grid.lookup_le(v, u)); picosat_add(m_ps, 0);  // Add v <= u
}

void gsat_icp::add_vector(vector<int> const & vec) {
    cerr << "add_vector:";
    for (int const l : vec) {
        cerr << " " << l;
        picosat_add(m_ps, l);
    }
    cerr << endl;
}

box gsat_icp::build_box_from_sat_model() {
    box b = m_initial_box;
    for (Enode * v : b.get_vars()) {
        set<double> const & point_row = m_grid.get_point_row(v);

        // Find an upperbound for a b[v]
        for (auto it = point_row.cbegin(); it != point_row.cend(); ++it) {
            // We scan the point row from left to right (i.e. 1, 2, 3, 4, 5)
            // Find the first element that holds in the SAT assignment.
            // This should give us an upperbound for a variable. For example,
            //
            // (v <= 1), (v <= 2), (v <= 3), (v <= 4), (v <= 5)
            //
            double const p = *it;
            int const l = m_grid.lookup_le(v, p);
            int const r = picosat_deref_partial(m_ps, l);
            if (r == 1) {
                b[v] = ibex::Interval(b[v].lb(), p);
                break;  // exit this for-loop
            }
        }
        // Find a lowerbound for a b[v]
        for (auto it = point_row.crbegin(); it != point_row.crend(); ++it) {
            // We scan the point row from right to left (i.e. 5, 4, 3, 2, 1)
            // Find the first element that holds in the SAT assignment.
            // This should give us a lowerbound for a variable. For example,
            //
            // (v >= 5), (v >= 4), (v >= 3), (v >= 2), (v >= 1)
            //
            double const p = *it;
            int const l = m_grid.lookup_ge(v, p);
            int const r = picosat_deref_partial(m_ps, l);
            if (r == 1) {
                b[v] = ibex::Interval(b[v].lb(), p);
                break;  // exit this for-loop
            }
        }
    }
    return b;
}

// Given a box `b`, Add `!b`
void gsat_icp::add_learned_clause(box const & b) {
    for (Enode * v : b.get_vars()) {
        double const l = b[v].lb();
        double const u = b[v].ub();
        // !(l <= v <= u) --> !(l <= v  /\   v <= u)
        //                --> !(l <= v) \/ !(v <= u)
        picosat_add(m_ps, -m_grid.lookup_le(l, v));
        picosat_add(m_ps, -m_grid.lookup_le(v, u));
    }
    picosat_add(m_ps, 0);
}

// Given boxes `b1` and `b2`, add `b1 => b2`, that is `!b1 \/ b2`
void gsat_icp::add_learned_clause(box const & b1, box const & b2) {
    // Step 1. Collect literals in b1
    vector<int> b1_lits;
    auto vars = b1.get_vars();
    for (Enode * v : vars) {
        double const l = b1[v].lb();
        double const u = b1[v].ub();
        // (l <= v <= u) --> (l <= v  /\   v <= u)
        b1_lits.push_back(m_grid.lookup_le(l, v));
        b1_lits.push_back(m_grid.lookup_le(v, u));
    }

    for (Enode * v : vars) {
        // !b1 --> (l <= v)
        for (int l : b1_lits) {
            picosat_add(m_ps, -l);
        }
        double const l = b2[v].lb();
        picosat_add(m_ps, m_grid.lookup_le(l, v));
        picosat_add(m_ps, 0);

        // !b1 --> (v <= u)
        for (int l : b1_lits) {
            picosat_add(m_ps, -l);
        }
        double const u = b2[v].ub();
        picosat_add(m_ps, m_grid.lookup_le(v, u));
        picosat_add(m_ps, 0);
    }
}

// scoped_vec<std::shared_ptr<constraint>> const & ctrs
box gsat_icp::solve(contractor & ctc, SMTConfig & config) {
    unordered_set<std::shared_ptr<constraint>> used_constraints;
    box b = m_initial_box;

    // ============== INITIALIZATION BEGIN ===============
    // Add initial box
    for (Enode * v : b.get_vars()) {
        double const l = b[v].lb();
        double const u = b[v].ub();
        add_interval(v, l, u);
    }
    vector<int> no_bounds = m_grid.get_push_nobounds_formula();
    add_vector(no_bounds);
    // ============== INITIALIZATION END =================

    while (true) {
        int const ret = picosat_sat(m_ps, -1);
        if (ret == PICOSAT_SATISFIABLE) {
            // ============== PRUNING BEGIN ===============
            DREAL_LOG_FATAL << "Found a Boolean assignment";
            b = build_box_from_sat_model();
            assert(!b.is_empty());
            DREAL_LOG_FATAL << "Current Box\n===========\n" << b;

            box old_box(b);
            try {
                ctc.prune(b, config);
                if (config.nra_use_stat) { config.nra_stat.increase_prune(); }
                auto const & this_used_constraints = ctc.used_constraints();
                used_constraints.insert(this_used_constraints.begin(), this_used_constraints.end());
            } catch (contractor_exception & e) {
                // Do nothing
            }
            if (b.is_empty()) {
                // Case 1: B ==> empty set
                //  - Learn clause !B
                add_learned_clause(old_box);
                continue;  // Go back to the beginning of the while loop, to get another assignment.
            } else {
                // Case 2: B ==> B'
                //  - Learn B => B'
                //  - Add Points in B'
                assert(b.is_subset(old_box));
                add_learned_clause(old_box, b);
            }
            // ============== PRUNING END   ===============
        } else if (ret == PICOSAT_UNSATISFIABLE) {
            // No satisfying Boolean Assignment
            break;  // Exit the while loop
        } else {
            assert(ret == PICOSAT_UNKNOWN);
            DREAL_LOG_FATAL << "Something is wrong with PICOSAT, which returns UNKNOWN.";
            abort();
        }
        // ============== BRANCHING BEGIN ===============
        assert(!b.is_empty());  // It should be the case B1 => B2 => ... => Bn and Bn is non-empty
        tuple<int, box, box> splits = b.bisect(config.nra_precision);
        int const i = get<0>(splits);  // i is branched dimension
        if (i >= 0) {
            if (config.nra_use_stat) { config.nra_stat.increase_branch(); }
            Enode * br_var = b.get_vars()[i];   // branched variable;
            double const br_pt = b[i].mid();  // branching point on br_var (for now, it's the mid point).
            m_grid.add_point(br_var, br_pt);
        } else {
            // i < 0, which indicates it's not possible to bisect.
            break;  // exit the while loop, since we have a small enough box (DELTA-SAT);
        }
        // ============== BRANCHING END ===============
    }  // end of while loop
    picosat_reset(m_ps);
    return b;
}
}  // namespace dreal
