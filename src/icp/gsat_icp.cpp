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
    cerr << "ADD VECTOR:";
    if (vec.empty()) {
        cerr << " NOTHING";
    } else {
        for (int const l : vec) {
            picosat_add(m_ps, l);
            cerr << " " << l;
        }
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
            // int const r = picosat_deref_partial(m_ps, l);
            int const r = picosat_deref(m_ps, l);
            DREAL_LOG_FATAL << "\t\t" << v << " <= " << p << " " << r;
            if (r == 1) {
                b[v] = ibex::Interval(b[v].lb(), p);
                DREAL_LOG_FATAL << v << " <= " << p;
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
            // int const r = picosat_deref_partial(m_ps, l);
            int const r = picosat_deref(m_ps, l);
            DREAL_LOG_FATAL << "\t\t" << v << " >= " << p << " " << r;
            if (r == 1) {
                b[v] = ibex::Interval(p, b[v].ub());
                DREAL_LOG_FATAL << p << " <= " << v;
                break;  // exit this for-loop
            }
        }
    }
    return b;
}

// Given a box `b`, Add `!b`
void gsat_icp::add_learned_clause(box const & b) {
    add_imply(b);
}

// Given boxes `b1` and `b2`, add `b1 => b2`, that is `!b1 \/ b2`
void gsat_icp::add_learned_clause(box const & b1, box const & b2) {
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
        double const lb = b2[v].lb();
        double const ub = b2[v].ub();
        m_grid.add_point(v, lb);
        m_grid.add_point(v, ub);

        // !b1 --> (lb <= v)
        for (int l : b1_lits) {
            picosat_add(m_ps, -l);
        }
        picosat_add(m_ps, m_grid.lookup_le(lb, v));
        picosat_add(m_ps, 0);

        // !b1 --> (v <= ub)
        for (int l : b1_lits) {
            picosat_add(m_ps, -l);
        }
        picosat_add(m_ps, m_grid.lookup_le(v, ub));
        picosat_add(m_ps, 0);
    }
}

// Add B => l1 \/ l2 \/ l3 \/ l4
void gsat_icp::add_imply(box const & b, int const l1, int const l2, int const l3, int const l4) {
    // TODO(soonhok): delete the following (it's for debugging purposes)
    vector<int> added_clause;
    for (Enode * v : b.get_vars()) {
        double const l = b[v].lb();
        double const u = b[v].ub();
        //     !((l <= v) /\  (v <= u))
        // -->  !(l <= v) \/ !(v <= u)
        picosat_add(m_ps, -m_grid.lookup_le(l, v));
        picosat_add(m_ps, -m_grid.lookup_le(v, u));
        added_clause.push_back(-m_grid.lookup_le(l, v));
        added_clause.push_back(-m_grid.lookup_le(v, u));
    }
    if(l1) {
        picosat_add(m_ps, l1);
        added_clause.push_back(l1);
        if (l2) {
            picosat_add(m_ps, l2);
            added_clause.push_back(l2);
            if (l3) {
                picosat_add(m_ps, l3);
                added_clause.push_back(l3);
                if (l4) {
                    picosat_add(m_ps, l4);
                    added_clause.push_back(l4);
                }
            }
        }
    }
    picosat_add(m_ps, 0);
    added_clause.push_back(0);
    cerr << "ADD_IMPLY: "; m_grid.debug_print_clause(added_clause);
}

// Add l1 \/ l2 \/ l3 \/ l4
void gsat_icp::add_imply(int const l1, int const l2, int const l3, int const l4) {
    assert(l1);
    picosat_add(m_ps, l1);
    if (l2) {
        picosat_add(m_ps, l2);
        if (l3) {
            picosat_add(m_ps, l3);
            if (l4) {
                picosat_add(m_ps, l4);
            }
        }
    }
    picosat_add(m_ps, 0);
}

void gsat_icp::add_branch(box const & b, Enode * v, double const p) {
    // Introduce a branch point, v = p
    m_grid.add_point(v, p);

    //  B => (p <= v) \/ (v <= p)
    // !B \/ (p <= v) \/ (v <= p)
    int const l1 = m_grid.lookup_le(p, v);
    int const l2 = m_grid.lookup_le(v, p);
    add_imply(b, l1, l2);
}

// scoped_vec<std::shared_ptr<constraint>> const & ctrs
box gsat_icp::solve(contractor & ctc, SMTConfig & config) {
    unordered_set<std::shared_ptr<constraint>> used_constraints;
    box b = m_initial_box;

    while (true) {
        DREAL_LOG_FATAL << "\n\n\n";
        vector<int> no_bounds = m_grid.get_push_nobounds_formula();
        add_vector(no_bounds);
        m_grid.debug_print();
        int const ret = picosat_sat(m_ps, -1);
        if (ret == PICOSAT_SATISFIABLE) {
            // ============== PRUNING BEGIN ===============
            DREAL_LOG_FATAL << "Found a Boolean assignment";
            b = build_box_from_sat_model();
            assert(!b.is_empty());
            DREAL_LOG_FATAL << "Current Box";
            DREAL_LOG_FATAL << "===========";
            DREAL_LOG_FATAL << b;

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
                DREAL_LOG_FATAL << "After pruning, we have an empty set";
                // Case 1: B ==> empty set
                //  - Learn clause !B
                add_learned_clause(old_box);
                continue;  // Go back to the beginning of the while loop, to get another assignment.
            } else {
                DREAL_LOG_FATAL << "After pruning, we have a non-empty set";
                DREAL_LOG_FATAL << b;
                // Case 2: B ==> B'
                //  - Learn B => B' (+ add_point in B')
                assert(b.is_subset(old_box));
                if (old_box != b) {
                    add_learned_clause(old_box, b);
                }
            }
            // ============== PRUNING END   ===============
        } else if (ret == PICOSAT_UNSATISFIABLE) {
            DREAL_LOG_FATAL << "No SAT Assignment";
            b.set_empty();
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
            DREAL_LOG_FATAL << "Branching on " << br_var << " = " << br_pt;
            add_branch(b, br_var, br_pt);
        } else {
            // i < 0, which indicates it's not possible to bisect.
            DREAL_LOG_FATAL << "Failed to branch (|b = " << b.max_diam() << ")";
            break;  // exit the while loop, since we have a small enough box (DELTA-SAT);
        }
        // ============== BRANCHING END ===============
    }  // end of while loop
    return b;
}
}  // namespace dreal
