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
using std::random_device;
using std::unordered_set;
using std::vector;
using std::numeric_limits;
using std::dynamic_pointer_cast;
using std::tuple;

namespace dreal {

void add_interval(PicoSAT * ps, Grid const & g, Enode * v, double const l, double const u) {
    // Add v >= l
    picosat_add(ps, g.lookup_ge(v, l)); picosat_add(ps, 0);
    // Add v <= u
    picosat_add(ps, g.lookup_le(v, u)); picosat_add(ps, 0);
}

void add_vector(PicoSAT * ps, vector<int> const & vec) {
    for (int const l : vec) {
        picosat_add(ps, l);
    }
}

box build_box_from_sat_model(PicoSAT * ps, Grid const & g, box b) {
    // TODO(soonhok): implement this
    return b;
}

// Given a box `b`, Add `!b`
void add_learned_clause(PicoSAT * ps, Grid const & g, box const & b) {
    for (Enode * v : b.get_vars()) {
        double const l = b[v].lb();
        double const u = b[v].ub();
        // !(l <= v <= u) --> !(l <= v  /\   v <= u)
        //                --> !(l <= v) \/ !(v <= u)
        picosat_add(ps, -g.lookup_le(l, v));
        picosat_add(ps, -g.lookup_le(v, u));
    }
    picosat_add(ps, 0);
}

// Given boxes `b1` and `b2`, add `b1 => b2`, that is `!b1 \/ b2`
void add_learned_clause(PicoSAT * ps, Grid const & g, box const & b1, box const & b2) {
    // Step 1. Collect literals in b1
    vector<int> b1_lits;
    auto vars = b1.get_vars();
    for (Enode * v : vars) {
        double const l = b1[v].lb();
        double const u = b1[v].ub();
        // (l <= v <= u) --> (l <= v  /\   v <= u)
        b1_lits.push_back(g.lookup_le(l, v));
        b1_lits.push_back(g.lookup_le(v, u));
    }

    for (Enode * v : vars) {
        // !b1 --> (l <= v)
        for (int l : b1_lits) {
            picosat_add(ps, -l);
        }
        double const l = b2[v].lb();
        picosat_add(ps, g.lookup_le(l, v));
        picosat_add(ps, 0);

        // !b1 --> (v <= u)
        for (int l : b1_lits) {
            picosat_add(ps, -l);
        }
        double const u = b2[v].ub();
        picosat_add(ps, g.lookup_le(v, u));
        picosat_add(ps, 0);
    }
}

// scoped_vec<std::shared_ptr<constraint>> const & ctrs
box gsat_icp::solve(box b, contractor & ctc, SMTConfig & config) {
    unordered_set<std::shared_ptr<constraint>> used_constraints;
    box initial_box(b);
    Grid g(b);
    PicoSAT * ps = picosat_init();

    // ============== INITIALIZATION BEGIN ===============
    // Add initial box
    for (Enode * v : b.get_vars()) {
        double const l = b[v].lb();
        double const u = b[v].ub();
        add_interval(ps, g, v, l, u);
    }
    // Add NO Bounds (Before PUSH)
    vector<int> no_bounds = g.get_push_nobounds_formula();
    no_bounds.pop_back();  // HACK
    add_vector(ps, no_bounds);
    picosat_push(ps);  // <------ PICOSAT_PUSH
    // Add Bounds (After PUSH)
    vector<int> bounds = g.get_push_nobounds_formula();
    bounds.pop_back();  // HACK
    add_vector(ps, bounds);
    // ============== INITIALIZATION END =================

    while (true) {
        int const ret = picosat_sat(ps, -1);
        if (ret == PICOSAT_SATISFIABLE) {
            // ============== PRUNING BEGIN ===============
            // Find a satisfying Boolean Assignment
            b = build_box_from_sat_model(ps, g, initial_box);
            assert(!b.is_empty());

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
                add_learned_clause(ps, g, old_box);
                continue;  // Go back to the beginning of the while loop, to get another assignment.
            } else {
                // Case 2: B ==> B'
                //  - Learn B => B'
                //  - Add Points in B'
                assert(b.is_subset(old_box));
                add_learned_clause(ps, g, old_box, b);
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
            g.add_point(br_var, br_pt);
        } else {
            // i < 0, which indicates it's not possible to bisect.
            break;  // exit the while loop, since we have a small enough box (DELTA-SAT);
        }
        // ============== BRANCHING END ===============
    }  // end of while loop
    picosat_reset(ps);
    return b;
}
}  // namespace dreal
