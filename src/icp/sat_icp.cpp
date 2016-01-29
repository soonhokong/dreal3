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
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include "icp/sat_icp.h"
#include "util/logging.h"
#include "util/stat.h"
#include "icp/picosat_wrapper.h"

using std::cerr;
using std::cout;
using std::endl;
using std::get;
using std::map;
using std::make_tuple;
using std::random_device;
using std::tuple;
using std::unordered_map;
using std::unordered_set;
using std::vector;
using std::numeric_limits;

namespace dreal {
box sat_icp::solve(box b, contractor & ctc, SMTConfig & config) {
    thread_local static std::unordered_set<std::shared_ptr<constraint>> used_constraints;
    used_constraints.clear();
    // Step 1. Initialize SAT Solver
    picosat_wrapper pw;
    pw.add_box(b);  // Add initial interval constraints
    box const initial_box(b);
    // Step 2. Main Part
    while (true) {
        DREAL_LOG_FATAL << "\n\n";
        // pw.debug_print();
        int ret = pw.check_sat();
        if (ret == PICOSAT_SATISFIABLE) {
            DREAL_LOG_INFO << "SAT solver found a satisfying Boolean assignment";
            // Case 1: SAT solver found a satisfying Boolean assignment
            // 1.1. Concretize the satisfying Boolean assignment into a conjunction of constraints
            // Check each Boolean Variable and if partially assigned to be true, shrink the interval
            // DREAL_LOG_INFO << "store.get_num_vars() = " << store.get_num_vars();
            b = pw.reduce_using_model(initial_box);
            DREAL_LOG_FATAL << "Current Box = " << "\n"
                            << b;
            if (b.is_empty()) {
                DREAL_LOG_FATAL << "SOMETHING IS WRONG";
                abort();
            }
            // 1.2. Apply pruning operators with box B until it reaches a fixed point B'
            box old_b(b);
            try {
                ctc.prune(b, config);
                auto const this_used_constraints = ctc.used_constraints();
                used_constraints.insert(this_used_constraints.begin(), this_used_constraints.end());
            } catch (contractor_exception & e) { /* Do nothing */ }
            if (config.nra_use_stat) { config.nra_stat.increase_prune(); }
            // Collect Used Variables
            unordered_set<Enode *> used_vars;
            for (auto const & used_ctr : ctc.used_constraints()) {
                auto const this_used_vars = used_ctr->get_vars();
                used_vars.insert(this_used_vars.begin(), this_used_vars.end());
            }
            DREAL_LOG_INFO << "|USED VARS| = " << used_vars.size();

            if (b.is_empty()) {
                DREAL_LOG_FATAL << "After Pruning, it became an empty set.";
                // Case i: Pruning returns an empty box.
                // Add Blocking Clause: !old_b

                // DREAL_LOG_INFO << "BEFORE ADD BLOCKING CLAUSE";
                // pw.debug_print();

                pw.add_generalized_blocking_box(old_b, used_vars);

                // DREAL_LOG_INFO << "AFTER ADD BLOCKING CLAUSE";
                // pw.debug_print();


            } else {
                DREAL_LOG_FATAL << "After Pruning, it became a non-empty set.\n" << b;
                // Case ii: Pruning returns a non-empty box b
                if (b.max_diam() < config.nra_precision) {
                    DREAL_LOG_FATAL << "Box is small enough to stop.";
                    // Box is small enough to stop => delta-SAT
                    break;
                } else {
                    DREAL_LOG_FATAL << "Box is big: width = " << b.max_diam();
                    // Box is big, and needs to be branched. We also need to learn a clause from this pruning
                    if (config.nra_use_stat) { config.nra_stat.increase_branch(); }

                    // Pick a branching variable and branching point
                    // TODO(soonhok): b.bisect is an overkill here
                    // since it returns two boxes which are not used
                    auto const bisect_result = b.bisect(config.nra_precision);
                    Enode * br_var = b.get_vars()[get<0>(bisect_result)];
                    double const br_point = b[br_var].mid();
                    pw.add_branching(b, br_var, br_point);
                    pw.add_generalized_blocking_box(old_b, b, used_vars);
                }
            }
        } else if (ret == PICOSAT_UNSATISFIABLE) {
            DREAL_LOG_FATAL << "SAT solver failed to find a satisfying Boolean assignment";
            // Case 2: SAT solver concludes UNSAT. Return UNSAT.
            b.set_empty();
            break;
        } else {
            assert(ret == PICOSAT_UNKNOWN);
            DREAL_LOG_FATAL << "SAT Solver failed.";
        }
    }
    ctc.set_used_constraints(used_constraints);
    return b;
}
}  // namespace dreal
