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
#include <map>
#include <random>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include "icp/sat_icp.h"
#include "util/logging.h"
#include "util/stat.h"

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

namespace dreal {
template <class T>
inline void hash_combine(std::size_t& seed, const T& v) {
    std::hash<T> hasher;
    seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}
}

namespace std {
template<>
struct hash<tuple<Enode*, double, double>> {
    size_t operator () (const tuple<Enode*, double, double> & v) const {
        std::size_t s = 23;
        dreal::hash_combine<Enode*>(s, get<0>(v));
        dreal::hash_combine<double>(s, get<1>(v));
        dreal::hash_combine<double>(s, get<2>(v));
        return s;
    }
};
template<>
struct equal_to<tuple<Enode*, double, double>> {
    bool operator() (const tuple<Enode*, double, double> & v1, const tuple<Enode*, double, double> & v2) const {
        return get<0>(v1) == get<0>(v2) &&
            get<1>(v1) == get<1>(v2) &&
            get<2>(v1) == get<2>(v2);
    }
};
}  // namespace std

namespace dreal {

class pred_abs {
private:
    int m_num_vars = 0;
    unordered_map<int, tuple<Enode*, double, double>> m_con_map;
    unordered_map<tuple<Enode*, double, double>, int> m_abs_map;

public:
    int add(Enode* e, double lb, double ub) {
        auto p = make_tuple(e, lb, ub);
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
    int lookup(Enode * e, double lb, double ub) const {
        return m_abs_map.at(make_tuple(e, lb, ub));
    }

    tuple<Enode*, double, double> lookup(int const n) const {
        return m_con_map.at(n);
    }
    void debug_print() const {
        for (int i = 1; i <= m_num_vars; ++i) {
            auto t = lookup(i);
            DREAL_LOG_FATAL << "B" << i << "\t <---> \t"
                            << get<0>(t)
                            << " [" << get<1>(t) << ", " << get<2>(t) << "]";
        }
    }
};

box sat_icp::solve(box b, contractor & ctc, SMTConfig & config) {
    // Step 1. Initialize SAT Solver
    //  - For each constraint p_i, introduce a corresponding Boolean variable b_i to SAT solver.
    //  - Add initial clauses $b_1, \dots, b_n$ to SAT solver.

    PicoSAT * psat = picosat_init();
    picosat_save_original_clauses(psat);  // to use picosat_deref_partial function
    pred_abs store;

    //  - Add initial box (unit) clauses
    vector<int> vars_of_B;
    auto vars = b.get_vars();
    for (unsigned i = 0; i < b.size(); ++i) {
        Enode * e = vars[i];
        double const lb = b[e].lb();
        double const ub = b[e].ub();
        picosat_add(psat, store.add(e, lb, ub));
        picosat_add(psat, 0);
    }

    // Step 2. Main Part
    while (true) {
        DREAL_LOG_FATAL << "\n\n\n\n";
        // Ask SAT Solver for a satisfying assignment
        DREAL_LOG_FATAL << "===============================";
        picosat_print(psat, stderr);
        DREAL_LOG_FATAL << "===============================";
        store.debug_print();
        DREAL_LOG_FATAL << "===============================";

        int ret = picosat_sat(psat, -1);
        if (ret == PICOSAT_SATISFIABLE) {
            DREAL_LOG_FATAL << "SAT solver found a satisfying Boolean assignment";
            // Case 1: SAT solver found a satisfying Boolean assignment
            // 1.1. Concretize the satisfying Boolean assignment into a conjunction of constraints
            // Check each Boolean Variable and if partially assigned to be true, shrink the interval
            DREAL_LOG_FATAL << "store.get_num_vars() = " << store.get_num_vars();

            for (int i = 1; i <= store.get_num_vars(); i++) {
                int const r = picosat_deref_partial(psat, i);
                if (r == 1) {
                    // Because i-th Boolean Variable is set
                    // Intersect the current interval with the corresponding interval of i-th Boolean variable
                    tuple<Enode*, double, double> intv_ctr = store.lookup(i);
                    Enode * e = get<0>(intv_ctr);
                    double const lb = get<1>(intv_ctr);
                    double const ub = get<2>(intv_ctr);
                    b[e] &= ibex::Interval(lb, ub);

                    auto t = store.lookup(i);
                    DREAL_LOG_FATAL << "SAT MODEL: " << i << "\t"
                                    << get<0>(t) << " [" << get<1>(t) << ", " << get<2>(t) << "]";
                } else {
                    DREAL_LOG_FATAL << "SAT MODEL: " << i << "\t" << "RESULT = " << r;
                }
            }
            DREAL_LOG_FATAL << "Current Box = " << "\n"
                            << b;
            // 1.2. Apply pruning operators with box B until it reaches a fixed point B'
            box old_b(b);
            ctc.prune(b, config);
            if (b.is_empty()) {
                DREAL_LOG_FATAL << "After Pruning, it became an empty set.";
                // Case i: Pruning returns an empty box.
                // 1.2.i.1. Collect Used Variables
                unordered_set<Enode *> used_vars;
                used_vars.clear();
                for (auto used_ctr : ctc.used_constraints()) {
                    auto this_used_vars = used_ctr->get_vars();
                    used_vars.insert(this_used_vars.begin(), this_used_vars.end());
                }
                if (used_vars.size() < b.size()) {
                    // 1.2.i.2. Add a blocking clause
                    for (Enode * var : used_vars) {
                        // Add the Boolean literal
                        ibex::Interval const & itv = b[var];
                        double const lb = itv.lb();
                        double const ub = itv.ub();
                        int lit = store.add(var, lb, ub);
                        picosat_add(psat, -lit);
                    }
                    picosat_add(psat, 0);
                } else {
                    assert(used_vars.size() == b.size());
                    // There is no generalization, there is no use of learning a clause since it will not be used later.
                }
                break;
            } else {
                DREAL_LOG_FATAL << "After Pruning, it became a non-empty set.";
                DREAL_LOG_FATAL << b;
                // Case ii: Pruning returns a non-empty box b
                if (b.max_diam() < config.nra_precision) {
                    DREAL_LOG_FATAL << "Box is small enough to stop.";
                    // Box is small enough to stop => delta-SAT
                    break;
                } else {
                    DREAL_LOG_FATAL << "Box is is big: " << b.max_diam();
                    // Box is big, and needs to be branched. We also need to learn a clause from this pruning
                    // ii.1. Pick a branching variable
                    // TODO(soonhok): b.bisect is an overkill here since it returns two boxes which are not necessary
                    auto const bisect_result = b.bisect(config.nra_precision);
                    Enode * branching_var = b.get_vars()[get<0>(bisect_result)];
                    // ii.2. Pick a branching point
                    ibex::Interval const I0 = b[branching_var];
                    double const branching_point = I0.mid();
                    double const lb = I0.lb();
                    double const ub = I0.ub();

                    DREAL_LOG_FATAL << "Branching Variable: " << branching_var
                                    << " [" << lb << ", " << branching_point << ", " << ub << "]";

                    // ii.3. Introduce boolean variables I1 and I2 which bisect the original interval I0
                    int var_i0 = store.add(branching_var, lb, ub);
                    int var_i1 = store.add(branching_var, lb, branching_point);
                    int var_i2 = store.add(branching_var, branching_point, ub);
                    // ii.4. Add B => (I1 xor I2) ==> (B => I1 \/ I2) /\ (B => \/ !I1 \/ !I2)
                    //                            ==> (!B \/ I1 \/ I2) /\ (!B \/ !I1 \/ !I2)
                    vector<int> vars_of_B;
                    auto vars = b.get_vars();
                    for (unsigned i = 0; i < b.size(); ++i) {
                        Enode * e = vars[i];
                        double const lb = b[e].lb();
                        double const ub = b[e].ub();
                        vars_of_B.push_back(store.add(e, lb, ub));
                    }
                    // Add (!B \/ I1 \/ I2)
                    for (int var : vars_of_B) {
                        picosat_add(psat, -var);
                    }
                    picosat_add(psat, var_i1);
                    picosat_add(psat, var_i2);
                    picosat_add(psat, 0);
                    // Add (!B \/ !I1 \/ !I2)
                    for (int var : vars_of_B) {
                        picosat_add(psat, -var);
                    }
                    picosat_add(psat, -var_i1);
                    picosat_add(psat, -var_i2);
                    picosat_add(psat, 0);
                    // ii.5. Add partial ordering among I0, I1, and I2
                    //           I1 => I0    ==>   !I1 \/ I0
                    //           I2 => I0    ==>   !I2 \/ I0
                    // Add !I1 \/ I0
                    picosat_add(psat, -var_i1);
                    picosat_add(psat, var_i0);
                    picosat_add(psat, 0);
                    // Add !I2 \/ I0
                    picosat_add(psat, -var_i2);
                    picosat_add(psat, var_i0);
                    picosat_add(psat, 0);

                    // ii.6. Add ordering between old and new interval in each dimension.
                    for (unsigned i = 0; i < b.size(); i++) {
                        auto const & old_intv = old_b[i];
                        auto const & intv = b[i];
                        if (old_intv != intv) {
                            int const old_var = store.add(vars[i], old_intv.lb(), old_intv.ub());
                            int const var     = store.add(vars[i],     intv.lb(),     intv.ub());
                            //  old_var => var
                            // !olb_var \/ var
                            picosat_add(psat, -old_var);
                            picosat_add(psat,      var);
                            picosat_add(psat,        0);
                        }
                    }
                }
            }
        } else if (ret == PICOSAT_UNSATISFIABLE) {
            DREAL_LOG_FATAL << "SAT solver failed to find a satisfying Boolean assignment";
            // Case 2: SAT solver concludes UNSAT. Return UNSAT.
            b.set_empty();
            // TODO(soonhok): need to set up explanation for this UNSAT
            break;
        } else {
            assert(ret == PICOSAT_UNKNOWN);
            DREAL_LOG_FATAL << "SAT Solver failed.";
        }
    }
    return b;
}
}  // namespace dreal
