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
template <class T>
inline void hash_combine(std::size_t& seed, const T& v) {
    std::hash<T> hasher;
    seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}
}

namespace std {
template<>
struct hash<tuple<Enode*, double>> {
    size_t operator () (const tuple<Enode*, double> & v) const {
        std::size_t s = 23;
        dreal::hash_combine<Enode*>(s, get<0>(v));
        dreal::hash_combine<double>(s, get<1>(v));
        return s;
    }
};
template<>
struct equal_to<tuple<Enode*, double>> {
    bool operator() (const tuple<Enode*, double> & v1, const tuple<Enode*, double> & v2) const {
        return get<0>(v1) == get<0>(v2) && get<1>(v1) == get<1>(v2);
    }
};
}  // namespace std

namespace dreal {

class pred_abs {
private:
    int m_num_vars = 0;
    unordered_map<int, tuple<Enode*, double>> m_con_map;
    unordered_map<tuple<Enode*, double>, int> m_abs_map;

public:
    int add(Enode* v, double bound) {
        auto p = make_tuple(v, bound);
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
    int lookup(Enode * v, double bound,) const {
        return m_abs_map.at(make_tuple(v, bound));
    }

    tuple<Enode*, double> lookup(int const n) const {
        return m_con_map.at(n);
    }
    void debug_print() const {
        for (int i = 1; i <= m_num_vars; ++i) {
            auto t = lookup(i);
            DREAL_LOG_INFO << "B" << i << "\t <---> \t"
                           << get<0>(t) " <= " << get<1>(t);
        }
    }
};

class picosat_wrapper {
private:
    PicoSAT * m_psat;
    pred_abs m_store;

public:
    picosat_wrapper() {
        m_psat = picosat_init();
        picosat_save_original_clauses(m_psat);  // to use picosat_deref_partial function
    }
    ~picosat_wrapper() {
        picosat_reset(m_psat);
    }
    // Add: v <= bound
    void add_le(Enode * v, double const bound) {
        picosat_add(m_psat, m_store.add(v, bound));
        picosat_add(m_psat, 0);
    }
    // Add: bound <= v
    void add_le(double const bound, Enode* v) {
        // Note: Strictly, !(x<=v) is (v < x) but we consider this as
        // (v <= x). This should be OK since we use numerical alg
        // anyway.
        picosat_add(m_psat, -m_store.add(v, bound));
        picosat_add(m_psat, 0);
    }
    // Add: lb <= v <= ub
    void add_intv(double const lb, Enode* v, double const ub) {
        add_le(lb, v); add_le(v, ub);
    }
    // Add: !(lb <= v /\ v <= ub)
    void add_neg_intv(double const lb, Enode* v, double const ub) {
        // That is, v <= lb \/ !(v <= ub)
        picosat_add(m_psat, m_store.add(v, lb));
        picosat_add(m_psat, -m_store.add(v, ub));
        picosat_add(m_psat, 0);
    }
    void add_box(box const & b) {
        auto const & vars = b.get_vars();
        for (unsigned i = 0; i < b.size(); ++i) {
            Enode * v = vars[i];
            double const lb = b[v].lb();
            double const ub = b[v].ub();
            add_intv(lb, v, ub);  // lb <= v <= ub
        }
    }
    // Add blocking clause ¬B, but generalize it using used_vars
    void add_generalized_blocking_box(box const & b, unordered_set<Enode *> const & used_vars) {
        for (Enode * v : used_vars) {
            double const lb = b[v].lb();
            double const ub = b[v].ub();
            add_neg_intv(lb, v, ub);  // !(lb <= v <= ub)
        }
    }
    // Add blocking clause B1 => B2, but generalize it using used_vars
    void add_generalized_blocking_box(box const & b1, box const & b2, unordered_set<Enode *> const & used_vars) {
        // B1 => B2
        //
        // /\ I1_j => /\ I2_i
        //  j          i
        //
        // !(/\ I1_j) \/ /\ I2_i
        //   j            i
        //
        //     \/ !I1_j \/ I2_1
        //     j
        // /\      ...
        //     \/ !I1_j \/ I2_n
        //     j
        for (Enode * v : used_vars) {
            double const lb = b[v].lb();
            double const ub = b[v].ub();
            add_neg_intv(lb, v, ub);  // !(lb <= v <= ub)
        }

        // TODO(soonhok): need to add ordering between each interval I1_i and I2_i
    }

    // Add B => l1
    void add_imply(box const & b, int l1) {
        add_imply(b, l1, 0);
    }
    // Add B => l1 \/ l2
    void add_imply(box const & b, int l1, int l2) {
        // Add !B
        for (Enode * v : b.get_vars()) {
            double const lb = b[v].lb();
            double const ub = b[v].ub();
            //    !(lb <= v /\ v <= ub)
            // => !(lb <= v) \/ !(v <= ub)
            // =>  (v <= lb) \/ !(v <= ub)
            picosat_add(m_psat, m_store.add(v, lb));   //  (v <= lb)
            picosat_add(m_psat, -m_store.add(v, ub));  // !(v <= ub)
        }
        picosat_add(m_psat, l1);
        if (l2) {
            picosat_add(m_psat, l2);
        }
        picosat_add(m_psat, 0);
    }

    // Add B => (B[v].lb <= v <= m) xor (m <= v <= B[v].ub)
    void add_branching(box const & b, Enode * v, double m) {
        // TODO(soonhok): only do this if v m is not in the store
        double const lb = b[v].lb();
        double const ub = b[v].ub();
        // 1. In CNF Form
        //  1.1. B => (B[v].lb <= v <= m) \/ (m <= v <= B[v].ub)
        //  1.2. B => !(B[v].lb <= v <= m) /\ (m <= v <= B[v].ub)

        //  1.1 part:
        //   B => (l <= v) \/ (m <= v) --> B => (l <= v) --> B <= !(v <= l)
        add_imply(b, -m.store.add(v, lb));
        //   B => (l <= v) \/ (v <= u) --> B => !(v <= l) \/ (v <= u)
        add_imply(b, m.store.add(v, ), m.store.add(v, ));
        //   B => (v <= m) \/ (m <= v) --> B => True
        //   --> Nothing to ADD!
        //   B => (v <= m) \/ (v <= u) --> B => (v <= m) \/ (v <= u) --> B => (v <= u)
        add_imply(b, m.store.add(v, u));

        //  1.2 part:
        //   B => !(l <= v /\ v <= m) \/ !(m <= v /\ v <= u)
        //   B => !(l <= v) \/ !(v <= m) \/ !(m <= v) /\ !(v <= u)
        //   B =>  (v <= l) \/ !(v <= m) \/  (v <= m) /\ !(v <= u)
        //   B => True (because of "!(v <= m) \/ (v <= m)")
        //   --> Nothing to ADD!

        // 2. Need to provide ordering among (v <= lb), (v <= m), (v <= ub)
        // (v <= lb) => (v <= m) --> !(v <= lb) \/ (v <= m)
        picosat_add(m_psat, -m_store.add(v, lb));
        picosat_add(m_psat, m_store.add(v, m));
        picosat_add(m_psat, 0);
        // (v <= m) => (v <= ub) --> !(v <= m) \/ (v <= ub)
        picosat_add(m_psat, -m_store.add(v, m));
        picosat_add(m_psat, m_store.add(v, ub));
        picosat_add(m_psat, 0);

        // Debug Print
        DREAL_LOG_INFO << "Branching on: " << v << "\t"
                       << "[" << lb << ", " << m  << "], "
                       << "[" << m  << ", " << ub << "]";
    }

    int check_sat() {
        return picosat_sat(m_psat, -1);
    }
    // Precondition: check_sat() == PICOSAT_SATISFIABLE
    // Reduce the given box b into a smaller box using SAT model
    box reduce_using_model(box b) {
        for (int i = 1; i <= m_store.get_num_vars(); i++) {
            int const r = picosat_deref_partial(m_psat, i);
            if (r == 0) { continue;  /* UNKNOWN */ }
            // pred := v <= bound
            tuple<Enode*, double> pred = m_store.lookup(i);
            Enode * v = get<0>(pred);
            double const bound = get<1>(pred);
            if (r == 1) {
                // b_i = True
                // ==> (v <= bound)
                b[v] &= ibex::Interval(b[v].lb(), bound);
            } else {
                // b_i = False
                // ==> (bound <= v)
                assert(r == -1);
                b[v] &= ibex::Interval(bound, b[v].ub())
            }
        }
    }
}

box sat_icp::solve(box b, contractor & ctc, SMTConfig & config) {
    thread_local static std::unordered_set<std::shared_ptr<constraint>> used_constraints;
    used_constraints.clear();

    // Step 1. Initialize SAT Solver
    picosat_wrapper pw;
    //  - Add initial box (unit) clauses
    pw.add_box(b);
    box const initial_box(b);

    // Step 2. Main Part
    while (true) {
        DREAL_LOG_INFO << "\n\n\n\n";
        // Ask SAT Solver for a satisfying assignment
        // DREAL_LOG_INFO << "===============================";
        // picosat_print(psat, stderr);
        // DREAL_LOG_INFO << "===============================";
        // store.debug_print();
        // DREAL_LOG_INFO << "===============================";

        int ret = pw.check_sat();
        if (ret == PICOSAT_SATISFIABLE) {
            DREAL_LOG_INFO << "SAT solver found a satisfying Boolean assignment";
            // Case 1: SAT solver found a satisfying Boolean assignment
            // 1.1. Concretize the satisfying Boolean assignment into a conjunction of constraints
            // Check each Boolean Variable and if partially assigned to be true, shrink the interval
            DREAL_LOG_INFO << "store.get_num_vars() = " << store.get_num_vars();
            b = pw.reduce_using_model(initial_box);
            DREAL_LOG_INFO << "Current Box = " << "\n"
                            << b;
            // 1.2. Apply pruning operators with box B until it reaches a fixed point B'
            box old_b(b);
            try {
                ctc.prune(b, config);
                auto const this_used_constraints = ctc.used_constraints();
                used_constraints.insert(this_used_constraints.begin(), this_used_constraints.end());
            } catch (contractor_exception & e) {
                // Do nothing
            }
            if (config.nra_use_stat) { config.nra_stat.increase_prune(); }
            if (b.is_empty()) {
                DREAL_LOG_INFO << "After Pruning, it became an empty set.";
                // Case i: Pruning returns an empty box.
                // 1.2.i.1. Collect Used Variables
                unordered_set<Enode *> used_vars;
                used_vars.clear();
                for (auto used_ctr : ctc.used_constraints()) {
                    auto this_used_vars = used_ctr->get_vars();
                    used_vars.insert(this_used_vars.begin(), this_used_vars.end());
                }
                // 1.2.i.2. Add a blocking clause
                pw.add_generalized_blocking_box(old_b, used_vars);
            } else {
                DREAL_LOG_INFO << "After Pruning, it became a non-empty set.";
                DREAL_LOG_INFO << b;
                // Case ii: Pruning returns a non-empty box b
                if (b.max_diam() < config.nra_precision) {
                    DREAL_LOG_INFO << "Box is small enough to stop.";
                    // Box is small enough to stop => delta-SAT
                    break;
                } else {
                    DREAL_LOG_INFO << "Box is big: width = " << b.max_diam();
                    // Box is big, and needs to be branched. We also need to learn a clause from this pruning
                    // ii.1. Pick a branching variable
                    if (config.nra_use_stat) { config.nra_stat.increase_branch(); }

                    // TODO(soonhok): b.bisect is an overkill here
                    // since it returns two boxes which are not
                    // necessary
                    auto const bisect_result = b.bisect(config.nra_precision);
                    Enode * br_var = b.get_vars()[get<0>(bisect_result)];
                    // ii.2. Pick a branching point
                    double const br_point = b[br_var].mid();
                    pw.add_branching(b, br_var, br_point);

                    // // ii.6. Add ordering between old and new interval in each dimension.
                    // for (unsigned i = 0; i < b.size(); i++) {
                    //     auto const & old_intv = old_b[i];
                    //     auto const & intv = b[i];
                    //     if (old_intv != intv) {
                    //         int const old_var = store.add(vars[i], old_intv.lb(), old_intv.ub());
                    //         int const var     = store.add(vars[i],     intv.lb(),     intv.ub());
                    //         //  var => old_var
                    //         // !var \/ old_var
                    //         picosat_add(psat,    -var);
                    //         picosat_add(psat, old_var);
                    //         picosat_add(psat,       0);
                    //         DREAL_LOG_INFO << "Add to SAT Solver: " << -var << " " << old_var << " 0";
                    //     }
                    // }
                    // // ii.7. Add learned clause
                    // //
                    // // Old Box => New Box
                    // // !(/\ Old_Intv_j) \/ /\ New_Intv_i
                    // //    j                i
                    // //
                    // // In CNF:
                    // //
                    // //     \/ Old_Intv_j \/ New_Intv_1
                    // //      j
                    // // /\  ...
                    // //     \/ Old_Intv_j \/ New_Intv_n
                    // //      j
                    // //
                    // for (unsigned i = 0; i < b.size(); i++) {
                    //     for (int old_var : vars_of_old_B) { picosat_add(psat, -old_var); }
                    //     auto const & intv = b[i];
                    //     int const var     = store.add(vars[i],     intv.lb(),     intv.ub());
                    //     picosat_add(psat, var);
                    //     picosat_add(psat,   0);
                    // }
                }
            }
        } else if (ret == PICOSAT_UNSATISFIABLE) {
            DREAL_LOG_INFO << "SAT solver failed to find a satisfying Boolean assignment";
            // Case 2: SAT solver concludes UNSAT. Return UNSAT.
            b.set_empty();
            break;
        } else {
            assert(ret == PICOSAT_UNKNOWN);
            DREAL_LOG_INFO << "SAT Solver failed.";
        }
    }
    picosat_reset(psat);
    ctc.set_used_constraints(used_constraints);
    return b;
}
}  // namespace dreal
