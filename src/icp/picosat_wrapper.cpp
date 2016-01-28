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

#include <unordered_set>
#include <unordered_map>
#include <tuple>
#include "icp/picosat_wrapper.h"

using std::unordered_set;
using std::tuple;
using std::get;

namespace dreal {
    picosat_wrapper::picosat_wrapper() {
        m_psat = picosat_init();
        picosat_save_original_clauses(m_psat);  // to use picosat_deref_partial function
    }
    picosat_wrapper::~picosat_wrapper() {
        picosat_reset(m_psat);
    }

    // Add: v <= bound
    void picosat_wrapper::add_le(Enode * v, double const bound) {
        picosat_add(m_psat, m_store.add(v, bound));
        picosat_add(m_psat, 0);
        DREAL_LOG_FATAL << "PICOSAT WRAPPER: ADD - (" << v << " <= " << bound << ")";
    }
    // Add: bound <= v
    void picosat_wrapper::add_le(double const bound, Enode* v) {
        // Note: Strictly, !(x<=v) is (v < x) but we consider this as
        // (v <= x). This should be OK since we use numerical alg
        // anyway.
        picosat_add(m_psat, -m_store.add(v, bound));
        picosat_add(m_psat, 0);
        DREAL_LOG_FATAL << "PICOSAT WRAPPER: ADD - (" << bound << " <= " << v << ")";
    }
    // Add: lb <= v <= ub
    void picosat_wrapper::add_intv(double const lb, Enode* v, double const ub) {
        add_le(lb, v); add_le(v, ub);
    }
    // Add: !(lb <= v /\ v <= ub)
    void picosat_wrapper::add_neg_intv(double const lb, Enode* v, double const ub) {
        // That is, v <= lb \/ !(v <= ub)
        picosat_add(m_psat, m_store.add(v, lb));
        picosat_add(m_psat, -m_store.add(v, ub));
        picosat_add(m_psat, 0);
        DREAL_LOG_FATAL << "PICOSAT WRAPPER: ADD - (" << v << " <= " << lb << ") \\/ "
                        << "(" << ub << " <= " << v << ")";
    }
    void picosat_wrapper::add_box(box const & b) {
        auto const & vars = b.get_vars();
        for (unsigned i = 0; i < b.size(); ++i) {
            Enode * v = vars[i];
            double const lb = b[v].lb();
            double const ub = b[v].ub();
            add_intv(lb, v, ub);  // lb <= v <= ub
        }
    }

    // Add blocking clause ¬B, but generalize it using used_vars
    void picosat_wrapper::add_generalized_blocking_box(box const & b, unordered_set<Enode *> const & used_vars) {
        for (Enode * v : used_vars) {
            double const lb = b[v].lb();
            double const ub = b[v].ub();
            add_neg_intv(lb, v, ub);  // !(lb <= v <= ub)
        }
    }
    // Add blocking clause B1 => B2, but generalize it using used_vars
    void picosat_wrapper::add_generalized_blocking_box(box const & b1, box const & b2, unordered_set<Enode *> const & used_vars) {
        DREAL_LOG_FATAL << "picosat_wrapper::add_generalized_blocking_box";
        DREAL_LOG_FATAL << "box1 = " << b1;
        DREAL_LOG_FATAL << "box2 = " << b2;
        if (used_vars.empty()) {
            DREAL_LOG_FATAL << "used var = empty!";
        } else {
            for (Enode* v : used_vars) {
                DREAL_LOG_FATAL << "used var = " << v;
            }
        }
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
        for (Enode * v_i : used_vars) {
            // \/ !I1_j \/ (b2[v_i].lb <= v_i)
            //  j
            for (Enode * v_j : used_vars) {
                double const i1_lb = b1[v_j].lb();
                double const i1_ub = b1[v_j].ub();
                picosat_add(m_psat, m_store.add(v_j, i1_lb));
                picosat_add(m_psat, -m_store.add(v_j, i1_ub));
            }
            double const i2_lb = b2[v_i].lb();
            picosat_add(m_psat, -m_store.add(v_i, i2_lb));
            picosat_add(m_psat, 0);
            // \/ !I1_j \/ (v_i <= b2[v_i].ub)
            //  j
            for (Enode * v_j : used_vars) {
                double const i1_lb = b1[v_j].lb();
                double const i1_ub = b1[v_j].ub();
                picosat_add(m_psat, m_store.add(v_j, i1_lb));
                picosat_add(m_psat, -m_store.add(v_j, i1_ub));
            }
            double const i2_ub = b2[v_i].ub();
            picosat_add(m_psat, m_store.add(v_i, i2_ub));
            picosat_add(m_psat, 0);
        }

        // Add ordering between each interval I1_i and I2_i
        for (Enode * v : used_vars) {
            double const i1_lb = b1[v].lb();
            double const i1_ub = b1[v].ub();
            double const i2_lb = b2[v].lb();
            double const i2_ub = b2[v].ub();
            // [               i1                  ]
            //           [     i2      ]
            // i1_lb   i2_lb         i2_ub       i1_ub
            //
            // (v <= i1_lb) => (v <= i2_lb)
            if (i1_lb != i2_lb) {
                picosat_add(m_psat, -m_store.add(v, i1_lb));
                picosat_add(m_psat,  m_store.add(v, i2_lb));
                picosat_add(m_psat,  0);
            }
            // (v <= i2_ub) => (v <= i1_ub)
            if (i2_ub != i1_ub) {
                picosat_add(m_psat, -m_store.add(v, i2_ub));
                picosat_add(m_psat,  m_store.add(v, i1_ub));
                picosat_add(m_psat,  0);
            }
        }
    }

    // Add B => l1 \/ l2 \/ l3 \/ l4
    void picosat_wrapper::add_imply(box const & b, int const l1, int const l2, int const l3, int const l4) {
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
            if (l3) {
                picosat_add(m_psat, l3);
                if (l4) {
                    picosat_add(m_psat, l4);
                }
            }
        }
        picosat_add(m_psat, 0);
        DREAL_LOG_FATAL << "PICOSAT WRAPPER: ADD - !B => " << l1 << " " << l2;
    }

    // Add B => (B[v].lb <= v <= m) xor (m <= v <= B[v].ub)
    void picosat_wrapper::add_branching(box const & b, Enode * v, double const m) {
        // TODO(soonhok): only do this if v m is not in the store
        double const lb = b[v].lb();
        double const ub = b[v].ub();
        // 1. In CNF Form
        //  1.1. B => (B[v].lb <= v <= m) \/ (m <= v <= B[v].ub)
        //  1.2. B => !(B[v].lb <= v <= m) /\ (m <= v <= B[v].ub)

        //  1.1 part:
        //   B => (l <= v) \/ (m <= v) --> B => (l <= v) --> B <= !(v <= l)
        add_imply(b, -m_store.add(v, lb), -m_store.add(v, m));
        //   B => (l <= v) \/ (v <= u) --> B => !(v <= l) \/ (v <= u)
        add_imply(b, -m_store.add(v, lb), m_store.add(v, ub));
        //   B => (v <= m) \/ (m <= v) --> B => True
        add_imply(b,  m_store.add(v, m), -m_store.add(v, m));
        //   B => (v <= m) \/ (v <= u) --> B => (v <= m) \/ (v <= u) --> B => (v <= u)
        add_imply(b,  m_store.add(v, m),  m_store.add(v, ub));

        //  1.2 part:
        //   B => !(l <= v /\ v <= m) \/ !(m <= v /\ v <= u)
        //   B => !(l <= v) \/ !(v <= m) \/ !(m <= v) /\ !(v <= u)
        //   B =>  (v <= l) \/ !(v <= m) \/  (v <= m) /\ !(v <= u)
        add_imply(b, m_store.add(v, lb), -m_store.add(v, m), m_store.add(v, m), -m_store.add(v, ub));

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
        DREAL_LOG_FATAL << "Branching on: " << v << "\t"
                       << "[" << lb << ", " << m  << "], "
                       << "[" << m  << ", " << ub << "]";
    }

    int picosat_wrapper::check_sat() {
        return picosat_sat(m_psat, -1);
    }
    // Precondition: check_sat() == PICOSAT_SATISFIABLE
    // Reduce the given box b into a smaller box using SAT model
    box picosat_wrapper::reduce_using_model(box b) const {
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
                DREAL_LOG_FATAL << "b[" << v << "] : " << b[v] << " /\\ " << ibex::Interval(b[v].lb(), bound);
                b[v] &= ibex::Interval(b[v].lb(), bound);
                DREAL_LOG_FATAL << " = " << b[v];
            } else {
                // b_i = False
                // ==> (bound <= v)
                assert(r == -1);
                DREAL_LOG_FATAL << "b[" << v << "] : " << b[v] << " /\\ " << ibex::Interval(bound, b[v].ub());
                b[v] &= ibex::Interval(bound, b[v].ub());
                DREAL_LOG_FATAL << " = " << b[v];
            }
        }
        return b;
    }

    void picosat_wrapper::debug_print() const {
        DREAL_LOG_FATAL << "======================";
        for (int i = 1; i <= m_store.get_num_vars(); i++) {
            tuple<Enode*, double> pred = m_store.lookup(i);
            Enode * v = get<0>(pred);
            double const bound = get<1>(pred);
            DREAL_LOG_FATAL << "b" << i << " := "
                            << "(" << v << " <= " << bound << ")";
        }
        DREAL_LOG_FATAL << "~~~~~~~~~~~~~~~~~~~~~~";
        picosat_print(m_psat, stderr);
        DREAL_LOG_FATAL << "======================";
    }

}  // namespace dreal
