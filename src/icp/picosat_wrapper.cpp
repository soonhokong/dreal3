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
using std::cerr;
using std::endl;

namespace dreal {
    picosat_wrapper::picosat_wrapper() {
        m_psat = picosat_init();
        picosat_save_original_clauses(m_psat);  // to use picosat_deref_partial function
    }
    picosat_wrapper::~picosat_wrapper() {
        picosat_reset(m_psat);
    }

    // Given a variable `v` and two constants `l` and `u` where l <= u holds.
    void picosat_wrapper::add_ordering(box const & b, Enode *v, double const l, double const u) {
        if (l < u) {
            // Linear ordering on "<=": B =>  (v <= lb) => (v <= ub)
            //                           --> !(v <= lb) \/ (v <= ub)
            add_imply(b, -m_store.add(v, l, true), m_store.add(v, u, true));
            // Linear ordering on ">=":       (v >= ub) => (v >= lb)
            //                           --> !(v >= ub) \/ (v >= lb)
            add_imply(b, -m_store.add(v, u, false), m_store.add(v, l, false));
            //                          B =>  (v <= lb) => !(v >= ub)
            //                           --> !(v <= lb) \/ !(v >= ub)
            add_imply(b, -m_store.add(v, l, true), -m_store.add(v, u, false));
            //                          B =>  (v >= ub) => !(v <= lb)
            //                           --> !(v >= ub) \/ !(v <= lb)
            add_imply(b, -m_store.add(v, u, false), -m_store.add(v, l, true));
        }

        // // Compatibility between lower and upper bounds: B =>  (v <= lb) => !(v >= lb)
        // //                                                --> !(v <= lb) \/ !(v >= lb)
        // add_imply(b, -m_store.add(v, l, true), -m_store.add(v, l, false));
        // // Compatibility between lower and upper bounds: B =>  (v <= ub) => !(v >= ub)
        // //                                                --> !(v <= ub) \/ !(v >= ub)
        // add_imply(b, -m_store.add(v, u, true), -m_store.add(v, u, false));

        // (v >= l) \/ (v <= l)
        add_imply(b, m_store.add(v, l, false), m_store.add(v, l, true));
        // (v >= u) \/ (v <= u)
        add_imply(b, m_store.add(v, u, false), m_store.add(v, u, true));
    }

    // Add: v <= bound
    void picosat_wrapper::add_le(Enode * v, double const bound) {
        picosat_add(m_psat, m_store.add(v, bound, true));
        picosat_add(m_psat, 0);
        DREAL_LOG_INFO << "PICOSAT WRAPPER: ADD - (" << v << " <= " << bound << ")";
    }

    // Add: bound <= v  <-->  (v >= bound)
    void picosat_wrapper::add_le(double const bound, Enode* v) {
        add_ge(v, bound);
    }

    // Add: v >= bound
    void picosat_wrapper::add_ge(Enode * v, double const bound) {
        picosat_add(m_psat, m_store.add(v, bound, false));
        picosat_add(m_psat, 0);
        DREAL_LOG_INFO << "PICOSAT WRAPPER: ADD - (" << v << " >= " << bound << ")";
    }

    // Add: bound >= v  <-->  (v <= bound)
    void picosat_wrapper::add_ge(double const bound, Enode* v) {
        add_le(v, bound);
    }

    // // Add: lb <= v <= ub
    // void picosat_wrapper::add_intv(double const lb, Enode* v, double const ub) {
    // }

    void picosat_wrapper::add_box(box const & b) {
        auto const & vars = b.get_vars();
        for (unsigned i = 0; i < b.size(); ++i) {
            Enode * v = vars[i];
            double const lb = b[v].lb();
            double const ub = b[v].ub();
            // lb <= v <= ub
            add_le(lb, v);
            add_le(v, ub);
            add_ordering(b, v, lb, ub);
        }
    }

    // Add blocking clause ¬B, but generalize it using used_vars
    void picosat_wrapper::add_generalized_blocking_box(box const & b, unordered_set<Enode *> const & used_vars) {
        for (Enode * v : used_vars) {
            double const lb = b[v].lb();
            double const ub = b[v].ub();
            add_ordering(b, v, lb, ub);
        }
        for (Enode * v : used_vars) {
            double const lb = b[v].lb();
            double const ub = b[v].ub();
            // !(lb <= v <= ub) --> !(lb <= v  /\   v <= ub)
            //                  --> !(lb <= v) \/ !(v <= ub)
            //                  --> !(v >= lb) \/ !(v <= ub)
            picosat_add(m_psat, -m_store.add(v, lb, false));
            picosat_add(m_psat, -m_store.add(v, ub, true));
        }
        picosat_add(m_psat, 0);
    }

    // Add blocking clause B1 => B2, but generalize it using used_vars
    void picosat_wrapper::add_generalized_blocking_box(box const & b1, box const & b2, unordered_set<Enode *> const & used_vars) {
        DREAL_LOG_INFO << "picosat_wrapper::add_generalized_blocking_box";
        DREAL_LOG_INFO << "box1 = " << b1;
        DREAL_LOG_INFO << "box2 = " << b2;
        if (used_vars.empty()) {
            DREAL_LOG_INFO << "used var = empty!";
        } else {
            for (Enode* v : used_vars) {
                DREAL_LOG_INFO << "used var = " << v;
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
            //     \/ !I1_j \/ (b2[v_i].lb <= v_i)
            //     j
            // --> \/ !I1_j \/ (v_i >= b2[v_i].lb)
            //     j
            for (Enode * v_j : used_vars) {
                double const i1_lb = b1[v_j].lb();
                double const i1_ub = b1[v_j].ub();
                //     !((i1.lb <= v_j) /\  (v_j <= i1.ub))
                // -->  !(i1.lb <= v_j) \/ !(v_j <= i1.ub)
                // -->  !(v_j >= i1.lb) \/ !(v_j <= i1.ub)
                picosat_add(m_psat, -m_store.add(v_j, i1_lb, false));
                picosat_add(m_psat, -m_store.add(v_j, i1_ub, true));
            }
            double const i2_lb = b2[v_i].lb();
            picosat_add(m_psat, m_store.add(v_i, i2_lb, false));  //i2_lb <= v_i
            picosat_add(m_psat, 0);

            // \/ !I1_j \/ (v_i <= b2[v_i].ub)
            //  j
            for (Enode * v_j : used_vars) {
                double const i1_lb = b1[v_j].lb();
                double const i1_ub = b1[v_j].ub();
                //     !((i1.lb <= v_j) /\  (v_j <= i1.ub))
                // -->  !(i1.lb <= v_j) \/ !(v_j <= i1.ub)
                // -->  !(v_j >= i1.lb) \/ !(v_j <= i1.ub)
                picosat_add(m_psat, -m_store.add(v_j, i1_lb, false));
                picosat_add(m_psat, -m_store.add(v_j, i1_ub, true));
            }
            double const i2_ub = b2[v_i].ub();
            picosat_add(m_psat, m_store.add(v_i, i2_ub, true));  // v_i <= i2_ub
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
            add_ordering(b1, v, i1_lb, i2_lb);
            add_ordering(b1, v, i2_lb, i2_ub);
            add_ordering(b1, v, i2_ub, i1_ub);
        }
    }

    // Add B => l1 \/ l2 \/ l3 \/ l4
    void picosat_wrapper::add_imply(box const & b, int const l1, int const l2, int const l3, int const l4) {
        // Add !B
        // std::cerr << "PICOSAT WRAPPER: Add Imply: ";
        for (Enode * v : b.get_vars()) {
            double const lb = b[v].lb();
            double const ub = b[v].ub();
            //     !((lb <= v) /\  (v <= ub))
            // --> !(lb <= v)  \/ !(v <= ub)
            // --> !(v >= lb)  \/ !(v <= ub)
            picosat_add(m_psat, -m_store.add(v, lb, false));
            picosat_add(m_psat, -m_store.add(v, ub, true));
            // std::cerr << m_store.add(v, lb) << " " << -m_store.add(v, ub) << " ";
        }
        picosat_add(m_psat, l1);
        // std::cerr << l1 << " ";
        if (l2) {
            picosat_add(m_psat, l2);
            // std::cerr << l2 << " ";
            if (l3) {
                picosat_add(m_psat, l3);
                // std::cerr << l3 << " ";
                if (l4) {
                    picosat_add(m_psat, l4);
                    // std::cerr << l4 << " ";
                }
            }
        }
        picosat_add(m_psat, 0);
        // std::cerr << 0 << endl;
        // DREAL_LOG_INFO << "PICOSAT WRAPPER: ADD - !B => "
        //                 << l1 << " " << l2 << " " << l3 << " " << l4;
    }

    // Add B => (B[v].lb <= v < m) xor (m <= v <= B[v].ub)
    void picosat_wrapper::add_branching(box const & b, Enode * v, double const m) {
        // TODO(soonhok): only do this if v m is not in the store
        double const lb = b[v].lb();
        double const ub = b[v].ub();
        DREAL_LOG_INFO << "ADD_BRANCHING on "
                        << v << "[" << lb << ", " << m << ", " << ub << "]\n"
                        << b;
        // 1. In CNF Form
        //  1.1. B => (B[v].lb <= v <= m) \/ (m <= v <= B[v].ub)
        //   B => (l <= v) \/ (m <= v) --> B => (v >= l) \/ (v >= m)
        add_imply(b, m_store.add(v, lb, false), m_store.add(v, m, false));
        //   B => (l <= v) \/ (v <= u) --> B => (v >= l) \/ (v <= u)
        add_imply(b, m_store.add(v, lb, false), m_store.add(v, ub, true));
        //   B => (v <= m) \/ (m <= v) --> B => (v <= m) \/ (v >= m)
        add_imply(b, m_store.add(v, m, true),   m_store.add(v, m, false));  // TODO(soonhok): this is a duplicate from the ordering
        //   B => (v <= m) \/ (v <= u) --> B => (v <= m) \/ (v <= u)
        add_imply(b, m_store.add(v, m, true),   m_store.add(v, ub, true));

        //  1.2. B => !(B[v].lb <= v <= m) /\ !(m <= v <= B[v].ub)
        //   B => !(l <= v /\ v <= m) \/ !(m <= v /\ v <= u)
        //   B => !(l <= v) \/ !(v <= m) \/ !(m <= v) /\ !(v <= u)
        //   B => !(v >= l) \/ !(v <= m) \/ !(v >= m) /\ !(v <= u)
        add_imply(b, -m_store.add(v, lb, false)
                   , -m_store.add(v, m,  true)
                   , -m_store.add(v, m,  false)
                   , -m_store.add(v, ub, true));

        // 2. Need to provide ordering among lb, m, and ub.
        add_ordering(b, v, lb, m);
        add_ordering(b, v, m, ub);

        // Debug Print
        DREAL_LOG_WARNING << "Branching on: " << v << "\t"
                        << "[" << lb << ", " << m  << "], "
                        << "[" << m  << ", " << ub << "]";
    }

    int picosat_wrapper::check_sat() {
        int const ret = picosat_sat(m_psat, -1);
        // picosat_simplify(m_psat);
        return ret;
    }
    // Precondition: check_sat() == PICOSAT_SATISFIABLE
    // Reduce the given box b into a smaller box using SAT model
    box picosat_wrapper::reduce_using_model(box b) const {
        // TODO(soonhok): this can be a bottleneck. Consider optimization.
        for (int i = 1; i <= m_store.get_num_vars(); i++) {
            int const r = picosat_deref_partial(m_psat, i);
            // pred := v <= bound
            tuple<Enode*, double, bool> pred = m_store.lookup(i);
            Enode * v = get<0>(pred);
            double const bound = get<1>(pred);
            bool const le = get<2>(pred);

            DREAL_LOG_WARNING << "b" << i << "\t"
                            << (r == 1 ? "+" : (r == 0 ? "0 " : "! "))
                            << v
                            << (le ? " <= " : " >= ")
                            << bound;

            if (r == 0) { continue;  /* UNKNOWN */ }

            if (r == 1 && le) {
                // (v <= bound)
                DREAL_LOG_WARNING << "b[" << v << "] : "
                                << b[v] << " /\\ " << ibex::Interval(b[v].lb(), bound)
                                << " [" << b[v].lb() << ", " << bound << "]";
                b[v] &= ibex::Interval(b[v].lb(), bound);
                DREAL_LOG_WARNING << " = " << b[v];
            } else if (r == 1 && !le) {
                // (v >= bound)
                DREAL_LOG_WARNING << "b[" << v << "] : "
                                << b[v] << " /\\ " << ibex::Interval(bound, b[v].ub())
                                << " [" << bound << ", " << b[v].ub() << "]";
                b[v] &= ibex::Interval(bound, b[v].ub());
                DREAL_LOG_WARNING << " = " << b[v];
            }//  else if (r == -1 && le) {
            //     // !(v <= bound) --> (v > bound)
            //     if (bound == b[v].ub()) {
            //         b[v].set_empty();
            //     } else {
            //         DREAL_LOG_WARNING << "b[" << v << "] : "
            //                         << b[v] << " /\\ " << ibex::Interval(bound, b[v].ub())
            //                         << " (" << bound << ", " << b[v].ub() << "]";
            //         b[v] &= ibex::Interval(bound, b[v].ub());
            //         DREAL_LOG_WARNING << " = " << b[v];
            //     }
            // } else if (r == -1 && !le) {
            //     // !(v >= bound) --> (v < bound)
            //     if (bound == b[v].lb()) {
            //         b[v].set_empty();
            //     } else {
            //         DREAL_LOG_WARNING << "b[" << v << "] : "
            //                         << b[v] << " /\\ " << ibex::Interval(b[v].lb(), bound)
            //                         << " [" << b[v].lb() << ", " << bound << ")";
            //         b[v] &= ibex::Interval(b[v].lb(), bound);
            //         DREAL_LOG_WARNING << " = " << b[v];
            //     }
            // }
            if (b[v].is_empty()) {
                DREAL_LOG_WARNING << "SOMETHING IS WRONG, WE GOT AN EMPTY INTERVAL HERE";
                abort();
                b.set_empty();
                break;
            }
        }
        return b;
    }

    void picosat_wrapper::debug_print() const {
        DREAL_LOG_WARNING << "======================";
        for (int i = 1; i <= m_store.get_num_vars(); i++) {
            tuple<Enode*, double, bool> pred = m_store.lookup(i);
            Enode * v = get<0>(pred);
            double const bound = get<1>(pred);
            bool const le = get<2>(pred);
            DREAL_LOG_WARNING << "b" << i << " := "
                            << "(" << v
                           << (le ? " <= " : " >= ")
                            << std::setprecision(16) << bound << ")";
        }
        DREAL_LOG_WARNING << "~~~~~~~~~~~~~~~~~~~~~~";
        picosat_print(m_psat, stderr);
        DREAL_LOG_WARNING << "======================";
    }

}  // namespace dreal
