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

#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <tuple>
#include "icp/picosat_wrapper.h"

using std::unordered_set;
using std::tuple;
using std::get;
using std::cerr;
using std::endl;
using std::vector;

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
        assert(l <= u);
        if (l < u) {
            auto const & var_vec = b.get_vars();
            unordered_set<Enode *> const used_vars(var_vec.begin(), var_vec.end());
            return add_ordering(b, used_vars, v, l, u);
        }
    }
    void picosat_wrapper::add_ordering(box const & b, unordered_set<Enode *> const & used_vars, Enode *v, double const l, double const u) {
        assert(l <= u);
        if (l < u) {
            // B =>  (v <= lb) =>  (v <= ub)
            //  --> !(v <= lb) \/  (v <= ub)
            add_imply(b, used_vars, -m_store.add(v, l, true),   m_store.add(v, u, true));
            // B =>  (v >= ub) =>  (v >= lb)
            //  --> !(v >= ub) \/  (v >= lb)
            //  -->  (v >= lb) \/ !(v >= ub)  -- comm
            add_imply(b, used_vars, m_store.add(v, l, false), -m_store.add(v, u, false));
            // B =>  (v <= lb) =>  !(v >= ub)
            //  --> !(v <= lb) \/  !(v >= ub)
            add_imply(b, used_vars, -m_store.add(v, l, true),  -m_store.add(v, u, false));
        }
    }

    void picosat_wrapper::add_axiom(box const & b, Enode * v, double const c) {
        // (v <= c) \/ (v >= c)
        auto const & var_vec = b.get_vars();
        unordered_set<Enode *> const used_vars(var_vec.begin(), var_vec.end());
        return add_axiom(b, used_vars, v, c);
    }

    void picosat_wrapper::add_axiom(box const & b, unordered_set<Enode*> const & used_vars, Enode * v, double const c) {
        // (v <= c) \/ (v >= c)
        add_imply(b, used_vars, m_store.add(v, c, true), m_store.add(v, c, false));
    }

    // Add: v <= bound
    void picosat_wrapper::add_le(Enode * v, double const bound) {
        picosat_add(m_psat, m_store.add(v, bound, true));
        picosat_add(m_psat, 0);
        DREAL_LOG_WARNING << "PICOSAT WRAPPER: ADD - (" << v << " <= " << bound << ")";
    }

    // Add: bound <= v  <-->  (v >= bound)
    void picosat_wrapper::add_le(double const bound, Enode* v) {
        add_ge(v, bound);
    }

    // Add: v >= bound
    void picosat_wrapper::add_ge(Enode * v, double const bound) {
        picosat_add(m_psat, m_store.add(v, bound, false));
        picosat_add(m_psat, 0);
        DREAL_LOG_WARNING << "PICOSAT WRAPPER: ADD - (" << v << " >= " << bound << ")";
    }

    // Add: bound >= v  <-->  (v <= bound)
    void picosat_wrapper::add_ge(double const bound, Enode* v) {
        add_le(v, bound);
    }

    // // Add: lb <= v <= ub
    // void picosat_wrapper::add_intv(double const lb, Enode* v, double const ub) {
    // }

    void picosat_wrapper::add_box(box const & b) {
        for (Enode * v : b.get_vars()) {
            double const lb = b[v].lb();
            double const ub = b[v].ub();
            // lb <= v <= ub
            add_le(lb, v);
            add_le(v, ub);
        }
        for (Enode * v : b.get_vars()) {
            double const lb = b[v].lb();
            double const ub = b[v].ub();
            assert(lb <= ub); add_ordering(b, v, lb, ub);
            add_axiom(b, v, lb);
            if (lb < ub) {
                add_axiom(b, v, ub);
            }
        }
    }

    // Add blocking clause ¬B, but generalize it using used_vars
    void picosat_wrapper::add_generalized_blocking_box(box const & b, unordered_set<Enode *> const & vars) {
        // assert(used_vars.size() > 0);
        // assert(vars.size() > 0);
        unordered_set<Enode*> used_vars(b1.get_vars().begin(), b1.get_vars().end());
        DREAL_LOG_WARNING << "The following BOX is UNSAT and blocked by add_generalized_blocking_box function";
        DREAL_LOG_WARNING << b;
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

        // // TODO(soonhok): Safety Net
        // block_current_assignment();
    }

    // Add blocking clause B1 => B2, but generalize it using used_vars
    void picosat_wrapper::add_generalized_blocking_box(box const & b1, box const & b2, unordered_set<Enode *> const & vars) {
        // assert(used_vars.size() > 0);
        assert(vars.size() > 0);
        unordered_set<Enode*> used_vars(b1.get_vars().begin(), b1.get_vars().end());
        assert(b1 != b2);
        assert(b1.is_superset(b2));
        DREAL_LOG_WARNING << "picosat_wrapper::add_generalized_blocking_box";
        DREAL_LOG_WARNING << "box1 = " << b1;
        DREAL_LOG_WARNING << "box2 = " << b2;
        // B1 => B2
        //
        //   /\ I1_j => /\ I2_i  --- (1)
        //    j          i
        //
        // !(/\ I1_j) \/ /\ I2_i --- (2)
        //   j            i
        //
        //     \/ !I1_j \/ I2_1
        //     j
        // /\      ...           --- (3)
        //     \/ !I1_j \/ I2_n
        //     j
        vector<int> b1_lits;
        for (Enode * v1 : used_vars) {
            double const i1_lb = b1[v1].lb();
            double const i1_ub = b1[v1].ub();
            b1_lits.push_back(m_store.add(v1, i1_lb, false));
            b1_lits.push_back(m_store.add(v1, i1_ub, true));
        }

        for (Enode * v_i : used_vars) {
            //     \/ !I1_j \/ (b2[v_i].lb <= v_i)
            //     j
            // --> \/ !I1_j \/ (v_i >= b2[v_i].lb)
            //     j
            for (int b1_lit : b1_lits) {
                picosat_add(m_psat, -b1_lit);
            }
            double const i2_lb = b2[v_i].lb();
            picosat_add(m_psat, m_store.add(v_i, i2_lb, false));  //i2_lb <= v_i --> v_i >= i2_lb;
            picosat_add(m_psat, 0);

            // \/ !I1_j \/ (v_i <= b2[v_i].ub)
            //  j
            for (int b1_lit : b1_lits) {
                picosat_add(m_psat, -b1_lit);
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
            assert(i1_lb <= i2_lb);
            if (i1_lb < i2_lb) {
                add_ordering(b1, used_vars, v, i1_lb, i2_lb);
                add_axiom(b2, used_vars, v, i2_lb);
            }

            assert(i2_lb <= i2_ub);
            if (i2_lb < i2_ub) {
                add_ordering(b2, used_vars, v, i2_lb, i2_ub);
                if (i2_ub < i1_ub) {
                    add_axiom(b2, used_vars, v, i2_ub);
                }
            }
            assert(i2_ub <= i1_ub);
            if (i2_ub < i1_ub) {
                add_ordering(b1, used_vars, v, i2_ub, i1_ub);
            }
        }
    }

    // Add B => l1 \/ l2 \/ l3 \/ l4
    void picosat_wrapper::add_imply(box const & b, vector<Enode *> const & used_var_vec, int const l1, int const l2, int const l3, int const l4) {
        unordered_set<Enode *> const used_vars(used_var_vec.begin(), used_var_vec.end());
        return add_imply(b, used_vars, l1, l2, l3, l4);
    }

    // Add B => l1 \/ l2 \/ l3 \/ l4
    void picosat_wrapper::add_imply(box const & b, unordered_set<Enode *> const & used_vars, int const l1, int const l2, int const l3, int const l4) {
        // Add !B
        vector<int> c;
        for (Enode * v : used_vars) {
            double const lb = b[v].lb();
            double const ub = b[v].ub();
            //     !((lb <= v) /\  (v <= ub))
            // --> !(lb <= v)  \/ !(v <= ub)
            // --> !(v >= lb)  \/ !(v <= ub)
            picosat_add(m_psat, -m_store.add(v, lb, false));
            picosat_add(m_psat, -m_store.add(v, ub, true));
            c.push_back(-m_store.add(v, lb, false));
            c.push_back(-m_store.add(v, ub, true));
        }
        picosat_add(m_psat, l1);
        c.push_back(l1);
        if (l2) {
            picosat_add(m_psat, l2);
            c.push_back(l2);
            if (l3) {
                picosat_add(m_psat, l3);
                c.push_back(l3);
                if (l4) {
                    picosat_add(m_psat, l4);
                    c.push_back(l4);
                }
            }
        }
        picosat_add(m_psat, 0);
        c.push_back(0);
    }

    // Add B => (v <= m) or (m <= v)
    void picosat_wrapper::add_branching(box const & b, Enode * v, double const m) {
        // TODO(soonhok): only do this if v m is not in the store
        double const lb = b[v].lb();
        double const ub = b[v].ub();
        DREAL_LOG_WARNING << "ADD_BRANCHING on "
                        << v << "[" << lb << ", " << m << ", " << ub << "]\n"
                        << b;
        assert(lb <= m);  add_ordering(b, v, lb, m);

        assert(m  <= ub); add_ordering(b, v, m, ub);
        if (lb < m && m < ub) {
            add_axiom(b, v, m);
        }

        assert(lb <= ub); add_ordering(b, v, lb, ub);
    }

    int picosat_wrapper::check_sat() {
        int const ret = picosat_sat(m_psat, -1);
        return ret;
    }
    // Precondition: check_sat() == PICOSAT_SATISFIABLE
    // Reduce the given box b into a smaller box using SAT model
    box picosat_wrapper::reduce_using_model(box b) {
        // TODO(soonhok): this can be a bottleneck. Consider an optimization.
        m_pmodel.clear();
        for (int i = 1; i <= m_store.get_num_vars(); i++) {
            int const r = picosat_deref_partial(m_psat, i);
            // pred := v <= bound
            tuple<Enode*, double, bool> pred = m_store.lookup(i);
            Enode * v = get<0>(pred);
            double const bound = get<1>(pred);
            bool const le = get<2>(pred);
            DREAL_LOG_FATAL << "b" << i << "\t"
                            << (r == 1 ? "+" : (r == 0 ? "0 " : "! "))
                            << v
                            << (le ? " <= " : " >= ")
                            << std::setprecision(16) << bound;

            if (r != 1) { continue;  /* UNKNOWN */ }
            m_pmodel.push_back(r * i);
            if (r == 1 && le) {
                // (v <= bound)
                DREAL_LOG_FATAL << "b[" << v << "] : "
                                  << b[v] << " /\\ " << ibex::Interval(b[v].lb(), bound)
                                  << " [" << b[v].lb() << ", " << bound << "]";
                b[v] &= ibex::Interval(b[v].lb(), bound);
                DREAL_LOG_FATAL << " = " << b[v];
            } else if (r == 1 && !le) {
                // (v >= bound)
                DREAL_LOG_FATAL << "b[" << v << "] : "
                                  << b[v] << " /\\ " << ibex::Interval(bound, b[v].ub())
                                  << " [" << bound << ", " << b[v].ub() << "]";
                b[v] &= ibex::Interval(bound, b[v].ub());
                DREAL_LOG_FATAL << " = " << b[v];
            } else if (r == -1 && le) {
                // !(v <= bound) --> (v > bound)
                if (bound == b[v].ub()) {
                    // b =  [                  ]
                    //                         |
                    //                       bound
                    DREAL_LOG_FATAL << "b[" << v << "] = " << b[v] << " intersect with "
                                    << v << " > " << bound << " = empty";
                    b[v].set_empty();
                } else {
                    // b = [                    ]
                    //                  |
                    //                  +----------------
                    DREAL_LOG_FATAL << "b[" << v << "] : "
                                    << b[v] << " /\\ " << ibex::Interval(bound, b[v].ub())
                                    << " (" << bound << ", " << b[v].ub() << "]";
                    b[v] &= ibex::Interval(bound, b[v].ub());
                    DREAL_LOG_FATAL << " = " << b[v];
                }
            } else if (r == -1 && !le) {
                // !(v >= bound) --> (v < bound)
                if (bound == b[v].lb()) {
                    // b =     [                     ]
                    //         |
                    //       bound
                    DREAL_LOG_FATAL << "b[" << v << "] = " << b[v] << " intersect with "
                                    << v << " < " << bound << " = empty";
                    b[v].set_empty();
                } else {
                    // b =     [                     ]
                    //                 |
                    //         --------+
                    DREAL_LOG_FATAL << "b[" << v << "] : "
                                    << b[v] << " /\\ " << ibex::Interval(b[v].lb(), bound)
                                    << " [" << b[v].lb() << ", " << bound << ")";
                    b[v] &= ibex::Interval(b[v].lb(), bound);
                    DREAL_LOG_FATAL << " = " << b[v];
                }
            } else {
                DREAL_LOG_FATAL << "?? The return value of picosat should be in {-1, 0, 1}";
                abort();
            }
            if (b[v].is_empty()) {
                b.set_empty();
                DREAL_LOG_FATAL << "SOMETHING IS WRONG, WE GOT AN EMPTY INTERVAL HERE";
                break;
            }
        }
        return b;
    }

    void picosat_wrapper::block_current_assignment() {
        cerr << "BLOCK: ";
        for (int const l : m_pmodel) {
            picosat_add(m_psat, -l);
            cerr << " " << -l;
        }
        picosat_add(m_psat, 0);
        cerr << " 0\n";
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
