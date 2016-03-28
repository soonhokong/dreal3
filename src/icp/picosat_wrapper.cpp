/*********************************************************************
Author: Soonho Kong <soonhok@cs.cmu.edu>

dReal -- Copyright (C) 2013 - 2016, the dReal Team

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

#include <initializer_list>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <vector>
#include "icp/picosat_wrapper.h"
#include "icp/point_map.h"

using std::cerr;
using std::endl;
using std::get;
using std::initializer_list;
using std::set;
using std::tuple;
using std::unordered_set;
using std::vector;

namespace dreal {
    picosat_wrapper::picosat_wrapper(gbox const & b) : m_init_box(b) {
        m_psat = picosat_init();
        // picosat_save_original_clauses(m_psat);  // to use picosat_deref_partial function
        add_box(b);
    }
    picosat_wrapper::~picosat_wrapper() {
        picosat_reset(m_psat);
    }

    // Add: v <= c
    void picosat_wrapper::add_le(Enode * const v, double const c) {
        add_clause({m_store.add_le(v, c)});
        DREAL_LOG_WARNING << "PICOSAT WRAPPER: ADD - (" << v << " <= " << c << ")";
    }
    // Add: c <= v  <-->  (v >= c)
    void picosat_wrapper::add_le(double const c, Enode * const v) {
        add_ge(v, c);
    }
    // Add: v >= c
    void picosat_wrapper::add_ge(Enode * const v, double const c) {
        add_clause({m_store.add_ge(v, c)});
        DREAL_LOG_WARNING << "PICOSAT WRAPPER: ADD - (" << v << " >= " << c << ")";
    }
    // Add: c >= v  <-->  (v <= c)
    void picosat_wrapper::add_ge(double const c, Enode * const v) {
        add_le(v, c);
    }
    // Add: v < c
    void picosat_wrapper::add_lt(Enode * const v, double const c) {
        add_clause({m_store.add_lt(v, c)});
        DREAL_LOG_WARNING << "PICOSAT WRAPPER: ADD - (" << v << " < " << c << ")";
    }
    // Add: c < v  <-->  (v > c)
    void picosat_wrapper::add_lt(double const c, Enode * const v) {
        add_gt(v, c);
    }
    // Add: v > c
    void picosat_wrapper::add_gt(Enode * const v, double const c) {
        add_clause({m_store.add_gt(v, c)});
        DREAL_LOG_WARNING << "PICOSAT WRAPPER: ADD - (" << v << " > " << c << ")";
    }
    // Add: c > v  <-->  (v < c)
    void picosat_wrapper::add_gt(double const c, Enode * const v) {
        add_lt(v, c);
    }
    void picosat_wrapper::add_clause(initializer_list<int const> const & c) {
        add_clause(begin(c), end(c));
    }
    void picosat_wrapper::add_imply(int const l1, int const l2) {
        add_clause({-l1, l2});  // Add: l1 => l2, that is, (!l1 \/ l2)
    }

    // Given a variable `v` and two constants `l` and `u` where l <= u holds.
    void picosat_wrapper::add_ordering(Enode * const v, double const l, double const u) {
        assert(l <= u);
        if (l < u) {
            DREAL_LOG_WARNING << "ADD_ORDERING: " << v << " " << l << " " << u;
            //  (v <= lb) =>  (v <= ub)
            add_imply(m_store.add_le(v, l), m_store.add_le(v, u));
            //  (v >= ub) =>  (v >= lb)
            add_imply(m_store.add_ge(v, u), m_store.add_ge(v, l));
            //  (v <= lb) => !(v >= ub)
            add_imply(m_store.add_le(v, l), -m_store.add_ge(v, u));
        }
    }

    // (v <= c) \/ (v >= c)
    void picosat_wrapper::add_axiom(Enode * const v, double const c) {
        add_clause({m_store.add_le(v, c), m_store.add_ge(v, c)});
    }

    void picosat_wrapper::add_box(gbox const & gb) {
        for (unsigned i = 0; i < gb.size(); ++i) {
            if (!gb.get_gbit(i)) {
                Enode * const v = gb.get_var(i);
                double const lb = gb.get_value(i).lb();
                double const ub = gb.get_value(i).ub();
                assert(lb <= ub);
                // lb <= v <= ub
                add_le(lb, v);
                add_le(v, ub);
            }
        }
    }

    bool picosat_wrapper::check_sat() {
        int const ret = picosat_sat(m_psat, -1);
        if (ret == PICOSAT_SATISFIABLE) {
            return true;
        } else if (ret == PICOSAT_UNSATISFIABLE) {
            return false;
        }
        assert(ret == PICOSAT_UNKNOWN);
        DREAL_LOG_FATAL << "picosat_wrapper::check_sat: PICOSAT returns PICOSAT_UNKNOWN";
        abort();
    }

    // Interpret the current Boolean sat result into a box
    gbox picosat_wrapper::interpret_model(point_map const & pmap) {
        assert(picosat_res(m_psat) == PICOSAT_SATISFIABLE);
        gbox gb(m_init_box);
        for (unsigned i = 0; i < gb.size(); ++i) {
            Enode * const v = gb.get_var(i);
            if (!pmap.has_key(v)) {
                continue;
            }
            set<double> const & s = pmap.at(v);
            bool need_to_update = false;
            // find upperbound of v, search smaller value first
            double ub = gb.get_value(i).ub();  // current value
            for (auto it = s.cbegin(); it != s.cend(); ++it) {
                double const c = *it;
                int const i = m_store.lookup(v, c, true);  // i <-> (v <= c)
                if (picosat_deref(m_psat, i) == 1) {
                    ub = c;
                    need_to_update = true;
                    break;
                }
            }
            // find lowerbound of v, search bigger value first
            double lb = gb.get_value(i).lb();  // current value
            for (auto it = s.crbegin(); it != s.crend(); ++it) {
                double const c = *it;
                int const i = m_store.lookup(v, c, false);  // i <-> (v >= c)
                if (picosat_deref(m_psat, i) == 1) {
                    lb = c;
                    need_to_update = true;
                    break;
                }
            }
            if (need_to_update) {
                gb.set_value(i, lb, ub);
            }
        }
        return gb;
    }

    void picosat_wrapper::debug_print(int const l) const {
        if (l == 0) {
            cerr << "0"; return;
        }
        bool const is_neg = l < 0;
        auto const p = m_store.lookup(is_neg ? -l : l);
        Enode * const v = get<0>(p);
        double const bound = get<1>(p);
        bool const le = get<2>(p);
        if (is_neg) {
            cerr << "!";
        }
        if (le) {
            cerr << "(" << v << " <= " << bound << ")";
        } else {
            cerr << "(" << v << " >= " << bound << ")";
        }
    }

    void picosat_wrapper::debug_print(vector<int> const & c) const {
        for (int const l : c) {
            debug_print(l);
            cerr << " ";
        }
        cerr << endl;
    }

    void picosat_wrapper::debug_print(bool const print_picosat) const {
        DREAL_LOG_FATAL << "======================";
        for (int i = 1; i <= m_store.get_num_vars(); ++i) {
            tuple<Enode *, double, bool> pred = m_store.lookup(i);
            Enode * const v = get<0>(pred);
            double const bound = get<1>(pred);
            bool const le = get<2>(pred);
            DREAL_LOG_FATAL << "b" << i << " := "
                            << "(" << v
                           << (le ? " <= " : " >= ")
                            << std::setprecision(16) << bound << ")";
        }
        if (print_picosat) {
            DREAL_LOG_FATAL << "~~~~~~~~~~~~~~~~~~~~~~";
            picosat_print(m_psat, stderr);
        }
        DREAL_LOG_FATAL << "======================";
    }

    void picosat_wrapper::add_conflict(gbox const & gb) {
        thread_local static vector<int> c;
        c.clear();
        for (unsigned i = 0; i < gb.size(); ++i) {
            if (!gb.get_gbit(i)) {
                Enode * const v = gb.get_var(i);
                double const lb = gb[i].lb();
                double const ub = gb[i].ub();
                // !((lb <= v) /\  (v <= ub)) --> !(lb <= v)  \/ !(v <= ub)
                c.push_back(-m_store.add_le(lb, v));
                c.push_back(-m_store.add_le(v, ub));
            }
        }
        assert(c.size() > 0);
        add_clause(c.begin(), c.end());
    }

    void picosat_wrapper::add_imply(gbox const & b, std::initializer_list<int const> const & c) {
        add_imply(b, begin(c), end(c));
    }

    void picosat_wrapper::add_imply(gbox const & gb1, gbox const & gb2) {
        assert(gb1 != gb2);
#ifndef NDEBUG
        if (!gb1.get_values().is_superset(gb2.get_values())) {
            DREAL_LOG_FATAL << "B1 => B2, but B1 is not a superset of B2:" << endl
                            << "B1 = " << gb1 << endl << endl
                            << "B2 = " << gb2 << endl << endl;
            assert(false);
        }
#endif
        assert(gb1.is_superset(gb2));
        // B1 -> B2   ==>  B1 -> /\ B2[i]
        for (unsigned i = 0; i < gb2.size(); ++i) {
            if (!gb2.get_gbit(i)) {
                Enode * const v = gb2.get_var(i);
                assert(gb1.get_var(i) == v);
                double const i1_lb = gb1.get_value(i).lb();
                double const i1_ub = gb1.get_value(i).ub();
                double const i2_lb = gb2.get_value(i).lb();
                double const i2_ub = gb2.get_value(i).ub();
                assert(i1_lb <= i2_lb); assert(i2_lb <= i2_ub); assert(i2_ub <= i1_ub);
                //     B1 -> (b2[v_i].lb <= v_i)
                if (i1_lb < i2_lb) {
                    add_imply(gb1, {m_store.add_le(i2_lb, v)});  // i2_lb <= v_i
                }
                if (i2_ub < i1_ub) {
                    //     B1 -> (v_i <= b2[v_i].ub)
                    add_imply(gb1, {m_store.add_le(v, i2_ub)});  // v_i <= i2_ub
                }
            }
        }
    }
}  // namespace dreal
