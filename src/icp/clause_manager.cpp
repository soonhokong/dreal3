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

#include <algorithm>
#include <chrono>
#include <ios>
#include <iostream>
#include <set>
#include <vector>
#include "icp/reduced_box_set.h"
#include "icp/clause_manager.h"
#include "icp/picosat_wrapper.h"
#include "util/box.h"
#include "icp/gbox.h"
#include "icp/point_map.h"

using std::boolalpha;
using std::chrono::duration;
using std::chrono::steady_clock;
using std::endl;
using std::milli;
using std::remove_if;
using std::set;
using std::vector;
using std::unordered_set;
using std::shared_ptr;

namespace dreal {

clause_manager::clause_manager(box const & b)
    : m_init_box(b),
      m_found_next_box(false),
      m_next_box(m_init_box) { }

// Add b1 => b2 to the clause stack
void clause_manager::add_imply(box const & b1, box const & b2) {
    if (!check_invariant()) {
        DREAL_LOG_FATAL << __func__ << " Before add_imply(box, box), invariant failed";
        abort();
    }
    assert(b1 != b2);
    // ------- MAIN BEGIN -------
    m_clauses.emplace_back(gbox{b1}, gbox{b2});
    simplify();
    // ------- MAIN END -------
    if (!check_invariant()) {
        DREAL_LOG_FATAL << __func__ << " After add_imply(box, box), invariant failed";
        abort();
    }
}

void clause_manager::add_imply(box const & b1, box const & b2, unordered_set<shared_ptr<constraint>> const & used_ctrs) {
#ifndef NDEBUG
    if (!check_invariant()) {
        DREAL_LOG_FATAL << __func__ << " Before add_imply(box, box, unordered_set), invariant failed";
        abort();
    }
#endif
    assert(b1 != b2);
    assert(b1.is_superset(b2));
    gbox const gb1(b1, used_ctrs, m_init_box.get_values());
    gbox const gb2(b2, used_ctrs, m_init_box.get_values());
    assert(gb1.is_superset(gb2));
    // ------- MAIN BEGIN -------
    if (gb1 != gb2) {
        m_clauses.emplace_back(gb1, gb2);
    }
    simplify();
    // ------- MAIN END -------
#ifndef NDEBUG
    if (!check_invariant()) {
        DREAL_LOG_FATAL << __func__ << " After add_imply(box, box, unordered_set), invariant failed";
        abort();
    }
#endif
}

// Add !b into the conflict stack
void clause_manager::add_conflict(box const & b) {
    if (!check_invariant()) {
        DREAL_LOG_FATAL << __func__ << " Before add_conflict(box const &), invariant failed";
        abort();
    }
    // ------- MAIN BEGIN -------
    m_conflict_boxes.add(gbox{b});
    simplify();
    // ------- MAIN END -------
    if (!check_invariant()) {
        DREAL_LOG_FATAL << __func__ << " After add_conflict(box const &), invariant failed";
        abort();
    }
}

void clause_manager::add_conflict(box const & b, unordered_set<shared_ptr<constraint>> const & used_ctrs) {
    if (!check_invariant()) {
        DREAL_LOG_FATAL << __func__ << " Before add_conflict(box const &, unordered_set), invariant failed";
        abort();
    }
    // ------- MAIN BEGIN -------
    m_conflict_boxes.add(gbox{b, used_ctrs, m_init_box.get_values()});
    simplify();
    // ------- MAIN END -------
    if (!check_invariant()) {
        DREAL_LOG_FATAL << __func__ << " After add_conflict(box const &, unordered_set), invariant failed";
        abort();
    }
}

// Add b => (b /\ (v >= c)) \/ (b /\ (v <= c)) to the cluase stack
void clause_manager::add_branch(box const & b, box const & b1, box const & b2) {
    if (!check_invariant()) {
        DREAL_LOG_FATAL << __func__ << " Before branch, invariant failed";
        abort();
    }
    // ------- MAIN BEGIN -------
    m_clauses.emplace_back(gbox{b}, gbox{b1}, gbox{b2});
    simplify();
    // ------- MAIN END -------
    if (!check_invariant()) {
        DREAL_LOG_FATAL << __func__ << " After branch, invariant failed";
        abort();
    }
}

// Check if there is a box to traverse for an ICP. If so, save it to
// m_next_box and set m_found_next_box flag on.
bool clause_manager::check_next_box() {
    if (!check_invariant()) {
        DREAL_LOG_FATAL << __func__ << " Before check_next_box, invariant failed";
        abort();
    }
    simplify();
    // 1. Pass a query to SAT solver
    gbox init_gbox(m_init_box);
    picosat_wrapper pw(init_gbox);
    point_map pmap;

    // add clauses
    for (auto const & c : m_clauses) {
        gbox const & ant = c.m_antecedent;
        gbox const & con = c.m_consequent[0];
        pw.add_imply(ant, con);
        pmap.add(ant);
        pmap.add(con);
    }
    // add conflict boxes
    for (auto const & b : m_conflict_boxes) {
        pw.add_conflict(b);
        pmap.add(b);
    }
    // add points
    for (auto const & p : pmap) {
        Enode * const v = p.first;
        set<double> const & s = p.second;
        // add axiom
        for (double const c : s) {
            pw.add_axiom(v, c);
        }
        // add ordering
        if (s.size() < 2) {
            continue;
        }
        auto it = s.cbegin();
        double c1 = *(it++);
        while (it != s.cend()) {
            double const c2 = *(it++);
            pw.add_ordering(v, c1, c2);
            c1 = c2;
        }
    }
    // 2. Check SAT
    bool const check_result = pw.check_sat();
    DREAL_LOG_FATAL << "|boxes| = " << m_conflict_boxes.size() << "\t"
                    << "|clauses| = " << m_clauses.size();
    if (check_result) {
        // SAT case, interpret the result
        m_next_box.set_values(pw.interpret_model(pmap).get_values());
        m_found_next_box = true;
    } else {
        // UNSAT case
        m_found_next_box = false;
        assert(sanity_check());
    }
    if (!check_invariant()) {
        DREAL_LOG_FATAL << __func__ << " After check_next_box, invariant failed";
        abort();
    }
    return m_found_next_box;
}

bool clause_manager::check_invariant() const {
    bool ret = true;
    for (auto const & cl : m_clauses) {
        gbox const & ant = cl.m_antecedent;
        for (gbox const & gb : cl.m_consequent) {
            if (m_conflict_boxes.includes(gb)) {
                DREAL_LOG_FATAL << "We detect an unresolved clause.";
                DREAL_LOG_FATAL << cl;
                DREAL_LOG_FATAL << "where the following box in the consequent list" << endl
                                << "is included in the current set of conflict boxes";
                DREAL_LOG_FATAL << gb;
                ret = false;
                DREAL_LOG_FATAL << "CLAUSE MANAGER: INVARIANT TYPE1";
            }
        }
        if (cl.m_consequent.size() == 1) {
            gbox const & con = cl.m_consequent[0];
            if (ant == con) {
                DREAL_LOG_FATAL << "We detect a cluase which is of form ant => con where ant = con holds.";
                DREAL_LOG_FATAL << ant;
                DREAL_LOG_FATAL << "================";
                DREAL_LOG_FATAL << con;
                DREAL_LOG_FATAL << "================";
                ret = false;
                DREAL_LOG_FATAL << "CLAUSE MANAGER: INVARIANT TYPE2";
            }

        }
    }
    return ret;
}

bool clause_manager::sanity_check() {
    gbox init_gbox(m_init_box);
    bool ret = true;
    ret &= check_invariant();
    for (auto const & cl : m_clauses) {
        if (cl.m_consequent.size() != 1) {
            DREAL_LOG_FATAL << "When this sanity function is called, it should be the case that all" << endl
                            << "clauses have a single box on their consequents." << endl
                            << "But we find a case where this invariant does not hold (size = " << cl.m_consequent.size() << ")";
            DREAL_LOG_FATAL << "================";
            DREAL_LOG_FATAL << cl;
            DREAL_LOG_FATAL << "================";
            auto const & gb1 = cl.m_consequent[0];
            auto const & gb2 = cl.m_consequent[1];
            DREAL_LOG_FATAL << "B1 is included in conflict boxes?: " << m_conflict_boxes.includes(gb1);
            DREAL_LOG_FATAL << "B2 is included in conflict boxes?: " << m_conflict_boxes.includes(gb2);
            DREAL_LOG_FATAL << "SANITY CHECK TYPE1";
            ret = false;
        }
    }
    // By the time we have UNSAT, it should be the case that it's
    // already shown that there is no solution in the initial search
    // space.
    if (!m_conflict_boxes.includes(init_gbox)) {
        ret = false;
        DREAL_LOG_FATAL << "===============";
        DREAL_LOG_FATAL << "Conflict Boxes";
        DREAL_LOG_FATAL << "===============";
        DREAL_LOG_FATAL << m_conflict_boxes;
        DREAL_LOG_FATAL << "===============";
        DREAL_LOG_FATAL << "Clauses";
        DREAL_LOG_FATAL << "===============";
        for (auto const & c : m_clauses) {
            DREAL_LOG_FATAL << c;
        }
        cerr << "|conflict_boxes| = " << m_conflict_boxes.size() << endl;
        cerr << "|clauses|        = " << m_clauses.size() << endl;
        DREAL_LOG_FATAL << "SANITY CHECK TYPE2";
        ret = false;
    }
    ret &= m_conflict_boxes.check_invariant();
    return ret;
}

// Return the box found by check_next_box
box clause_manager::get_next_box() const {
    assert(m_found_next_box);
    DREAL_LOG_WARNING << "clause_manager::get_next_box() = " << endl
                      << m_next_box << endl;
    return m_next_box;
}

// Resolve `c` using `b` and return true if `c` is changed.
bool clause_manager::resolve(clause & c) {
    DREAL_LOG_WARNING << "clause_manager::resolve()" << endl;
    // Say that the consequent of c is ` => b_1 \/ ... \/ b_n`
    // we remove b_i if we know that b_i doesn't have a solution in it.
    // Return true if `c` is updated, return false otherwise.
    auto & con = c.m_consequent;
    auto const len = con.size();
    con.erase(remove_if(con.begin(), con.end(),
                        [this](gbox const & gb) {
                            return this->m_conflict_boxes.includes(gb);
                        }),
              con.end());
    return con.size() != len;
}

void clause_manager::simplify() {
    for (int ii = 0; ii < 10; ++ii) {
    bool fixedpoint_reached = true;
    do {
        fixedpoint_reached = true;
        unsigned i = 0;
        while (i < m_clauses.size()) {
            clause & c = m_clauses[i];
            // First resolve c. That is, try to remove a box in the
            // consequent of c by checking whether the box is already
            // included in the set of conflict boxes or not.
            bool const resolved = resolve(c);
            if (resolved) { fixedpoint_reached = false; }
            switch (c.m_consequent.size()) {
            case 0:
                // c is a form of `b => false`, that is we have !b.
                m_conflict_boxes.add(c.m_antecedent);
                // pop c from m_clauses by moving the last element
                // to the current position (i-th)
                c = m_clauses.back();
                m_clauses.pop_back();
                continue;  // need to revisit the current i again
            case 1:
                if (c.m_antecedent.generalized()) {
                    assert(c.m_consequent[0].generalized());
                    // c is of form b => b', and b and b' are generalized.
                    std::vector<gbox> gen_queue;
                    for (gbox const & gb : m_conflict_boxes) {
                        gbox const & ant = c.m_antecedent;
                        gbox const & con = c.m_consequent[0];
                        if (gb.match(ant, con)) {
                            gbox const gen = generalize(gb, ant, con);
                        }
                    }
                    for (auto const & gen : gen_queue) {
                        m_conflict_boxes.add(gen);
                    }
                    fixedpoint_reached |= gen_queue.size() > 0;
                }
                break;
            case 2:
                // The only case that we have more than one box in
                // consequent is in branching. In that case, every box
                // in a cluase is not generalized.
#ifndef NDEBUG
                assert(!c.m_antecedent.generalized());
                for (gbox const & gb : c.m_consequent) {
                    assert(!gb.generalized());
                }
#endif
                break;
            default:
                DREAL_LOG_FATAL << "shouldn't be reachable";
                break;
            }
            ++i;
        }
    } while (!fixedpoint_reached);
    }
}
}  // namespace dreal
