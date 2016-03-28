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
#include <iostream>
#include <unordered_set>
#include <limits>
#include "icp/reduced_box_set.h"
#include "util/logging.h"
#include "util/ibex_helper.h"

using std::endl;
using std::min;
using std::max;
using std::ostream;
using std::unordered_set;

namespace dreal {

void reduced_box_set::add_and_simplify(gbox const & b) {
    // If there exists b_ in set which is included by b, remove b_
    auto it = m_set.begin();
    while (it != m_set.end()) {
        if (it->is_subset(b)) {
            // erase returns an iterator pointing to the position
            // immediately following the last of the elements erased
            it = m_set.erase(it);
        } else {
            ++it;
        }
    }
    // 1. for all `b_` in the set,
    for (gbox const & b_ : m_set) {
        // If b_ includes b, don't add b
        if (b_.is_superset(b)) {
            return;
        }
        // if it is possible to merge `b` and `b_`, pop(b_) from the
        // set, call `add_and_simplify(merge(b, b_))` again.
        int const merge_dim = find_merge_dim(b, b_);
        if (merge_dim >= 0) {
            gbox const b__ = b_;
            m_set.erase(b_);
            add_and_simplify(merge(b, b__, merge_dim));
            add_and_simplify(merge(b__, b, merge_dim));
            return;
        }
    }
    // b is not included in the set, and it's not mergable.
    // let's add it to the set
    m_set.emplace(b);
}

void reduced_box_set::add(gbox const & gb) {
    if (exist(gb)) { return; }
    if (!check_invariant()) {
        DREAL_LOG_FATAL << __func__ << " Before add, invariant failed";
        assert(false);
    }
    add_and_simplify(gb);
    if (!check_invariant()) {
        DREAL_LOG_FATAL << __func__ << " After add, invariant failed";
        assert(false);
    }
}

bool trivially_included(gbox const & gb, unordered_set<gbox> const & s) {
    for (gbox const & gb_ : s) {
        if (gb.is_subset(gb_)) { return true; }
    }
    return false;
}

bool trivially_not_included(gbox const & gb, unordered_set<gbox> const & s) {
    for (gbox const & gb_ : s) {
        if (gb.is_strict_overlap(gb_)) { return false; }
    }
    return true;
}

bool includes_core(gbox const & gb, unordered_set<gbox> const & s) {
    if (s.size() == 1) {
        return gb.is_subset(*s.begin());
    }
    if (trivially_included(gb, s)) {
        return true;
    }

    for (gbox const & gb_i : s) {
        if (gb.is_strict_overlap(gb_i)) {
            unordered_set<gbox> const diff_boxes = gb.set_minus(gb_i);
            // s.erase(gb_i);
            for (gbox const & gb_j : diff_boxes) {
                if (!includes_core(gb_j, s)) {
                    return false;
                }
            }
            return true;
        }
    }
    return false;
}

bool reduced_box_set::includes(gbox const & gb) const {
    if (trivially_included(gb, m_set)) {
        return true;
    }
    if (trivially_not_included(gb, m_set)) {
        return false;
    }
    // unordered_set<gbox> s;
    // for (gbox const & gb_i : m_set) {
    //     if (gb_i.is_strict_overlap(gb)) {
    //         s.insert(gb_i);
    //     }
    // }
    if (m_set.size() == 1) {
        return gb.is_subset(*m_set.begin());
    }
    bool const ret = includes_core(gb, m_set);
    return ret;
}

ostream & operator<<(ostream & out, reduced_box_set const & rs) {
    for (gbox const & b : rs.m_set) {
        out << b << endl;
        out << "---------------------" << endl;
    }
    return out;
}

bool reduced_box_set::check_invariant() const {
    bool ret = true;
    for (auto it1 = cbegin(); it1 != cend(); ++it1) {
        for (auto it2 = cbegin(); it2 != cend(); ++it2) {
            if (it1 != it2) {
                if (it1->is_subset(*it2)) {
                    DREAL_LOG_FATAL << "There is subset ordering between the following two boxes, but both of" << endl
                                    << "them are found in the reduced box set, which shouldn't be the case.";
                    DREAL_LOG_FATAL << *it1;
                    DREAL_LOG_FATAL << "===============";
                    DREAL_LOG_FATAL << *it2;
                    DREAL_LOG_FATAL << "===============";
                    ret = false;
                    DREAL_LOG_FATAL << "REDUCED BOX SET: INVARIANT TYPE1";
                }
                int const mergable_dim = find_merge_dim(*it1, *it2);
                if (mergable_dim >= 0) {
                    DREAL_LOG_FATAL << "There are mergable boxes in the set of conflict boxes.";
                    DREAL_LOG_FATAL << *it1;
                    DREAL_LOG_FATAL << "===========";
                    DREAL_LOG_FATAL << *it2;
                    DREAL_LOG_FATAL << "===========";
                    ret = false;
                    DREAL_LOG_FATAL << "REDUCED BOX SET: INVARIANT TYPE2";
                }
            }
        }
    }
    return ret;
}
}  // namespace dreal
