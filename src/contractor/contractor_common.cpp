/*********************************************************************
Author: Soonho Kong <soonhok@cs.cmu.edu>

dReal -- Copyright (C) 2013 - 2015, the dReal Team

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
#include <functional>
#include <initializer_list>
#include <iterator>
#include <limits>
#include <map>
#include <memory>
#include <queue>
#include <random>
#include <set>
#include <sstream>
#include <stack>
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
#include "constraint/constraint.h"
#include "contractor/contractor_common.h"
#include "contractor/contractor_basic.h"
#include "ibex/ibex.h"
#include "opensmt/egraph/Enode.h"
#include "util/box.h"
#include "util/logging.h"
#include "util/proof.h"
#include "util/interruptible_thread.h"

using std::back_inserter;
using std::cerr;
using std::cout;
using std::dynamic_pointer_cast;
using std::endl;
using std::function;
using std::initializer_list;
using std::make_pair;
using std::make_shared;
using std::ostream;
using std::ostringstream;
using std::pair;
using std::queue;
using std::set;
using std::shared_ptr;
using std::unordered_map;
using std::unordered_set;
using std::vector;

namespace dreal {
std::ostream & operator<<(std::ostream & out, ode_direction const & d) {
    switch (d) {
    case ode_direction::FWD:
        out << "FWD";
        break;
    case ode_direction::BWD:
        out << "BWD";
        break;
    }
    return out;
}

ostream & operator<<(ostream & out, contractor_cell const & c) {
    return c.display(out);
}

box generalize_box(box b, box const & init_box, unordered_set<Enode *> const & used_vars) {
    for (unsigned i = 0 ; i < b.size(); ++i) {
        Enode * const v = b.get_var(i);
        if (used_vars.find(v) == used_vars.end()) {
            // If i-th var is not used, relax i-th dim in the box
            assert(b[i].is_subset(init_box[i]));
            b.relax(i, init_box[i]);
        }
    }
    return b;
}

void contractor::prune(box & b, SMTConfig & config, clause_manager * const cm_ptr) {
    if (m_ptr) {
        // by default, clear output vector and used constraints.
        m_ptr->clear_output();
        m_ptr->clear_used_constraints();
        if (cm_ptr) {
            box old_b(b);
            m_ptr->prune(b, config, cm_ptr);
            if (b.is_empty()) {
                bool const generalize = true;
                if (generalize) {
                    cm_ptr->add_conflict(old_b, m_ptr->used_constraints());
                } else {
                    cm_ptr->add_conflict(old_b);
                }
            } else {
                if (old_b == b) { return; }
                bool const generalize = true;
                if (generalize) {
                    cm_ptr->add_imply(old_b, b, m_ptr->used_constraints());
                } else {
                    cm_ptr->add_imply(old_b, b);
                }
            }
        } else {
            m_ptr->prune(b, config, cm_ptr);
        }
    }
}

void contractor::prune_with_assert(box & b, SMTConfig & config, clause_manager * const cm_ptr) {
    assert(m_ptr != nullptr);
    thread_local static box old_box(b);
    old_box = b;
    m_ptr->prune(b, config, cm_ptr);
    if (!b.is_subset(old_box)) {
        cerr << "Pruning Violation: " << (*m_ptr) << endl;
        cerr << "Old Box" << endl
             << "==============" << endl
             << old_box << endl;
        cerr << "New Box" << endl
             << "==============" << endl
             << b << endl;
        cerr << "==============" << endl;
        display_diff(cerr, old_box, b);
        cerr << "==============" << endl;
        exit(1);
    }
}
}  // namespace dreal
