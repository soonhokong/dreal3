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

#include <functional>
#include <iostream>
#include <limits>
#include <map>
#include <random>
#include <tuple>
#include <unordered_set>
#include <vector>
#include "icp/sat_icp.h"
#include "util/logging.h"
#include "util/stat.h"
#include "icp/clause_manager.h"

using std::cerr;
using std::cout;
using std::endl;
using std::get;
using std::map;
using std::numeric_limits;
using std::tuple;
using std::unordered_set;
using std::vector;

namespace dreal {
void sat_icp::solve(contractor & ctc, contractor_status & cs) {
    DREAL_LOG_FATAL << "INITIAL BOX" << endl
                    << cs.m_box;
    thread_local static std::unordered_set<std::shared_ptr<constraint>> used_constraints;
    used_constraints.clear();

    // TODO(soonhok): need to add clause_manager inside of contractor_status
    clause_manager cm(cs.m_box);
    box old_b(cs.m_box);
    while (cm.check_next_box()) {
        cs.m_box = cm.get_next_box();
        assert(!cs.m_box.is_empty());
        while (true) {
            // Pruning -- Begin
            try {
                old_b = cs.m_box;
                ctc.prune(cs);
                if (cs.m_config.nra_use_stat) { cs.m_config.nra_stat.increase_prune(); }
            } catch (contractor_exception & e) { /* do nothing */ }
            // Pruning -- End
            if (cs.m_box.is_empty()) {
                break;  // exit inner while-loop
            }
            // Branch
            tuple<int, box, box> const splits = cs.m_box.bisect(cs.m_config.nra_precision);
            int const i = get<0>(splits);  // branched at i-th dim
            if (i >= 0) {
                if (cs.m_config.nra_use_stat) { cs.m_config.nra_stat.increase_branch(); }
                box const & b1 = get<1>(splits);
                box const & b2 = get<2>(splits);
                cm.add_branch(cs.m_box, b1, b2);
                cs.m_box = b1;
            } else {
                // i < 0, which indicates it's not possible to bisect.
                DREAL_LOG_FATAL << "SAT-ICP: return delta-sat (|b| = " << cs.m_box.max_diam() << ")";
                return;
            }
        }  // while(true)
    }  // while(cm.check_next_box())
    cs.m_box.set_empty();
    DREAL_LOG_FATAL << "SAT-ICP: return unsat";
}
}  // namespace dreal
