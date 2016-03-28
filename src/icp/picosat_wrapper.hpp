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

#pragma once

#include <functional>
#include <initializer_list>
#include <iostream>
#include <limits>
#include <map>
#include <random>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include "contractor/contractor.h"
#include "icp/pred_abs.h"
#include "opensmt/smtsolvers/SMTConfig.h"
#include "picosat/picosat.h"
#include "util/box.h"
#include "util/logging.h"
#include "util/stat.h"
#include "icp/gbox.h"
#include "icp/picosat_wrapper.h"

using std::cerr;
using std::endl;

namespace dreal {
template<class IT> void picosat_wrapper::add_clause(IT first, IT last) {
    assert(first != last);
    while(first != last) {
        int const l = *(first++);
        assert(l);
        picosat_add(m_psat, l);
    }
    picosat_add(m_psat, 0);
}

template<class IT> void picosat_wrapper::add_imply(gbox const & gb, IT l_first, IT const l_last) {
    std::vector<int> c;
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
    while (l_first != l_last) {
        int const l = *l_first;
        if (l) {
            c.push_back(l);
        } else {
            break;
        }
        ++l_first;
    }
    add_clause(c.cbegin(), c.cend());
}
}  // namespace dreal
