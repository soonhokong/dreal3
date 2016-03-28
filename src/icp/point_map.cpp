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

#include <vector>
#include "icp/point_map.h"

using std::vector;

namespace dreal {
void point_map::add(gbox const & gb) {
    for (unsigned i = 0; i < gb.size(); ++i) {
        if (!gb.get_gbit(i)) {
            Enode * const v = gb.get_var(i);
            double const lb = gb[i].lb();
            double const ub = gb[i].ub();
            add_intv(v, lb, ub);
        }
    }
}
void point_map::add(vector<gbox> const & v) {
    for (gbox const & gb : v) {
        add(gb);
    }
}
}  // namespace dreal
