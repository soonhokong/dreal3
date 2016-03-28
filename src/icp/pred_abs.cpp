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

#include "icp/pred_abs.h"
#include "util/logging.h"

using std::make_tuple;
using std::get;

namespace dreal {
    int pred_abs::add(Enode * const v, double const bound, bool const le) {
        auto const p = make_tuple(v, bound, le);
        auto const it = m_abs_map.find(p);
        if (it == m_abs_map.end()) {
            ++m_num_vars;
            m_con_map.emplace(m_num_vars, p);
            m_abs_map.emplace(p, m_num_vars);
            return m_num_vars;
        } else {
            return it->second;
        }
    }
    void pred_abs::debug_print() const {
        for (int i = 1; i <= m_num_vars; ++i) {
            auto const t = lookup(i);
            DREAL_LOG_FATAL << "B" << i << "\t <---> \t"
                            << get<0>(t)
                            << (get<2>(t) ? " <= " : " >= ")
                            << get<1>(t);
        }
    }


}  // namespace dreal
