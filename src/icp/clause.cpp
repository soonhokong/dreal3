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

#include <iostream>
#include "icp/clause.h"

using std::ostream;
using std::endl;

namespace dreal {
ostream & operator<<(ostream & out, clause const & c) {
    out << "ant:" << endl
        << "====" << endl
        << c.m_antecedent << endl;
    out << "con:" << endl;
    for (auto const & cl : c.m_consequent) {
        out << cl << endl
            << "---" << endl;
    }
    return out;
}
}  // namespace dreal
