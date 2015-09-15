/*********************************************************************
Author: Soonho Kong <soonhok@cs.cmu.edu>

dReal -- Copyright (C) 2013 - 2015, Soonho Kong, Sicun Gao, and Edmund Clarke

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
#include <ios>
#include <chrono>
#include "tools/dop/stat.h"

using std::chrono::duration;
using std::chrono::duration_cast;
using std::chrono::steady_clock;
using std::cout;
using std::endl;
using std::fixed;
using std::milli;

namespace dop {

stat::stat(std::ostream & out) : m_out(out) {
    m_start = std::chrono::steady_clock::now();
}

void stat::print() const {
    std::chrono::duration<double> diff = std::chrono::steady_clock::now() - m_start;
    m_out << "Running Time: "
          << fixed
          << diff.count()
          << " sec" << endl;

}
}  // namespace dop
