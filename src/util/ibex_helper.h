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

#include <unordered_set>
#include "util/box.h"
#include "ibex/ibex.h"

namespace dreal {
bool interval_overlaps(ibex::Interval const & b1, ibex::Interval const & b2);
int interval_complementary(ibex::Interval const & i, ibex::Interval& c1, ibex::Interval& c2);
int interval_diff(ibex::Interval const & i1, ibex::Interval const & i2, ibex::Interval & c1, ibex::Interval & c2);
std::unordered_set<ibex::IntervalVector> interval_vector_diff(ibex::IntervalVector v1, ibex::IntervalVector const & v2);
}  // namespace dreal
