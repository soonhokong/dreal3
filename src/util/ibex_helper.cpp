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

#include <limits>
#include <iostream>
#include <unordered_set>
#include "ibex/ibex.h"
#include "util/ibex_helper.h"
#include "util/ibex_interval_hash.h"
#include "util/logging.h"

using std::cerr;
using std::endl;
using std::numeric_limits;
using std::unordered_set;

namespace dreal {

bool interval_overlaps(ibex::Interval const & b1, ibex::Interval const & b2) {
    // Note: ibex::Interval::overlaps([0,1], [1,2]) returns false
    // while this method returns true.
    double const b1_lb = b1.lb();
    double const b1_ub = b1.ub();
    double const b2_lb = b2.lb();
    double const b2_ub = b2.ub();
    // [    b1    ]
    //        [    b2    ]
    // or
    //        [    b2    ]
    // [    b1    ]
    return (b1_lb <= b2_ub && b2_lb <= b1_ub) ||
        (b2_lb <= b1_ub && b1_lb <= b2_ub);
}

int interval_complementary(ibex::Interval const & i, ibex::Interval& c1, ibex::Interval& c2) {
    if (i.lb() == i.ub() && i.lb() != -numeric_limits<double>::infinity() && i.ub() != numeric_limits<double>::infinity()) {
        c1 = ibex::Interval(-numeric_limits<double>::infinity(), i.lb());
        c2 = ibex::Interval(i.ub(), +numeric_limits<double>::infinity());
        return 2;
    }
    if (i.is_empty() || i.is_degenerated()) {  // x.is_empty() should not happen if called from compl()
        c1 = ibex::Interval::ALL_REALS;
        c2 = ibex::Interval::EMPTY_SET;
        return 1;
    } else {
        if (i.lb() > -numeric_limits<double>::infinity()) {
            c1 = ibex::Interval(-numeric_limits<double>::infinity(), i.lb());
            if (i.ub() < +numeric_limits<double>::infinity()) {
                c2 = ibex::Interval(i.ub(), +numeric_limits<double>::infinity());
                return 2;
            } else {
                c2 = ibex::Interval::EMPTY_SET;
                return 1;
            }
        } else if (i.ub() < +numeric_limits<double>::infinity()) {
            c1 = ibex::Interval(i.ub(), +numeric_limits<double>::infinity());
            c2 = ibex::Interval::EMPTY_SET;
            return 1;
        } else {
            c1 = c2 = ibex::Interval::EMPTY_SET;
            return 0;
        }
    }
}

int interval_diff(ibex::Interval const & i1, ibex::Interval const & i2, ibex::Interval & c1, ibex::Interval & c2) {
    interval_complementary(i2, c1, c2);
    c1 &= i1;
    int res = 2;
    if (c1.is_degenerated()) { c1=ibex::Interval::EMPTY_SET; res--; }
    c2 &= i1;
    if (c2.is_degenerated()) { c2=ibex::Interval::EMPTY_SET; res--; }
    if (c1.is_empty()) {
        c1 = c2;
        c2 = ibex::Interval::EMPTY_SET;
        assert(res < 2);
    }
    if (res == 2) {
        assert(!c2.is_empty());
    }
    if (res >= 1) {
        assert(!c1.is_empty());
    }
    return res;
}

// return v1 \ v2
unordered_set<ibex::IntervalVector> interval_vector_diff(ibex::IntervalVector v1, ibex::IntervalVector const & v2) {
    ibex::IntervalVector org(v1);
    unordered_set<ibex::IntervalVector> ret;
    assert(v1.size() == v2.size());
    ibex::Interval tmp1;
    ibex::Interval tmp2;
    for (int i = 0; i < v1.size(); ++i) {
        auto const & i1 = v1[i];
        auto const & i2 = v2[i];
        if (i1 != i2) {
            int const diff_num = interval_diff(i1, i2, tmp1, tmp2);
            if (diff_num > 0) {
                v1[i] &= v2[i];
                ibex::IntervalVector new_v(v1);
                new_v[i] = tmp1;
                ret.insert(new_v);
                if (diff_num > 1) {
                    new_v[i] = tmp2;
                    ret.insert(new_v);
                }
            }
        }
    }
    return ret;
}
}  // namespace dreal
