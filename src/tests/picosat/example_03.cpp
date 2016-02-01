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
#include "picosat/picosat.h"
#include "icp/picosat_wrapper.h"
#include "api/OpenSMTContext.h"

using namespace dreal;

int main() {
    OpenSMTContext ctx;
    ctx.SetLogic(QF_NRA);
    Snode * real_sort = ctx.mkSortReal();
    ctx.DeclareFun("x", real_sort);
    ctx.DeclareFun("y", real_sort);
    Enode * x = ctx.mkVar("x");
    Enode * y = ctx.mkVar("y");
    x->setDomainLowerBound(-1);
    x->setDomainUpperBound(5);
    y->setDomainLowerBound(-1);
    y->setDomainUpperBound(4);

    box b0({x, y});
    b0["x"] = ibex::Interval(0, 4);
    b0["y"] = ibex::Interval(0, 3);
    DREAL_LOG_FATAL << "b0:\n" << b0 << "\n=====================\n\n";
    picosat_wrapper pw;


    pw.add_box(b0);
    int ret = pw.check_sat();
    auto const b1 = pw.reduce_using_model(b0);
    // b1 : [0, 4] x [0, 3]
    DREAL_LOG_FATAL << "b1:\n" << b1 << "\n=====================\n\n";


    // PRUNE [0, 10] x [0, 10] => [2,10] x [0, 10]
    box b2(b0);
    b2["x"] = ibex::Interval(1, 4);
    b2["y"] = ibex::Interval(0, 2);
    DREAL_LOG_FATAL << "b2:\n" << b2 << "\n=====================\n\n";
    DREAL_LOG_FATAL << "b1 => b2\n";
    pw.add_generalized_blocking_box(b1, b2, {x, y});
    ret = pw.check_sat();
    auto b3 = pw.reduce_using_model(b0);
    DREAL_LOG_FATAL << "b3:\n" << b3 << "\n=====================\n\n";

    pw.add_branching(b3, x, 3);
    ret = pw.check_sat();
    auto b4 = pw.reduce_using_model(b0);
    DREAL_LOG_FATAL << "b4:\n" << b4 << "\n=====================\n\n";

    box b5(b4);
    b5["x"] = ibex::Interval(1, 2);
    b5["y"] = ibex::Interval(1, 2);
    DREAL_LOG_FATAL << "b5:\n" << b5 << "\n=====================\n\n";

    DREAL_LOG_FATAL << "b4 => b5\n";
    pw.add_generalized_blocking_box(b4, b5, {x, y});

    while (pw.check_sat() == PICOSAT_SATISFIABLE) {
        box b = pw.reduce_using_model(b0);
        pw.add_generalized_blocking_box(b, {x, y});
        DREAL_LOG_FATAL << "block the following:\n" << b << "\n=====================\n\n";
    }

    return 0;
}
