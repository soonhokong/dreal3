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
    x->setDomainLowerBound(0);
    x->setDomainUpperBound(10);
    y->setDomainLowerBound(0);
    y->setDomainUpperBound(10);

    box b0({x, y});
    DREAL_LOG_FATAL << "b0:\n" << b0 << "\n=====================\n\n";
    picosat_wrapper pw;


    pw.add_box(b0);
    int ret = pw.check_sat();
    auto const b1 = pw.reduce_using_model(b0);
    // b1 : [0, 10] x [0, 10]
    DREAL_LOG_FATAL << "b1:\n" << b1 << "\n=====================\n\n";


    // PRUNE [0, 10] x [0, 10] => [2,10] x [0, 10]
    box b2(b0);
    b2["x"] = ibex::Interval(2, 9);
    b2["y"] = ibex::Interval(0, 10);
    pw.add_generalized_blocking_box(b1, b2, {x, y});
    ret = pw.check_sat();
    auto b3 = pw.reduce_using_model(b0);
    DREAL_LOG_FATAL << "b3:\n" << b3 << "\n=====================\n\n";


    // PRUNE [0, 10] x [0, 10] => [2,10] x [0, 10]
    box b4(b0);
    b4["x"] = ibex::Interval(2, 9);
    b4["y"] = ibex::Interval(5, 9);
    pw.add_generalized_blocking_box(b3, b4, {x, y});
    ret = pw.check_sat();
    auto b5 = pw.reduce_using_model(b0);
    DREAL_LOG_FATAL << "b5:\n" << b5 << "\n=====================\n\n";


    // BRANCH on x, at x = 6
    pw.add_branching(b5, x, 6);
    pw.add_branching(b5, y, 6);
    // ret = pw.check_sat();
    // auto b6 = pw.reduce_using_model(b0);
    // // b4 : [0, 5] x [0, 10]
    // DREAL_LOG_FATAL << "b6:\n" << b6 << "\n=====================\n\n";

    while (pw.check_sat() == PICOSAT_SATISFIABLE) {
        box b = pw.reduce_using_model(b0);
        DREAL_LOG_FATAL << "b:\n" << b << "\n=====================\n\n";
        pw.add_generalized_blocking_box(b, {x, y});
    }
    // // PRUNE [0, 5] x [0, 10] ==> [2, 3] x [3, 6]
    // box b21(b2);
    // b21["x"] = ibex::Interval(0, 3);
    // b21["y"] = ibex::Interval(6, 10);
    // pw.add_generalized_blocking_box(b2, b21, {x, y});
    // ret = pw.check_sat();
    // auto b3 = pw.reduce_using_model(b0);
    // DREAL_LOG_FATAL << "b3:\n" << b3 << "\n=====================\n\n";


    // // PRUNE [2,3] x [3,6] ==> EMPTY
    // pw.add_generalized_blocking_box(b3, {x, y});
    // ret = pw.check_sat();
    // auto b4 = pw.reduce_using_model(b0);
    // // b4 : [5, 10] x [0, 10]
    // DREAL_LOG_FATAL << "b4:\n" << b4 << "\n=====================\n\n";


    // pw.add_generalized_blocking_box(b4, {x, y}); // prune out b4
    // ret = pw.check_sat();
    // DREAL_LOG_FATAL << "ret = " << (ret == PICOSAT_UNSATISFIABLE);


    return 0;
}
