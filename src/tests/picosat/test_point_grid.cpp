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
#include "icp/point_grid.h"
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
    b0["x"] = ibex::Interval(0, 10);
    b0["y"] = ibex::Interval(0, 10);
    DREAL_LOG_FATAL << "========== INITIAL BOX =============";
    DREAL_LOG_FATAL << b0;
    Grid g(b0);
    g.debug_print();

    DREAL_LOG_FATAL << "\n\n========== ADD_POINT(x, 3)==========";
    g.add_point(x, 3);
    g.debug_print();


    DREAL_LOG_FATAL << "\n\n========== ADD_POINT(x, 6)==========";
    g.add_point(x, 6);
    g.debug_print();


    // DREAL_LOG_FATAL << "\n\n========== ADD_POINT(y, 5)==========";
    // g.add_point(y, 5);
    // g.debug_print();

    DREAL_LOG_FATAL << "\n\n========== PUSH NO BOUNDS FORMULA ==========";
    g.debug_print_clause(g.get_push_nobounds_formula());

    DREAL_LOG_FATAL << "\n\n========== PUSH BOUNDS ONLY FORMULA ==========";
    g.debug_print_clause(g.get_push_bounds_only_formula());

    DREAL_LOG_FATAL << "\n\n========== PUSH FORMULA ==========";
    g.debug_print_clause(g.get_push_formula());

    return 0;
}
