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
#include <vector>
#include "opensmt/api/OpenSMTContext.h"
#include "opensmt/egraph/Enode.h"
#include "opensmt/egraph/Egraph.h"
#include "util/box.h"
#include "icp/gbox.h"
#include "icp/reduced_box_set.h"
#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include "catch/catch.hpp"

using std::cerr;
using std::endl;
using std::vector;

namespace dreal {
bool g_initialized = false;
OpenSMTContext * g_context_ptr;
Snode * g_real_sort;
vector<Enode *> g_vars;

void init() {
    if (!g_initialized) {
        g_context_ptr = new OpenSMTContext();
        g_real_sort = g_context_ptr->mkSortReal();
        g_context_ptr->DeclareFun("x", g_real_sort);
        g_context_ptr->DeclareFun("y", g_real_sort);
        auto x = g_context_ptr->mkVar("x");
        auto y = g_context_ptr->mkVar("y");
        g_vars.push_back(x);
        g_vars.push_back(y);
        g_initialized = true;
    }
}

gbox make_gbox(double const x_l, double const x_u, double const y_l, double const y_u) {
    init();
    gbox gb(g_vars);
    gb.set_value(0, x_l, x_u);
    gb.set_value(1, y_l, y_u);
    return gb;
}

TEST_CASE("!=") {
    gbox const b1 = make_gbox(0, 10, 0, 1);  // [0, 10] x [0, 1]
    gbox const b2 = make_gbox(5 , 7, 1, 2);  // [5,  7] x [1, 2]
    REQUIRE(b1 != b2);
}

TEST_CASE("==") {
    gbox const b1 = make_gbox(0, 10, 0, 1);  // [0, 10] x [0, 1]
    gbox const b2 = make_gbox(0, 10, 0, 1);  // [0, 10] x [0, 1]
    REQUIRE(b1 == b2);
}

}  // namespace dreal
