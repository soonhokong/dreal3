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
#include <limits>
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
using std::numeric_limits;

namespace dreal {
bool g_initialized = false;
OpenSMTContext * g_context_ptr;
Snode * g_real_sort;
vector<Enode *> g_vars;

static double constexpr inf = numeric_limits<double>::infinity();

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

TEST_CASE("== and !=") {
    gbox const b1 = make_gbox(0, 10, 0, 1);  // [0, 10] x [0, 1]
    gbox const b2 = make_gbox(5, 7, 1, 2);   // [5,  7] x [1, 2]
    gbox const b3 = make_gbox(0, 10, 0, 1);  // [0, 10] x [0, 1]
    REQUIRE(b1 == b1);
    REQUIRE(b2 == b2);
    REQUIRE(b1 != b2);
    REQUIRE(b3 == b3);
}

TEST_CASE("get_gbit, set_gbit") {
    gbox b1 = make_gbox(0, 10, 0, 1);  // [0, 10] x [0, 1]
    gbox b2 = make_gbox(0, 10, 1, 2);  // [0, 10] x [1, 2]
    REQUIRE(b1 != b2);

    REQUIRE(!b1.get_gbit(0));
    REQUIRE(!b1.get_gbit(1));
    REQUIRE(!b2.get_gbit(0));
    REQUIRE(!b2.get_gbit(1));

    b1.set_gbit(1, true);
    b2.set_gbit(1, true);

    REQUIRE(b1.get_gbit(1));
    REQUIRE(b2.get_gbit(1));

    REQUIRE(b1 == b2);
    REQUIRE(b1.hash() == b2.hash());
}

TEST_CASE("set minus - no gen") {
    cerr << "set minus -1\n";
    gbox const b1 = make_gbox(0, 5, 0, 5);  // [0, 5] x [0, 5]
    gbox const b2 = make_gbox(1, 4, 1, 4);  // [1, 4] x [1, 4]
    auto const & ret = b1.set_minus(b2);

    REQUIRE(ret.size() == 4);
    REQUIRE(ret.count(make_gbox(0, 1, 0, 5)) == 1);
    REQUIRE(ret.count(make_gbox(4, 5, 0, 5)) == 1);
    REQUIRE(ret.count(make_gbox(1, 4, 0, 1)) == 1);
    REQUIRE(ret.count(make_gbox(1, 4, 4, 5)) == 1);
}

TEST_CASE("set minus - gen b1(0)") {
    cerr << "set minus - gen b1(0)";
    gbox b1 = make_gbox(0, 5, 0, 5);        //      [0, 5] x [0, 5]
    gbox const b2 = make_gbox(1, 4, 1, 4);  //      [1, 4] x [1, 4]
    b1.set_gbit(1, true);                   // b1 = [0, 5] x [ENTIRE]
    auto const & ret = b1.set_minus(b2);
    for (auto const & b : ret) {
        cerr << b << endl;
        cerr << b[1].lb() << " " << b1[1].ub() <<endl;
        cerr << "~~~~~~~~~~~~\n";
    }
    REQUIRE(ret.size() == 4);
    auto ret1 = make_gbox(0, 1, -inf, inf);
    auto ret2 = make_gbox(4, 5, -inf, inf);
    auto ret3 = make_gbox(1, 4, -inf, 1);
    auto ret4 = make_gbox(1, 4, 4, inf);
    ret1.set_gbit(1, true);
    ret2.set_gbit(1, true);
    REQUIRE(ret.count(ret1) == 1);
    REQUIRE(ret.count(ret2) == 1);
    REQUIRE(ret.count(ret3) == 1);
    REQUIRE(ret.count(ret4) == 1);
}

// TEST_CASE("set minus - gen b1(1)") {
//     cerr << "set minus - 2\n";
//     gbox b1 = make_gbox(0, 5, 0, 5);  // [0, 5] x [0, 5]
//     gbox b2 = make_gbox(1, 4, 1, 4);  // [1, 4] x [1, 4]
//     b1.set_gbit(1, true);
//     auto const & ret = b1.set_minus(b2);

//     for (auto const & b : ret) {
//         cerr << b << endl;
//         cerr << "~~~~~~~~~~~~\n";
//     }
// }

// TEST_CASE("set minus - gen b1(0) and b1(1)") {
//     cerr << "set minus - 2\n";
//     gbox b1 = make_gbox(0, 5, 0, 5);  // [0, 5] x [0, 5]
//     gbox b2 = make_gbox(1, 4, 1, 4);  // [1, 4] x [1, 4]
//     b1.set_gbit(1, true);
//     auto const & ret = b1.set_minus(b2);

//     for (auto const & b : ret) {
//         cerr << b << endl;
//         cerr << "~~~~~~~~~~~~\n";
//     }
// }

// // TEST_CASE("set minus - 3") {
// //     cerr << "set minus - 3\n";
// //     gbox b1 = make_gbox(0, 5, 0, 5);  // [0, 5] x [0, 5]
// //     gbox b2 = make_gbox(1, 4, 1, 4);  // [1, 4] x [1, 4]
// //     b2.set_gbit(1, true);
// //     auto const & ret = b1.set_minus(b2);

// //     for (auto const & b : ret) {
// //         cerr << "~~~~~~~~~~~~" << endl;
// //         cerr << b << endl;
// //     }
// // }

// // TEST_CASE("set minus - 4") {
// //     cerr << "set minus - 4\n";
// //     gbox b1 = make_gbox(0, 5, 0, 5);  // [0, 5] x [0, 5]
// //     gbox b2 = make_gbox(1, 4, 1, 4);  // [1, 4] x [1, 4]
// //     b1.set_gbit(0, true);
// //     b1.set_gbit(1, true);
// //     auto const & ret = b1.set_minus(b2);

// //     for (auto const & b : ret) {
// //         cerr << "~~~~~~~~~~~~" << endl;
// //         cerr << b << endl;
// //     }
// // }

// // TEST_CASE("set minus - 5") {
// //     cerr << "set minus - 4\n";
// //     gbox b1 = make_gbox(0, 5, 0, 5);  // [0, 5] x [0, 5]
// //     gbox b2 = make_gbox(1, 4, 1, 4);  // [1, 4] x [1, 4]
// //     b1.set_gbit(0, true);
// //     b1.set_gbit(1, true);
// //     auto const & ret = b1.set_minus(b2);

// //     for (auto const & b : ret) {
// //         cerr << "~~~~~~~~~~~~" << endl;
// //         cerr << b << endl;
// //     }
// // }

}  // namespace dreal
