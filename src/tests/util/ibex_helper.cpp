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
#include "util/ibex_helper.h"
#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include "catch/catch.hpp"

using std::cerr;
using std::endl;
using std::vector;
using std::numeric_limits;

namespace dreal {
TEST_CASE("interval overlaps 1") {
    ibex::Interval const i1(0, 1);
    ibex::Interval const i2(0, 1);
    REQUIRE(interval_overlaps(i1, i2));
}
TEST_CASE("interval overlaps 2") {
    ibex::Interval const i1(0, 1);
    ibex::Interval const i2(0.5, 1.5);
    REQUIRE(interval_overlaps(i1, i2));
}
TEST_CASE("interval overlaps 3") {
    ibex::Interval const i1(0, 1);
    ibex::Interval const i2(1, 2);
    REQUIRE(interval_overlaps(i1, i2));
}
TEST_CASE("interval overlaps 4") {
    ibex::Interval const i1(0, 1);
    ibex::Interval const i2(0);
    REQUIRE(interval_overlaps(i1, i2));
}
TEST_CASE("interval overlaps 5") {
    ibex::Interval const i1(0, 1);
    ibex::Interval const i2(1);
    REQUIRE(interval_overlaps(i1, i2));
}

TEST_CASE("interval complementary 1") {
    ibex::Interval const i(0, 1);
    ibex::Interval r1;
    ibex::Interval r2;
    int const r = interval_complementary(i, r1, r2);
    REQUIRE(r == 2);
    REQUIRE(r1.lb() == -numeric_limits<double>::infinity());
    REQUIRE(r1.ub() == 0);
    REQUIRE(r2.lb() == 1);
    REQUIRE(r2.ub() == numeric_limits<double>::infinity());
}

TEST_CASE("interval complementary 2") {
    ibex::Interval const i(0);
    ibex::Interval r1;
    ibex::Interval r2;
    int const r = interval_complementary(i, r1, r2);
    REQUIRE(r == 2);
    REQUIRE(r1.lb() == -numeric_limits<double>::infinity());
    REQUIRE(r1.ub() == 0);
    REQUIRE(r2.lb() == 0);
    REQUIRE(r2.ub() == numeric_limits<double>::infinity());
}

TEST_CASE("interval diff 1") {
    // [0, 2] \ [1, 3] = [0, 1]
    ibex::Interval const i1(0, 2);
    ibex::Interval const i2(1, 3);
    ibex::Interval r1;
    ibex::Interval r2;
    int const r = interval_diff(i1, i2, r1, r2);
    REQUIRE(r == 1);
    REQUIRE(r1.lb() == 0);
    REQUIRE(r1.ub() == 1);
}

TEST_CASE("interval diff 2") {
    // [-2, 2] \ [-1, 1] = {[-2, -1], [1, 2]}
    ibex::Interval const i1(-2, 2);
    ibex::Interval const i2(-1, 1);
    ibex::Interval r1;
    ibex::Interval r2;
    int const r = interval_diff(i1, i2, r1, r2);
    REQUIRE(r == 2);
    REQUIRE(r1.lb() == -2);
    REQUIRE(r1.ub() == -1);
    REQUIRE(r2.lb() ==  1);
    REQUIRE(r2.ub() ==  2);
}

TEST_CASE("interval diff 3") {
    // [-2, 2] \ -1 = {[-2, -1], [-1, 2]}
    ibex::Interval const i1(-2, 2);
    ibex::Interval const i2(-1);
    ibex::Interval r1;
    ibex::Interval r2;
    int const r = interval_diff(i1, i2, r1, r2);
    REQUIRE(r == 2);
    REQUIRE(r1.lb() == -2);
    REQUIRE(r1.ub() == -1);
    REQUIRE(r2.lb() == -1);
    REQUIRE(r2.ub() ==  2);
}

TEST_CASE("interval diff 4") {
    // [-2, 10] \ -2 = [-2, 10]
    ibex::Interval const i1(-2, 10);
    ibex::Interval const i2(-2);
    ibex::Interval r1;
    ibex::Interval r2;
    int const r = interval_diff(i1, i2, r1, r2);
    REQUIRE(r == 1);
    REQUIRE(r1.lb() == -2);
    REQUIRE(r1.ub() == 10);
}

TEST_CASE("interval diff 5") {
    // [-2, 10] \ 10 = [-2, 10]
    ibex::Interval const i1(-2, 10);
    ibex::Interval const i2(10);
    ibex::Interval r1;
    ibex::Interval r2;
    int const r = interval_diff(i1, i2, r1, r2);
    REQUIRE(r == 1);
    REQUIRE(r1.lb() == -2);
    REQUIRE(r1.ub() == 10);
}

TEST_CASE("interval diff 6") {
    // [-2, 10] \  = [-2, 10]
    ibex::Interval const i1(-2, 10);
    ibex::Interval const i2(10);
    ibex::Interval r1;
    ibex::Interval r2;
    int const r = interval_diff(i1, i2, r1, r2);
    REQUIRE(r == 1);
    REQUIRE(r1.lb() == -2);
    REQUIRE(r1.ub() == 10);
}

TEST_CASE("interval diff 7") {
    // [0, 2] \ [-1, 1] = [1, 2]
    ibex::Interval const i1(0, 2);
    ibex::Interval const i2(-1, 1);
    ibex::Interval r1;
    ibex::Interval r2;
    int const r = interval_diff(i1, i2, r1, r2);
    REQUIRE(r == 1);
    REQUIRE(r1.lb() == 1);
    REQUIRE(r1.ub() == 2);
}

TEST_CASE("interval diff 8") {
    // [0, 2] \ [3, 4] = [0, 2]
    ibex::Interval const i1(0, 2);
    ibex::Interval const i2(3, 4);
    ibex::Interval r1;
    ibex::Interval r2;
    int const r = interval_diff(i1, i2, r1, r2);
    REQUIRE(r == 1);
    REQUIRE(r1.lb() == 0);
    REQUIRE(r1.ub() == 2);
}

TEST_CASE("interval diff 9") {
    // [0, 2] \ [-4, -3] = [0, 2]
    ibex::Interval const i1(0, 2);
    ibex::Interval const i2(-4, -3);
    ibex::Interval r1;
    ibex::Interval r2;
    int const r = interval_diff(i1, i2, r1, r2);
    REQUIRE(r == 1);
    REQUIRE(r1.lb() == 0);
    REQUIRE(r1.ub() == 2);
}

TEST_CASE("interval vector diff 1") {
    ibex::IntervalVector const v1({{0, 10}, {0, 10}});
    ibex::IntervalVector const v2({{5, 15}, {5, 15}});
    auto const diff_result = interval_vector_diff(v1, v2);

    // expected output
    REQUIRE(diff_result.size() == 2);
    ibex::IntervalVector r1({{0,  5}, {0, 10}});
    ibex::IntervalVector r2({{5, 10}, {0,  5}});
    REQUIRE(diff_result.count(r1) == 1);
    REQUIRE(diff_result.count(r2) == 1);
}

TEST_CASE("interval vector diff 2") {
    ibex::IntervalVector const v1({{ 0, 30}, { 0, 30}});
    ibex::IntervalVector const v2({{10, 20}, {10, 20}});
    auto const diff_result = interval_vector_diff(v1, v2);

    // expected output
    REQUIRE(diff_result.size() == 4);
    ibex::IntervalVector r1({{ 0, 10}, { 0, 30}});
    ibex::IntervalVector r2({{20, 30}, { 0, 30}});
    ibex::IntervalVector r3({{10, 20}, { 0, 10}});
    ibex::IntervalVector r4({{10, 20}, {20, 30}});
    REQUIRE(diff_result.count(r1) == 1);
    REQUIRE(diff_result.count(r2) == 1);
    REQUIRE(diff_result.count(r3) == 1);
    REQUIRE(diff_result.count(r4) == 1);
}

TEST_CASE("interval vector diff 3") {
    ibex::IntervalVector const v1({{ 0, 20}, { 0, 20}});
    ibex::IntervalVector const v2({{10, 10}, {10, 10}});
    auto const diff_result = interval_vector_diff(v1, v2);
}

TEST_CASE("interval vector diff 4") {
    ibex::IntervalVector const v1({{ 0, 20}, { 0, 20}, { 0, 20}});
    ibex::IntervalVector const v2({{10, 10}, {10, 10}, {10, 10}});
    auto const diff_result = interval_vector_diff(v1, v2);
}
}  // namespace dreal
