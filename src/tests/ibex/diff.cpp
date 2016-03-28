/*********************************************************************
Author: Soonho Kong <soonhok@cs.cmu.edu>

dReal -- Copyright (C) 2013 - 2015, the dReal Team

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
#include "ibex/ibex.h"
#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include "catch/catch.hpp"

using ibex::Variable;
using ibex::Function;
using ibex::NumConstraint;
using ibex::CtcFwdBwd;
using ibex::IntervalVector;
using ibex::Interval;
using ibex::ExprSymbol;
using ibex::ExprCtr;
using ibex::ExprConstant;
using ibex::ExprNode;
using std::cout;
using std::endl;

TEST_CASE("box diff1") {
    cout << "============================== BOX DIFF 1 ========================" << endl;
    ibex::Interval iv1(0, 10);
    ibex::Interval iv2(5, 15);
    ibex::IntervalVector v1(2, iv1);
    ibex::IntervalVector v2(2, iv2);
    ibex::IntervalVector * diff_result = nullptr;
    int const diff_num = v1.diff(v2, diff_result);
    cout << "====== v1 ======" << endl;
    cout << v1 << endl;
    cout << "====== v2 ======" << endl;
    cout << v2 << endl << endl;
    for (int i = 0; i < diff_num; ++i) {
        cout << "====== diff " << i << " ======" << endl
             << diff_result[i] << endl;
    }
    delete[] diff_result;
}

TEST_CASE("box diff2") {
    cout << "============================== BOX DIFF 2 ========================" << endl;
    ibex::Interval iv1(0, 30);
    ibex::Interval iv2(10, 20);
    ibex::IntervalVector v1(2, iv1);
    ibex::IntervalVector v2(2, iv2);
    ibex::IntervalVector * diff_result = nullptr;
    int const diff_num = v1.diff(v2, diff_result);
    cout << "====== v1 ======" << endl;
    cout << v1 << endl;
    cout << "====== v2 ======" << endl;
    cout << v2 << endl << endl;
    for (int i = 0; i < diff_num; ++i) {
        cout << "====== diff " << i << " ======" << endl
             << diff_result[i] << endl;
    }
    delete[] diff_result;
}

/** b\b = emptyset */
TEST_CASE("TestIntervalVector::diff01()") {
    double _b[][2] = {{0, 1}, {0, 1}};
    IntervalVector b(2, _b);
    IntervalVector* c;
    int n = b.diff(b, c);

    REQUIRE(n == 0);

    REQUIRE(c[0].size() == 2);
    REQUIRE(c[0].is_empty());

    delete[] c;
}

/** b\emptyset = b */
TEST_CASE("TestIntervalVector::diff02()") {
    double _b[][2] = {{0, 1}, {0, 1}};
    IntervalVector b(2, _b);
    IntervalVector* c;

    int n = b.diff(IntervalVector::empty(2), c);

    REQUIRE(n == 1);
    REQUIRE(c[0].size() == 2);
    REQUIRE(c[0] == b);

    delete[] c;
}

/**
 * [-7,7]x[-7,7]  \ [-5,5]x[-5,5] =
 *  {  [-7,-5[x[-7,7] ; ]5,7]x[-7,7]  ; [-5,5]x[-7,-5[ ; [-5,5]x]5,7] }
 */
TEST_CASE("TestIntervalVector::diff03()") {
    double _b1[][2] = {{-7, 7}, {-7, 7}};
    double _b2[][2] = {{-5, 5}, {-5, 5}};
    IntervalVector b1(2, _b1);
    IntervalVector b2(2, _b2);
    IntervalVector* c;

    int n = b1.diff(b2, c);

    REQUIRE(n == 4);

    double _b3[][2] = {{-7, -5}, {-7, 7}};
    double _b4[][2] = {{ 5, 7}, {-7, 7}};
    double _b5[][2] = {{-5, 5}, {-7, -5}};
    double _b6[][2] = {{-5, 5}, { 5, 7}};
    IntervalVector b3(2, _b3);
    IntervalVector b4(2, _b4);
    IntervalVector b5(2, _b5);
    IntervalVector b6(2, _b6);

    REQUIRE(c[0] == b3);
    REQUIRE(c[1] == b4);
    REQUIRE(c[2] == b5);
    REQUIRE(c[3] == b6);

    delete[] c;
}

/**
 * [-7,7]x[-7,7]  \ [-7,7]x[-7,5] =
 *  {  [-7,7]x]5,7] }
 */
TEST_CASE("TestIntervalVector::diff04()") {
    double _b1[][2] = {{-7, 7}, {-7, 7}};
    double _b2[][2] = {{-7, 7}, {-7, 5}};
    IntervalVector b1(2, _b1);
    IntervalVector b2(2, _b2);
    IntervalVector* c;

    int n = b1.diff(b2, c);

    REQUIRE(n == 1);

    double _b3[][2] = {{-7, 7}, {5, 7}};
    IntervalVector b3(2, _b3);

    REQUIRE(c[0] == b3);

    delete[] c;
}

TEST_CASE("TestIntervalVector::diff05()") {
    double _b[][2] = {{0, 1}, {0, 1}};
    IntervalVector b(2, _b);
    IntervalVector* c;

    int n = IntervalVector::empty(2).diff(b, c);

    REQUIRE(n == 0);

    delete[] c;
}
