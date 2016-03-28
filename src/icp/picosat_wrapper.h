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

#include <functional>
#include <initializer_list>
#include <iostream>
#include <limits>
#include <map>
#include <random>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include "contractor/contractor.h"
#include "icp/pred_abs.h"
#include "icp/point_map.h"
#include "opensmt/smtsolvers/SMTConfig.h"
#include "picosat/picosat.h"
#include "util/box.h"
#include "icp/gbox.h"
#include "util/stat.h"

namespace dreal {

class picosat_wrapper {
private:
    gbox const & m_init_box;
    PicoSAT * m_psat;
    pred_abs m_store;

    void add_le(Enode * const v, double const c);  // Add: v <= c
    void add_le(double const c, Enode * const v);  // Add: c <= v
    void add_ge(Enode * const v, double const c);  // Add: v >= c
    void add_ge(double const c, Enode * const v);  // Add: c >= v
    void add_lt(Enode * const v, double const c);  // Add: v <  c
    void add_lt(double const c, Enode * const v);  // Add: c <  v
    void add_gt(Enode * const v, double const c);  // Add: v >  c
    void add_gt(double const c, Enode * const v);  // Add: c >  v

    void add_clause(std::initializer_list<int const> const & c);  // Add: c_1 \/ ... \/ c_n
    template<class IT> void add_clause(IT first, IT last);  // Add: first \/ ... \/ last - 1

    void add_imply(int const l1, int const l2);  // Add l1 => l2
    // Add b => l_1 \/ ... \/ l_n
    template<class IT> void add_imply(gbox const & b, IT l_first, IT const l_last);
    // Add b => l_1 \/ ... \/ l_n
    void add_imply(gbox const & b, std::initializer_list<int const> const & c);

    void add_box(gbox const & b);

public:
    explicit picosat_wrapper(gbox const & b);
    ~picosat_wrapper();

    void add_conflict(gbox const & b);  // Add: ¬b
    void add_imply(gbox const & gb1, gbox const & gb2);  // Add: ant => con
    // Given a variable `v` and two constants `l` and `u` where l < u
    // holds, add the followings:
    //  - (v <= lb) =>  (v <= ub)
    //  - (v >= ub) =>  (v >= lb)
    //  - (v <= lb) =>  !(v >= ub)
    void add_ordering(Enode * const v, double const l, double const u);
    void add_axiom(Enode *v, double const c);  // Add (v >= c) \/ (v <= c)
    bool check_sat();

    // Interpret a SAT model into a box
    gbox interpret_model(point_map const & pmap);

    // Debugging Functions
    void debug_print(int const l) const;
    void debug_print(std::vector<int> const & c) const;
    void debug_print(bool const print_picosat = false) const;
};

}  // namespace dreal

#include "icp/picosat_wrapper.hpp"
