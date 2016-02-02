/*********************************************************************
Author: Soonho Kong <soonhok@cs.cmu.edu>

dReal -- Copyright (C) 2013 - 2016, Soonho Kong, Sicun Gao, and Edmund Clarke

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

#include "contractor/contractor.h"
#include "icp/point_grid.h"
#include "opensmt/smtsolvers/SMTConfig.h"
#include "picosat/picosat.h"
#include "util/box.h"
#include "util/scoped_vec.h"
#include "util/stat.h"

namespace dreal {

class gsat_icp {
private:
    box const m_initial_box;
    Grid m_grid;
    PicoSAT * m_ps;

    // add (l <= v) /\ (v <= u)
    void add_interval(Enode * v, double const l, double const u);
    // add b => l1 \/ l2 \/ l3 \/ l4
    void add_imply(box const & b, int const l1 = 0, int const l2 = 0, int const l3 = 0, int const l4 = 0);
    void add_vector(std::vector<int> const & vec);
    box build_box_from_sat_model();
    // given a box `b`, Add `!b`
    void add_learned_clause(box const & b);
    // given boxes `b1` and `b2`, add `b1 => b2`, that is `!b1 \/ b2`
    void add_learned_clause(box const & b1, box const & b2);
    void add_branch(box const & b, Enode * v, double const p);

public:
    gsat_icp(box const & b);
    ~gsat_icp();

    // TODO(soonhok): later we need to have
    // `scoped_vec<std::shared_ptr<constraint>> const & ctrs` to have
    // fine-grained learning.
    box solve(contractor & ctc, SMTConfig & config);
};

}  // namespace dreal
