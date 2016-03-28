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
#include <memory>
#include <random>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "constraint/constraint.h"
#include "icp/clause.h"
#include "icp/reduced_box_set.h"
#include "util/box.h"
#include "util/stat.h"

namespace dreal {

class clause_manager {
private:
    box                 m_init_box;
    bool                m_found_next_box;
    box                 m_next_box;
    std::vector<clause> m_clauses;
    reduced_box_set     m_conflict_boxes;

private:
    bool resolve(clause & c);
    void simplify();
    bool sanity_check();
    bool check_invariant() const;  // for debugging purpose

public:
    explicit clause_manager(box const & b);
    void add_imply(box const & b1, box const & b2);
    void add_conflict(box const & b);
    void add_imply(box const & b1, box const & b2, std::unordered_set<std::shared_ptr<constraint>> const & used_ctrs);
    void add_conflict(box const & b, std::unordered_set<std::shared_ptr<constraint>> const & used_ctrs);
    void add_branch(box const & b, box const & b1, box const & b2);
    box get_next_box() const;
    bool check_next_box();
};
}  // namespace dreal
