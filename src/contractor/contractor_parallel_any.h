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

#pragma once
#include <algorithm>
#include <atomic>
#include <cassert>
#include <condition_variable>
#include <initializer_list>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>
#include "./config.h"
#include "constraint/constraint.h"
#include "contractor/contractor_common.h"
#include "opensmt/egraph/Enode.h"
#include "opensmt/smtsolvers/SMTConfig.h"
#include "util/box.h"

namespace dreal {

// contractor_parallel_any
// - Run C1, C2, ... , Cn in parallel
// - If Ci returns SAT/Exception, then cancel the rest threads and propagate the result from Ci
// - Otherwise, return UNSAT.
class contractor_parallel_any : public contractor_cell {
private:
    std::vector<contractor> m_vec;
    std::mutex m_mutex;
    std::condition_variable m_cv;
    int m_index;
    void setup_input();

public:
    explicit contractor_parallel_any(std::initializer_list<contractor> const & l);
    explicit contractor_parallel_any(std::vector<contractor> const & v);
    contractor_parallel_any(contractor const & c1, contractor const & c2);
    void prune(contractor_status & s);
    std::ostream & display(std::ostream & out) const;
};

contractor mk_contractor_parallel_any(std::initializer_list<contractor> const & l);
contractor mk_contractor_parallel_any(std::vector<contractor> const & v);
contractor mk_contractor_parallel_any(contractor const & c1, contractor const & c2);

}  // namespace dreal
