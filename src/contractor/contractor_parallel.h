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
#include "contractor/contractor_common.h"
#include "opensmt/egraph/Enode.h"
#include "opensmt/smtsolvers/SMTConfig.h"
#include "util/box.h"
#include "constraint/constraint.h"

namespace dreal {

enum class pruning_thread_status {
    READY, RUNNING, SAT, UNSAT, EXCEPTION, KILLED
        };

std::ostream & operator<<(std::ostream & out, pruning_thread_status const & s);

void parallel_helper_fn(unsigned const id, contractor & c, contractor_status & cs, pruning_thread_status & ts,
                        std::mutex & m, std::condition_variable & cv, int & index, std::atomic_int & tasks_to_run);
}  // namespace dreal
