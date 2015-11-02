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

#include <atomic>
#include <algorithm>
#include <chrono>
#include <exception>
#include <functional>
#include <initializer_list>
#include <iterator>
#include <future>
#include <limits>
#include <map>
#include <memory>
#include <queue>
#include <random>
#include <set>
#include <sstream>
#include <stack>
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
#include "constraint/constraint.h"
#include "contractor/contractor_common.h"
#include "contractor/contractor_basic.h"
#include "contractor/contractor_parallel.h"
#include "contractor/contractor_parallel_all.h"
#include "ibex/ibex.h"
#include "opensmt/egraph/Enode.h"
#include "util/box.h"
#include "util/logging.h"
#include "util/proof.h"
#include "util/interruptible_thread.h"

using std::async;
using std::back_inserter;
using std::cerr;
using std::condition_variable;
using std::cout;
using std::dynamic_pointer_cast;
using std::endl;
using std::function;
using std::future;
using std::initializer_list;
using std::make_pair;
using std::make_shared;
using std::mutex;
using std::min;
using std::max;
using std::ostream;
using std::ostringstream;
using std::pair;
using std::queue;
using std::ref;
using std::set;
using std::shared_ptr;
using std::unique_lock;
using std::lock_guard;
using std::unordered_map;
using std::unordered_set;
using std::vector;
using std::atomic_int;
using std::thread;
using std::exception;

namespace dreal {

contractor_parallel_all::contractor_parallel_all(initializer_list<contractor> const & l)
    : contractor_cell(contractor_kind::PARALLEL_ALL), m_vec(l) { }
contractor_parallel_all::contractor_parallel_all(vector<contractor> const & v)
    : contractor_cell(contractor_kind::PARALLEL_ALL), m_vec(v) { }
contractor_parallel_all::contractor_parallel_all(contractor const & c1, contractor const & c2)
    : contractor_cell(contractor_kind::PARALLEL_ALL), m_vec(1, c1) { m_vec.push_back(c2); }

void contractor_parallel_all::prune(box & b, SMTConfig & config, vector<box> &) {
    DREAL_LOG_DEBUG << "contractor_parallel_all::prune";
    DREAL_LOG_FATAL << "-------------------------------------------------------------";
    // TODO(soonhok): implement this
    if (m_vec.size() == 0) {
        // Do nothing for empty vec
        return;
    }

    // 1. Make n copies of box b
    vector<box> boxes(m_vec.size(), b);
    vector<vector<box>> bins(m_vec.size());
    vector<pruning_thread_status> statuses(m_vec.size(), pruning_thread_status::READY);
    m_index = -1;

    // DREAL_LOG_FATAL << "parallel: Boxes are copied";

    // 2. Trigger execution with each contractor and a copied box
    vector<interruptible_thread> threads;
    atomic_int tasks_to_run(m_vec.size());
    // DREAL_LOG_FATAL << "parallel: tasks to run = " << tasks_to_run.load();
    for (unsigned i = 0; i < m_vec.size(); ++i) {
        DREAL_LOG_FATAL << "parallel : thread " << i << " / " << (tasks_to_run.load() - 1)
                        << " spawning...";
        threads.emplace_back(parallel_helper_fn,
                             i,
                             m_vec[i],
                             boxes[i],
                             bins[i],
                             config,
                             statuses[i],
                             m_mutex,
                             m_cv,
                             m_index,
                             tasks_to_run);

        // TODO(soonhok): need to handle this bins[i] (side effect of this extension of prune interface)

        DREAL_LOG_FATAL << "parallel : thread " << i << " / " << (tasks_to_run.load() - 1)
                        << " spawned...";
    }
    DREAL_LOG_FATAL << "parallel : " << m_vec.size() << " thread(s) got created";

    while (true) {
        DREAL_LOG_FATAL << "parallel: waiting for the lock";
        unique_lock<mutex> lk(m_mutex);
        DREAL_LOG_FATAL << "parallel: get a lock. " << tasks_to_run.load() << " tasks to go";
        if (tasks_to_run.load() == 0) {
            break;
        }
        DREAL_LOG_FATAL << "parallel: WAIT for CV." << tasks_to_run.load() << " tasks to go";;
        m_index = -1;
        m_cv.wait(lk, [&]() { return m_index != -1; });
        DREAL_LOG_FATAL << "parallel: wake up" << tasks_to_run.load();
        pruning_thread_status const & s = statuses[m_index];
        // DREAL_LOG_FATAL << "parallel: thread " << m_index << " " << s;
        if (s == pruning_thread_status::UNSAT || s == pruning_thread_status::EXCEPTION) {
            // Interrupt all the rest threads
            for (unsigned i = 0; i < statuses.size(); i++) {
                if (i - m_index != 0 && (statuses[i] == pruning_thread_status::READY || statuses[i] == pruning_thread_status::RUNNING)) {
                    threads[i].interrupt();
                }
            }

            if (s == pruning_thread_status::UNSAT) {
                DREAL_LOG_FATAL << "parallel: " << m_index << " got UNSAT";
                b.set_empty();
                m_input.union_with(m_vec[m_index].input());
                m_output.union_with(m_vec[m_index].output());
                unordered_set<shared_ptr<constraint>> const & used_ctrs = m_vec[m_index].used_constraints();
                m_used_constraints.insert(used_ctrs.begin(), used_ctrs.end());
                lk.unlock();
                for (unsigned i = 0; i < m_vec.size(); i++) {
                    threads[i].join();
                }
                DREAL_LOG_FATAL << "parallel: return UNSAT";
                return;
            }
            if (s == pruning_thread_status::EXCEPTION) {
                DREAL_LOG_FATAL << "parallel: " << m_index << " got EXCEPTION";
                lk.unlock();
                for (unsigned i = 0; i < m_vec.size(); i++) {
                    threads[i].join();
                }
                DREAL_LOG_FATAL << "parallel: throw exception";
                throw contractor_exception("exception during parallel contraction");
            }

        } else {
            // if (s != pruning_thread_status::SAT) {
            //     // DREAL_LOG_FATAL << "parallel: " << m_index << " got " << s;
            //     // DREAL_LOG_FATAL << "parallel: " << m_index << " got " << statuses[m_index];
            assert(s == pruning_thread_status::SAT);
            // }
            // if (threads[m_index].joinable()) {
            // threads[m_index].join();
            // }
            // DREAL_LOG_FATAL << "parallel: " << m_index << " got SAT";
            // Why?
            //  - Not READY/RUNNING: It's a job already done.
            //  - Not UNSAT/EXCEPTION: already handled above.
            //  - Not KILLED: There must be one which kill the killed
            //                job, and this loop stops after handling
            //                the first one
        }
    }

    // Assertion: All of them got SAT
    // for (pruning_thread_status const & s : statuses) {
    //     assert(s == pruning_thread_status::SAT);
    // }

    // DREAL_LOG_FATAL << "All of them are SAT";
    b = boxes[0];
    for (unsigned i = 0; i < m_vec.size(); i++) {
        contractor const & c = m_vec[i];
        b.intersect(boxes[i]);
        m_input.union_with(c.input());
        m_output.union_with(c.output());
        unordered_set<shared_ptr<constraint>> const & used_ctrs = c.used_constraints();
        m_used_constraints.insert(used_ctrs.begin(), used_ctrs.end());
        if (b.is_empty()) {
            // DREAL_LOG_FATAL << "Found an empty while intersecting...";
            for (unsigned i = 0; i < m_vec.size(); i++) {
                // if (threads[i].joinable()) {
                // DREAL_LOG_FATAL << "Try to join " << i << "...";
                threads[i].join();
                // DREAL_LOG_FATAL << "Try to join " << i << "... done";
                // }
            }
            // DREAL_LOG_FATAL << "parallel: return UNSAT";
            return;
        }
    }
    // DREAL_LOG_FATAL << "Intersection is nonempty exiting...";
    for (unsigned i = 0; i < m_vec.size(); i++) {
        // if (threads[i].joinable()) {
        threads[i].join();
        // }
    }
    // DREAL_LOG_FATAL << "parallel: return SAT";
    return;
}
ostream & contractor_parallel_all::display(ostream & out) const {
    out << "contractor_parallel_all(";
    for (contractor const & c : m_vec) {
        out << c << ", ";
    }
    out << ")";
    return out;
}

contractor mk_contractor_parallel_all(initializer_list<contractor> const & l) {
    return contractor(make_shared<contractor_parallel_all>(l));
}
contractor mk_contractor_parallel_all(vector<contractor> const & v) {
    return contractor(make_shared<contractor_parallel_all>(v));
}
contractor mk_contractor_parallel_all(contractor const & c1, contractor const & c2) {
    return contractor(make_shared<contractor_parallel_all>(c1, c2));
}

}  // namespace dreal
