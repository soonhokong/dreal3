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

#pragma once

#include <chrono>
#include <iostream>

namespace dop {

class stat {
private:
    std::chrono::time_point<std::chrono::steady_clock> m_start;

    std::ostream & m_out;
    unsigned m_num_vars;
    unsigned m_num_overapprox_CEs;
    unsigned m_num_underapprox_CEs;
    unsigned m_num_branches;


public:
    stat(std::ostream & out);
    void print() const;
    unsigned get_num_vars() const { return m_num_vars; }
    unsigned get_num_overapprox_CEs() const { return m_num_overapprox_CEs; }
    unsigned get_num_underapprox_CEs() const { return m_num_underapprox_CEs; }
    unsigned get_num_branches() const { return m_num_branches; }

    void set_num_vars(unsigned n) { m_num_vars = n; }
    void set_num_overapprox_CEs(unsigned n) { m_num_overapprox_CEs = n; }
    void set_num_underapprox_CEs(unsigned n) { m_num_underapprox_CEs = n; }
    void set_num_branches(unsigned n) { m_num_branches = n; }
};

}  // namespace dop
