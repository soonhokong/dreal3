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

#include <unordered_set>
#include <iostream>
#include "icp/gbox.h"
#include "util/box.h"
#include "util/ibex_interval_hash.h"

namespace dreal {
class reduced_box_set {
private:
    std::unordered_set<gbox> m_set;
    void add_and_simplify(gbox const & gb);
    bool exist(gbox const & gb) const {
        return m_set.find(gb) != m_set.end();
    }

public:
    typedef typename std::unordered_set<gbox>::iterator iterator;
    typedef typename std::unordered_set<gbox>::const_iterator const_iterator;

    //  functions
    unsigned size() const { return m_set.size(); }
    iterator begin() { return m_set.begin(); }
    const_iterator begin() const { return m_set.begin(); }
    const_iterator cbegin() const { return m_set.cbegin(); }
    iterator end() { return m_set.end(); }
    const_iterator end() const { return m_set.end(); }
    const_iterator cend() const { return m_set.cend(); }

    void add(gbox const & gb);
    bool includes(gbox const & gb) const;
    bool check_invariant() const;

private:
    friend std::ostream & operator<<(std::ostream & out, reduced_box_set const & rs);
};

std::unordered_set<ibex::IntervalVector> interval_vector_diff(ibex::IntervalVector v1, ibex::IntervalVector const & v2);
std::ostream & operator<<(std::ostream & out, reduced_box_set const & rs);
}  // namespace dreal
