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

#include <set>
#include <unordered_map>
#include <vector>
#include "opensmt/egraph/Enode.h"
#include "util/box.h"
#include "icp/gbox.h"

namespace dreal {
class point_map {
private:
    std::unordered_map<Enode *, std::set<double>> m_map;
    void add_intv(Enode * const v, double const l, double const u) {
        assert(l <= u);
        m_map[v].insert(l);
        if (l < u) {
            m_map[v].insert(u);
        }
    }

public:
    typedef typename std::unordered_map<Enode *, std::set<double>>::iterator iterator;
    typedef typename std::unordered_map<Enode *, std::set<double>>::const_iterator const_iterator;

    void add(gbox const & gb);
    void add(std::vector<gbox> const & v);
    std::set<double> get_set(Enode * const v) const {
        return m_map.at(v);
    }
    unsigned size() const { return m_map.size(); }
    iterator begin() { return m_map.begin(); }
    const_iterator begin() const { return m_map.begin(); }
    const_iterator cbegin() const { return m_map.cbegin(); }
    iterator end() { return m_map.end(); }
    const_iterator end() const { return m_map.end(); }
    const_iterator cend() const { return m_map.cend(); }

    bool has_key(Enode * const v) const { return m_map.find(v) != m_map.end(); }
    std::set<double> const & at(Enode * const v) const { return m_map.at(v); }
    std::set<double> at(Enode * const v) { return m_map.at(v); }
};
}  // namespace dreal
