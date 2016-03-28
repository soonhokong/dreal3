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
#include <iostream>
#include <vector>
#include "util/hash_combine.h"
#include "util/box.h"
#include "icp/gbox.h"

namespace dreal {

class clause {
public:
    gbox m_antecedent;
    std::vector<gbox> m_consequent;

public:
    clause(gbox const & antecedent, gbox const & consequent)
        : m_antecedent(antecedent), m_consequent(1, consequent) { }
    clause(gbox const & antecedent, gbox const & consequent_1, gbox const & consequent_2)
        : m_antecedent(antecedent), m_consequent({consequent_1, consequent_2}) { }
    friend std::ostream & operator<<(std::ostream & out, clause const & c);
};
std::ostream & operator<<(std::ostream & out, clause const & c);

}  // namespace dreal


namespace std {
template<>
struct hash<dreal::clause> {
    size_t operator () (const dreal::clause & v) const {
        std::size_t s = 23;
        dreal::hash_combine<dreal::gbox>(s, v.m_antecedent);
        for (dreal::gbox const & b : v.m_consequent) {
            dreal::hash_combine<dreal::gbox>(s, b);
        }
        return s;
    }
};
template<>
struct equal_to<dreal::clause> {
    bool operator() (const dreal::clause & v1, const dreal::clause & v2) const {
        return v1.m_antecedent == v2.m_antecedent &&
            v1.m_consequent == v2.m_consequent;
    }
};
}  // namespace std
