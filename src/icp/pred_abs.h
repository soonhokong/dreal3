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
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include "opensmt/egraph/Enode.h"
#include "util/hash_combine.h"

namespace std {
template<>
struct hash<std::tuple<Enode *, double, bool>> {
    size_t operator () (const std::tuple<Enode *, double, bool> & v) const {
        std::size_t s = 23;
        dreal::hash_combine<Enode *>(s, std::get<0>(v));
        dreal::hash_combine<double>(s, std::get<1>(v));
        dreal::hash_combine<bool>(s, std::get<2>(v));
        return s;
    }
};
template<>
struct equal_to<std::tuple<Enode *, double, bool>> {
    bool operator() (const std::tuple<Enode *, double, bool> & v1, const std::tuple<Enode *, double, bool> & v2) const {
        return std::get<0>(v1) == std::get<0>(v2)
            && std::get<1>(v1) == std::get<1>(v2)
            && std::get<2>(v1) == std::get<2>(v2);
    }
};
}  // namespace std

namespace dreal {
class pred_abs {
private:
    int m_num_vars = 0;
    std::unordered_map<int, std::tuple<Enode *, double, bool>> m_con_map;
    std::unordered_map<std::tuple<Enode *, double, bool>, int> m_abs_map;
    int add(Enode * const v, double const bound, bool const le);

public:
    int add_le(Enode * const v, double const bound) { return add(v, bound, true); }
    int add_le(double const bound, Enode * const v) { return add_ge(v, bound); }
    int add_ge(Enode * const v, double const bound) { return add(v, bound, false); }
    int add_ge(double const bound, Enode * const v) { return add_le(v, bound); }
    int add_lt(Enode * const v, double const bound) { return -add_ge(v, bound); }
    int add_lt(double const bound, Enode * const v) { return -add_ge(bound, v); }
    int add_gt(Enode * const v, double const bound) { return -add_le(v, bound); }
    int add_gt(double const bound, Enode * const v) { return -add_le(bound, v); }
    int get_num_vars() const { return m_num_vars; }
    int lookup(Enode * const v, double const bound, bool const le) const { return m_abs_map.at(std::make_tuple(v, bound, le)); }
    std::tuple<Enode *, double, bool> lookup(int const n) const { return m_con_map.at(n); }
    void debug_print() const;
};
}  // namespace dreal
