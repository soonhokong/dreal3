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

#include <algorithm>
#include <iomanip>
#include <memory>
#include <vector>
#include "icp/gbox.h"
#include "util/ibex_helper.h"
#include "util/logging.h"
#include "icp/reduced_box_set.h"

using std::any_of;
using std::cerr;
using std::endl;
using std::ostream;
using std::setprecision;
using std::setw;
using std::shared_ptr;
using std::size_t;
using std::unordered_set;
using std::vector;

namespace dreal {

gbox::gbox(box const & b)
    : m_vars(b.get_vars()), m_vec(b.get_values()), m_gbits(m_vec.size(), false) {
}
gbox::gbox(vector<Enode *> const & vars)
    : m_vars(vars), m_vec(m_vars.size()), m_gbits(m_vec.size(), false) {
}
gbox::gbox(vector<Enode *> const & vars, ibex::IntervalVector const & v)
    : m_vars(vars), m_vec(v), m_gbits(m_vec.size(), false) {
}
gbox::gbox(vector<Enode *> const & vars, ibex::IntervalVector const & vec, vector<bool> const & gbits)
    : m_vars(vars), m_vec(vec), m_gbits(gbits) {
}
gbox::gbox(box const & b, unordered_set<shared_ptr<constraint>> const & used_ctrs, ibex::IntervalVector const & dom)
    : m_vars(b.get_vars()), m_vec(b.get_values()), m_gbits(m_vec.size(), true) {
    for (shared_ptr<constraint> const ctr_ptr : used_ctrs) {
        for (Enode * const v : ctr_ptr->get_vars()) {
            int const i = b.get_index(v);
            m_gbits[i] = false;
        }
    }
    for (unsigned i = 0; i < size(); ++i) {
        if (m_gbits[i]) {
            m_vec[i] = dom[i];
        }
    }
}
size_t gbox::hash() const {
    std::size_t seed = 23;
    for (unsigned i = 0; i < m_vars.size(); ++i) {
        hash_combine<Enode*>(seed, m_vars[i]);
    }
    hash_combine<ibex::IntervalVector>(seed, m_vec);
    hash_combine<std::vector<bool>>(seed, m_gbits);
    return seed;
}

bool gbox::operator==(gbox const & gb) const {
    return m_vars == gb.m_vars && m_vec == gb.m_vec && m_gbits == gb.m_gbits;
}

unordered_set<gbox> gbox::set_minus(gbox const & gb) const {
    return ::dreal::set_minus(*this, gb);
}

void gbox::set_value(unsigned const i, double const lb, double const ub) {
    set_value(i, ibex::Interval(lb, ub));
}

void gbox::set_value(unsigned const i, ibex::Interval const & iv) {
    m_vec[i] = iv;
    m_gbits[i] = false;
}

// Generalize this box using another box ant.
//
// If a dimension is not specified by ant, let's not touch it.
void gbox::generalize(gbox const & ant, gbox const & con) {
    DREAL_LOG_FATAL << "Generalize";
    DREAL_LOG_FATAL << "=================";
    DREAL_LOG_FATAL << *this;
    DREAL_LOG_FATAL << "=== with(ant) ===";
    DREAL_LOG_FATAL << ant;
    DREAL_LOG_FATAL << "=== with(con) ===";
    DREAL_LOG_FATAL << con;
    DREAL_LOG_FATAL << "=== result    ===";
    assert(size() == ant.size());
    assert(size() == con.size());
    for (unsigned i = 0; i < size(); ++i) {
        if (ant.get_gbit(i)) {
            continue;
        } else {
            assert(!get_gbit(i));
            m_vec[i] = ant.m_vec[i];
        }
    }
    DREAL_LOG_FATAL << *this << endl << endl;
}

bool gbox::is_subset(gbox const & gb) const {
    assert(size() == gb.size());
    for (unsigned i = 0; i < size(); ++i) {
        if (gb.get_gbit(i)) {
#ifndef NDEBUG
            auto const & i1 = m_vec[i];
            auto const & i2 = gb.m_vec[i];
            assert(i2.lb() <= i1.lb());
            assert(i1.ub() <= i2.ub());
#endif
            continue;
        } else {
            // both of them are not generalized.
            auto const & i1 = m_vec[i];
            auto const & i2 = gb.m_vec[i];
            if (!i1.is_subset(i2)) {
                assert(i2.lb() > i1.lb() || i1.ub() > i2.ub());
                return false;
            }
            assert(i2.lb() <= i1.lb());
            assert(i1.ub() <= i2.ub());
        }
    }
    return true;
}

bool gbox::match(gbox const & ant, gbox const & con) const {
    assert(size() == con.size());
    assert(size() == ant.size());
    for (unsigned i = 0; i < size(); ++i) {
        assert(con.get_gbit(i) == ant.get_gbit(i));
        if (con.get_gbit(i)) {
            continue;
        } else {
            if (get_gbit(i)) {
                auto const & intv = m_vec[i];
                auto const & ant_intv = ant.m_vec[i];
                auto const & con_intv = con.m_vec[i];
                if (ant_intv.is_superset(intv) && intv.is_superset(con_intv)) {
                    continue;
                } else {
                    return false;
                }
            } else {
                // both of them are not generalized.
                auto const & intv = m_vec[i];
                auto const & ant_intv = ant.m_vec[i];
                auto const & con_intv = con.m_vec[i];
                if (ant_intv.is_strict_superset(intv) && intv.is_superset(con_intv)) {
                    continue;
                } else {
                    return false;
                }
            }
        }
    }
    return true;
}

// // If b1 and b2 are only different in a dim `i`
// // and b1[i] and b2[i] are adjacent and disjoint,
// // return i. Otherwise, return -1.
// int find_merge_dim(gbox const & b1, gbox const & b2) {
//     int ret = -1;
//     assert(b1.size() == b2.size());
//     for (unsigned i = 0; i < b1.size(); ++i) {
//         auto const & intv1 = b1[i];
//         auto const & intv2 = b2[i];
//         if (intv1 == intv2) { continue; }
//         if (b1.get_gbit(i) && b2.get_gbit(i)) {
//             assert(b1.get_value(i) == b2.get_value(i));
//             continue;
//         }
//         if (b1.get_gbit(i) || b2.get_gbit(i)) {
// #ifndef NDEBUG
//             auto const & i1 = b1.get_value(i);
//             auto const & i2 = b2.get_value(i);
//             assert(i1.is_subset(i2) || i2.is_subset(i1));
// #endif
//             if (ret < 0) {
//                 ret = i;
//             } else {
//                 ret = -1;
//                 break;
//             }
//         } else {
//             // b1 and b2 are not generalized at i-th dim
//             if (interval_overlaps(intv1, intv2)) {
//                 if (ret < 0) {
//                     ret = i;
//                 } else {
//                     // ret >= 0. That is, there was another dim satisfying
//                     // the condition. We conclude that these two boxes are
//                     // not mergable.
//                     ret = -1;
//                     break;
//                 }
//             } else {
//                 // there is no overlap in `b1[i]` and `b2[i]`
//                 // so they are not mergable.
//                 ret = -1;
//                 break;
//             }
//         }
//     }
//     return ret;
// }

int find_merge_dim(gbox const & b1, gbox const & b2) {
    int ret = -1;
    assert(b1.size() == b2.size());
    for (unsigned i = 0; i < b1.size(); ++i) {
        if (b1.get_gbit(i) || b2.get_gbit(i)) {
            continue;
        }
        auto const & intv1 = b1[i];
        auto const & intv2 = b2[i];
        // b1 and b2 are not generalized at i-th dim
        if (interval_overlaps(intv1, intv2)) {
            if (ret < 0) {
                ret = i;
            } else {
                // ret >= 0. That is, there was another dim satisfying
                // the condition. We conclude that these two boxes are
                // not mergable.
                ret = -1;
                break;
            }
        }
    }
    return ret;
}

gbox merge(gbox b1, gbox const & b2, int const merge_dim) {
    b1.set_value(merge_dim, b1[merge_dim] | b2[merge_dim]);
    return b1;
}

gbox generalize(gbox gb, gbox const & ant, gbox const & con) {
    gb.generalize(ant, con);
    return gb;
}

unordered_set<gbox> set_minus(gbox gb1, gbox const & gb2) {
    assert(gb1.size() == gb2.size());
    assert(gb1.is_strict_overlap(gb2));
#ifndef NDEBUG
    gbox gb1_init(gb1);  // save the initial gb1
#endif
    unordered_set<gbox> ret;
    ibex::Interval tmp1, tmp2;  // tmp value for interval_diff in the loop below
    for (unsigned i = 0; i < gb1.size(); ++i) {
        if (gb2.get_gbit(i)) {
            // we will remove the whole interval in this dimension
            continue;
        } else {
            // i-dim is not generalized in gb.
            if (gb1.get_gbit(i)) {
                // if i-dim is generalized in gb1, we need to
                // instantiaged with its domain so that it can be
                // subtracted from gb2[i].
                gb1.set_gbit(i, false);
            }
            auto const & i1 = gb1.get_value(i);
            auto const & i2 = gb2.get_value(i);
            if (i1 != i2) {
                int const diff_num = interval_diff(i1, i2, tmp1, tmp2);
                if (diff_num > 0) {
                    gb1.set_value(i, i1 & i2);
                    ibex::IntervalVector new_v(gb1.get_values());
                    new_v[i] = tmp1;
                    ret.emplace(gb1.get_vars(), new_v, gb1.get_gbits());
                    if (diff_num > 1) {
                        new_v[i] = tmp2;
                        ret.emplace(gb1.get_vars(), new_v, gb1.get_gbits());
                    }
                }
            }
        }
    }
#ifndef NDEBUG
    reduced_box_set rs;
    rs.add(intersection(gb1_init, gb2));
    for (gbox const & gb : ret) {
        rs.add(gb);
        assert(gb.is_subset(gb1_init));
        assert(!gb.is_strict_overlap(gb2));
    }
    assert(rs.size() == 1);
    gbox gb_in_rs = *rs.begin();
    if (gb_in_rs.get_values() != gb1_init.get_values()) {
        DREAL_LOG_FATAL << gb_in_rs;
        DREAL_LOG_FATAL << "=========";
        DREAL_LOG_FATAL << gb1_init;
        assert(false);
    }
#endif
    return ret;
}

gbox intersection(gbox gb1, gbox const & gb2) {
    for (unsigned i = 0; i < gb1.size(); ++i) {
        gb1.set_gbit(i, gb1.get_gbit(i) && gb2.get_gbit(i));
        gb1.set_value(i, gb1.get_value(i) & gb2.get_value(i));
    }
    return gb1;
}

ostream& operator<<(ostream& out, gbox const & gb) {
    for (int i = 0; i < gb.m_vec.size(); ++i) {
        out << setprecision(16) << "["
            << setw(20) << gb.m_vec[i].lb() << ", "
            << setw(20) << gb.m_vec[i].ub() << "]\t"
            << gb.m_gbits[i] << "\t"
            << (gb[i].lb() == gb[i].ub() ? "P" : "") << endl;
    }
    return out;
}
}  // namespace dreal
