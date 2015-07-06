/*********************************************************************
Author: Soonho Kong <soonhok@cs.cmu.edu>
        Sicun Gao <sicung@cs.cmu.edu>
        Edmund Clarke <emc@cs.cmu.edu>

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

#include <vector>
#include <string>
#include <algorithm>
#include <iterator>
#include <unordered_map>
#include <unordered_set>
#include <initializer_list>
#include <iostream>
#include <utility>
#include "opensmt/egraph/Enode.h"
#include "util/constraint.h"
#include "util/flow.h"
#include "util/ibex_enode.h"
#include "util/logging.h"

using std::back_inserter;
using std::cerr;
using std::copy;
using std::endl;
using std::initializer_list;
using std::ostream;
using std::string;
using std::unordered_map;
using std::unordered_set;
using std::vector;
using std::to_string;

namespace dreal {

ostream & operator<<(ostream & out, constraint_type const & ty) {
    switch (ty) {
    case constraint_type::Nonlinear:
        out << "Nonlinear";
        break;
    case constraint_type::ODE:
        out << "ODE";
        break;
    case constraint_type::Integral:
        out << "Integral";
        break;
    case constraint_type::ForallT:
        out << "ForallT";
        break;
    case constraint_type::Forall:
        out << "Forall";
        break;
    case constraint_type::Exists:
        out << "Exists";
        break;
    }
    return out;
}


// ====================================================
// constraint
// ====================================================
constraint::constraint(constraint_type ty) : m_type(ty) {
}
constraint::constraint(constraint_type ty, Enode * const e)
    : m_type(ty), m_enodes(1, e), m_vars(e->get_vars()) {
}

unordered_set<Enode *> build_vars_from_enodes(initializer_list<vector<Enode *>> const & enodes_list) {
    unordered_set<Enode *> ret;
    for (auto const & enodes : enodes_list) {
        for (auto const e : enodes) {
            unordered_set<Enode *> const & s = e->get_vars();
            ret.insert(s.begin(), s.end());
        }
    }
    return ret;
}

constraint::constraint(constraint_type ty, vector<Enode *> const & enodes)
    : m_type(ty), m_enodes(enodes), m_vars(build_vars_from_enodes({enodes})) {
}
constraint::constraint(constraint_type ty, vector<Enode *> const & enodes_1, vector<Enode *> const & enodes_2)
    : m_type(ty), m_enodes(enodes_1), m_vars(build_vars_from_enodes({enodes_1, enodes_2})) {
    copy(enodes_2.begin(), enodes_2.end(), back_inserter(m_enodes));
}
ostream & operator<<(ostream & out, constraint const & c) {
    return c.display(out);
}

// ====================================================
// Nonlinear constraint
// ====================================================
nonlinear_constraint::nonlinear_constraint(Enode * const e, lbool p, std::unordered_map<Enode*, double> const & subst)
    : constraint(constraint_type::Nonlinear, e), m_exprctr(nullptr), m_numctr(nullptr), m_numctr_ineq(nullptr), m_subst(subst) {
    unordered_map<string, ibex::Variable const> var_map;
    bool is_ineq = (p == l_False && e->isEq());
    p = is_ineq ? true : p;

    m_exprctr = translate_enode_to_exprctr(var_map, e, p, m_subst);
    assert(m_exprctr);

    m_var_array.resize(var_map.size());
    unsigned i = 0;
    for (auto const p : var_map) {
        m_var_array.set_ref(i, p.second);
        i++;
    }

    if (is_ineq) {
        m_numctr_ineq = new ibex::NumConstraint(m_var_array, *m_exprctr);
    } else {
        m_numctr = new ibex::NumConstraint(m_var_array, *m_exprctr);
    }
    DREAL_LOG_INFO << "nonlinear_constraint: "<< *this;
}

nonlinear_constraint::~nonlinear_constraint() noexcept {
    delete m_numctr;
    delete m_numctr_ineq;
    delete m_exprctr;
}
ostream & nonlinear_constraint::display(ostream & out) const {
    out << "nonlinear_constraint ";
    if (m_numctr) {
        out << *m_numctr;
    } else {
        out << "!(" << *m_numctr_ineq << ")";
    }
    return out;
}

pair<lbool, ibex::Interval> nonlinear_constraint::eval(ibex::IntervalVector const & iv) const {
    lbool sat = l_Undef;
    ibex::Interval result;
    if (m_numctr) {
        result = m_numctr->f.eval(iv);
        switch (m_numctr->op) {
            case ibex::LT:
                if (result.ub() < 0) {
                    //    [         ]
                    // --------------- 0 --------------
                    sat = l_True;
                } else if (0 <= result.lb()) {
                    //                 [         ]
                    // --------------- 0 --------------
                    sat = l_False;
                }
                break;
            case ibex::LEQ:
                if (result.ub() <= 0) {
                    //       [         ]
                    // --------------- 0 --------------
                    sat = l_True;
                } else if (0 < result.lb()) {
                    //                   [         ]
                    // --------------- 0 --------------
                    sat = l_False;
                }
                break;
            case ibex::GT:
                if (0 < result.lb()) {
                    //                   [         ]
                    // --------------- 0 --------------
                    sat = l_True;
                } else if (result.ub() <= 0) {
                    //       [         ]
                    // --------------- 0 --------------
                    sat = l_False;
                }
                break;
            case ibex::GEQ:
                if (0 <= result.lb()) {
                    //                 [         ]
                    // --------------- 0 --------------
                    sat = l_True;
                } else if (result.ub() < 0) {
                    //     [         ]
                    // --------------- 0 --------------
                    sat = l_False;
                }
                break;
            case ibex::EQ:
                if (result.lb() == 0 && result.ub() == 0) {
                    //                 x
                    // --------------- 0 --------------
                    sat = l_True;
                } else if ((result.ub() < 0) || (0 < result.lb())) {
                    //  [            ] | [            ]
                    // --------------- 0 --------------
                    sat = l_False;
                }
                break;
        }
    } else {
        // Ineq case: lhs - rhs != 0
        assert(m_numctr_ineq);
        result = m_numctr_ineq->f.eval(iv);
        if ((result.lb() == 0) && (result.ub() == 0)) {
            //     [      lhs      ]
            //     [      rhs      ]
            sat = l_False;
        } else {
            sat = l_True;
        }
    }
    if (DREAL_LOG_DEBUG_IS_ON && sat == l_False) {
        DREAL_LOG_DEBUG << "nonlinear_constraint::eval: unsat detected";
        DREAL_LOG_DEBUG << "\t" << *this;
        DREAL_LOG_DEBUG << "input: " << iv;
        DREAL_LOG_DEBUG << "output:" << result;
    }
    return make_pair(sat, result);
}

pair<lbool, ibex::Interval> nonlinear_constraint::eval(box const & b) const {
    // Construct iv from box b
    ibex::IntervalVector iv(m_var_array.size());
    for (int i = 0; i < m_var_array.size(); i++) {
        iv[i] = b[m_var_array[i].name];
        DREAL_LOG_DEBUG << m_var_array[i].name << " = " << iv[i];
    }
    return eval(iv);
}

// ====================================================
// ODE constraint
// ====================================================
ode_constraint::ode_constraint(integral_constraint const & integral, vector<forallt_constraint> const & invs)
    : constraint(constraint_type::ODE, integral.get_enodes()), m_int(integral), m_invs(invs) {
    for (auto const & inv : invs) {
        copy(inv.get_enodes().begin(), inv.get_enodes().end(), back_inserter(m_enodes));
    }
}
ostream & ode_constraint::display(ostream & out) const {
    out << "ode_constraint" << endl;
    out << m_int << endl;
    for (forallt_constraint const & inv : m_invs) {
        out << inv << endl;
    }
    return out;
}

// ====================================================
// Integral constraint
// ====================================================
integral_constraint mk_integral_constraint(Enode * const e, unordered_map<string, flow> const & flow_map) {
    // nra_solver::inform: (integral 2 0 time_9 v_9_0 v_9_t x_9_0 x_9_t)
    Enode const * tmp = e->getCdr();
    unsigned flow_id = tmp->getCar()->getValue();
    tmp = tmp->getCdr();

    Enode * const time_0 = tmp->getCar();
    tmp = tmp->getCdr();

    Enode * const time_t = tmp->getCar();
    tmp = tmp->getCdr();

    string key = string("flow_") + to_string(flow_id);
    auto const it = flow_map.find(key);
    if (it == flow_map.end()) {
        throw std::logic_error(key + " is not in flow_map. Failed to create integral constraint");
    }
    flow const & _flow = it->second;
    vector<string> const & flow_vars = _flow.get_vars();
    vector<Enode *> const & flow_odes = _flow.get_odes();
    vector<Enode *> vars_0, vars_t, pars_0, pars_t;
    vector<string> par_lhs_names;
    vector<pair<string, Enode *>> odes;

    for (unsigned i = 0; i < flow_vars.size(); i++) {
        Enode * const var_0 = tmp->getCar();
        tmp = tmp->getCdr();
        Enode * const var_t = tmp->getCar();
        tmp = tmp->getCdr();
        string const & ode_var = flow_vars[i];
        Enode * const ode_rhs  = flow_odes[i];

        if (ode_rhs->isConstant() && ode_rhs->getValue() == 0.0) {
            // Parameter
            pars_0.push_back(var_0);
            pars_t.push_back(var_t);
            par_lhs_names.push_back(ode_var);
        } else {
            // Variable
            vars_0.push_back(var_0);
            vars_t.push_back(var_t);
            odes.emplace_back(ode_var, ode_rhs);
        }
    }
    assert(tmp->isEnil());

    return integral_constraint(e, flow_id, time_0, time_t, vars_0, pars_0, vars_t, pars_t, par_lhs_names, odes);
}

integral_constraint::integral_constraint(Enode * const e, unsigned const flow_id, Enode * const time_0, Enode * const time_t,
                                         vector<Enode *> const & vars_0, vector<Enode *> const & pars_0,
                                         vector<Enode *> const & vars_t, vector<Enode *> const & pars_t,
                                         vector<string> const & par_lhs_names,
                                         vector<pair<string, Enode *>> const & odes)
    : constraint(constraint_type::Integral, e),
      m_flow_id(flow_id), m_time_0(time_0), m_time_t(time_t),
      m_vars_0(vars_0), m_pars_0(pars_0), m_vars_t(vars_t), m_pars_t(pars_t),
      m_par_lhs_names(par_lhs_names), m_odes(odes) { }
ostream & integral_constraint::display(ostream & out) const {
    out << "integral_constraint = " << m_enodes[0] << endl;
    out << "\t" << "flow_id = " << m_flow_id << endl;
    out << "\t" << "time = [" << m_time_0 << "," << m_time_t << "]" << endl;
    for (Enode * par_0 : m_pars_0) { out << "\t" << "par_0 : " << par_0 << endl; }
    for (Enode * par_t : m_pars_t) { out << "\t" << "par_t : " << par_t << endl; }
    for (Enode * var_0 : m_vars_0) { out << "\t" << "var_0 : " << var_0 << endl; }
    for (Enode * var_t : m_vars_t) { out << "\t" << "var_t : " << var_t << endl; }
    for (auto const & ode : m_odes) {
        out << "\t" << "d/dt[" << ode.first << "] = " << ode.second << endl;
    }
    return out;
}


forallt_constraint mk_forallt_constraint(Enode * const e) {
    // nra_solver::inform: (forallt 2 0 time_9 (<= 0 v_9_t))
    Enode const * tmp = e->getCdr();
    unsigned const flow_id = tmp->getCar()->getValue();
    tmp = tmp->getCdr();

    Enode * const time_0 = tmp->getCar();
    tmp = tmp->getCdr();

    Enode * const time_t = tmp->getCar();
    tmp = tmp->getCdr();

    Enode * const inv = tmp->getCar();
    return forallt_constraint(e, flow_id, time_0, time_t, inv);
}

// ====================================================
// ForallT constraint
// ====================================================
forallt_constraint::forallt_constraint(Enode * const e, unsigned const flow_id, Enode * const time_0, Enode * const time_t, Enode * const inv)
    : constraint(constraint_type::ForallT, e), m_flow_id(flow_id), m_time_0(time_0), m_time_t(time_t), m_inv(inv) {
}
ostream & forallt_constraint::display(ostream & out) const {
    out << "forallt_constraint = " << m_enodes[0] << endl;
    out << "\t" << "flow_id = " << m_flow_id << endl;
    out << "\t" << "time = [" << m_time_0 << "," << m_time_t << "]" << endl;
    out << "\t" << "inv : " << m_inv << endl;
    return out;
}

// ====================================================
// Forall constraint
// ====================================================
forall_constraint::forall_constraint(Enode * const e, lbool const p)
    : constraint(constraint_type::Forall, e), m_forall_vars(e->get_forall_vars()), m_polarity(p) {
}
unordered_set<Enode *> forall_constraint::get_forall_vars() const {
    return m_forall_vars;
}
ostream & forall_constraint::display(ostream & out) const {
    out << "forall_constraint = "
        << (m_polarity == l_True ? "" : "!")
        << m_enodes[0] << endl;
    for (Enode * const var : m_forall_vars) {
        out << "quantified var = " << var << endl;
    }
    return out;
}

}  // namespace dreal
