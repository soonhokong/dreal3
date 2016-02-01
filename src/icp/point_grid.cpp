/*********************************************************************
Author: Sicun Gao <sicung@cs.cmu.edu>
        Soonho Kong <soonhok@cs.cmu.edu>

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

#include <initializer_list>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include "icp/point_grid.h"

using std::cerr;
using std::endl;
using std::get;
using std::initializer_list;
using std::make_pair;
using std::set;
using std::tuple;
using std::unordered_set;
using std::vector;

namespace dreal {
    Grid::Grid(box const & b) : full_lb_clauses({0}), full_ub_clauses({0}) {
        top_lit = 0;
        auto const & vars = b.get_vars();
        for (unsigned i = 0; i < b.size(); ++i) {
            //initialize storage
            Enode * v = vars[i];
            point_rows.emplace(v, initializer_list<double>{});
            // vector<int> vec;
            lb_lits.emplace(v, initializer_list<int>{});
            ub_lits.emplace(v, initializer_list<int>{});
            //add the points
            add_point(v,b[v].lb());
            add_point(v,b[v].ub());
        }
    }

    void Grid::add_box(box const & b) {
        //add the endpoints on each dimension to its point row
        auto const & vars = b.get_vars();
        for (unsigned i = 0; i < b.size(); ++i) {
            Enode * v = vars[i];
            add_point(v,b[v].lb());
            add_point(v,b[v].ub());
        }
    }

//generate a new literal for a point at position it
    void Grid::add_point(Enode * v, double const p) {
        assert(v);
        set<double> & row_v = point_rows[v];
        unsigned old_size = row_v.size();
        set<double>::iterator it = row_v.emplace(p).first;

        if ( old_size == row_v.size() ) return;

        //I'm gonna move around on the sequence and it became buggy.
        //so I'm temporarily putting things in a vector. can definitely be improved.
        vector<double> row_helper;
        int id = 0;
        //copy the first part of the vector
        for (set<double>::iterator itt = row_v.begin(); itt != it; itt++ ) {
            row_helper.push_back(*itt);
            id++;
        }//exits when id is for the new point
        //copy the rest of the set to the vector
        for (set<double>::iterator itt = it; itt!=row_v.end(); itt++ ) {
            row_helper.push_back(*itt);
        }

        //get the lists of lb literals and ub literals for v
        vector<int> & lb_list = lb_lits[v];
        vector<int> & ub_list = ub_lits[v];

        //new point always adds a new lower bound literal and a new upper bound literal
        int new_lb_lit = ++top_lit;
        lb_list.push_back(new_lb_lit); //lb lit is odd
        lb_lit_map.emplace(make_pair(v,p),new_lb_lit);

        int new_ub_lit = ++top_lit;
        ub_list.push_back(new_ub_lit); //ub lit is even
        assert(new_ub_lit == new_lb_lit + 1);//sanity check
        assert(top_lit%2 == 0); //should be even

        //update the lb and ub disjunctive constraints. note that 0 is poped first and then pushed back
        assert(full_lb_clauses.size() > 0);
        full_lb_clauses.pop_back();
        full_lb_clauses.push_back(new_lb_lit);
        full_lb_clauses.push_back(0);

        assert(full_ub_clauses.size() > 0);
        full_ub_clauses.pop_back();
        full_ub_clauses.push_back(new_ub_lit);
        full_ub_clauses.push_back(0);

        //clean up the cache for pushed_clauses
        push_linear_clauses.clear();
        push_lu_clauses.clear();

        /*compare the new point with its neighbors*/
        //if it's not the smallest, add constraints regarding the left side
        int pre = id-1;
        if (pre > 0) {
            //value of the previous point
            double pre_value = row_helper[pre];
            //index of the previous lb literal
            int pre_index = lb_lit_map[make_pair(v,pre_value)];

            //first, add a new linear clause for "x>it \implies x>it--"
            vector<int> lc_left_low;
            lc_left_low.push_back( -new_lb_lit ); //negation of x>it
            lc_left_low.push_back( pre_index ); //x>it--
            lc_left_low.push_back( 0 ); //end of clause

            linear_clauses.push_back(lc_left_low); //put it in the global storage
            push_linear_clauses.push_back(lc_left_low);

            //next, add a new linear clause for "x<it-- \implies x<it"
            //WARNING: note that all ub literals requires +1 on the index
            vector<int> lc_left_up;
            lc_left_up.push_back( - (pre_index+1) ); //negation of x<it--
            lc_left_up.push_back( new_ub_lit ); //x>it--
            lc_left_up.push_back( 0 ); //end of clause

            linear_clauses.push_back(lc_left_up); //put it in the global storage
            push_linear_clauses.push_back(lc_left_up);

            //now add a new clause saying "u_i \implies \not l_i-1"
            vector<int> lu_left;
            lu_left.push_back( - new_ub_lit );
            lu_left.push_back( - pre_index );
            lu_left.push_back( 0 );

            lu_clauses.push_back(lu_left);
            push_lu_clauses.push_back(lu_left);

        }//else, do nothing

        unsigned succ = id + 1;
        if (succ < row_helper.size()) {
            //value of the next point
            double succ_value = row_helper[succ];
            //index of the next lb literal
            int succ_index = lb_lit_map[make_pair(v,succ_value)];

            //first, add a new linear clause for "x>it++ \implies x>it"
            vector<int> lc_right_low;
            lc_right_low.push_back( -succ_index ); //negation of x>it
            lc_right_low.push_back( new_lb_lit ); //x>it--
            lc_right_low.push_back( 0 ); //end of clause

            linear_clauses.push_back(lc_right_low); //put it in the global storage
            push_linear_clauses.push_back(lc_right_low);

            //next, add a new linear clause for "x<it \implies x<it++"
            //WARNING: note that all ub literals requires +1 on the index
            vector<int> lc_right_up;
            lc_right_up.push_back( - new_ub_lit ); //negation of x<it
            lc_right_up.push_back( succ_index+1 ); //x>it++
            lc_right_up.push_back( 0 ); //end of clause

            linear_clauses.push_back(lc_right_up); //put it in the global storage
            push_linear_clauses.push_back(lc_right_up);

            //now add a new clause saying "l_i \implies \not u_i+1"
            vector<int> lu_right;
            lu_right.push_back( -new_lb_lit );
            lu_right.push_back( - (succ_index+1) );
            lu_right.push_back( 0 );

            lu_clauses.push_back(lu_right);
            push_lu_clauses.push_back(lu_right);

        }//else, do nothing
    }

    void Grid::debug_print_clause(vector<int> const & c) const {
        for (int const l : c) {
            cerr << l << " ";
        }
        cerr << 0 << endl;
    }

    void Grid::debug_print() const {
        cerr << "======== POINT ROWS ===========" << endl;
        for (auto it : point_rows) {
            cerr << it.first << " : ";
            for (double const c : it.second) {
                cerr << c << " ";
            }
            cerr << endl;
        }
        cerr << "======== VAR ENCODING ===========" << endl;
        for (auto it : point_rows) {
            Enode * v = it.first;
            for (double const c : it.second) {
                int const lb_lit = lb_lit_map.at(make_pair(v, c));
                int const ub_lit = lb_lit + 1;
                cerr << v << " <= " << c << " : " << lb_lit << endl;
                cerr << v << " >= " << c << " : " << ub_lit << endl;
            }
        }
        cerr << "======== LINEAR CLAUSES ===========" << endl;
        for (auto const & c : linear_clauses) {
            debug_print_clause(c);
        }
        cerr << "======== PUSH LINEAR CLAUSES ===========" << endl;
        for (auto const & c : push_linear_clauses) {
            debug_print_clause(c);
        }
        cerr << "======== LU CLAUSES ===========" << endl;
        for (auto const & c : lu_clauses) {
            debug_print_clause(c);
        }
        cerr << "======== PUSH LU CLAUSES ===========" << endl;
        for (auto const & c : push_lu_clauses) {
            debug_print_clause(c);
        }

    }

}
