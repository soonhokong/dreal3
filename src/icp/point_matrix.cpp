#include <unordered_set>
#include <unordered_map>
#include <tuple>
#include "icp/picosat_wrapper.h"
#include "icp/point_matrix.h"

using std::unordered_set;
using std::tuple;
using std::get;
using std::cerr;
using std::endl;
using std::vector;
using std::set;
using std::make_pair;

namespace dreal {

PtMatrix::PtMatrix(box const & b) {
	top_lit = 0;
	auto const & vars = b.get_vars();
        for (unsigned i = 0; i < b.size(); ++i) {
		//initialize storage
		Enode * v = vars[i];
		set<double> * row_v = new set<double>;
		vector<int> * lb_v = new vector<int>;
		vector<int> * ub_v = new vector<int>;
		point_rows.emplace(v,row_v);
		lb_lits.emplace(v,lb_v);
		ub_lits.emplace(v,ub_v);
		//add the point
		add_point(v,b[v].lb());
		add_point(v,b[v].ub());
        }
}

PtMatrix::~PtMatrix() {
	for (auto const & it : point_rows) delete it.second;
	for (auto const & it : lb_lits) delete it.second;
	for (auto const & it : ub_lits) delete it.second;
	for (auto const & it : linear_clauses) delete it;
	for (auto const & it : lu_clauses) delete it;
}

void PtMatrix::add_box(box const & b) {
	//add the endpoints on each dimension to its point row
        auto const & vars = b.get_vars();
        for (unsigned i = 0; i < b.size(); ++i) {
        	Enode * v = vars[i];
		add_point(v,b[v].lb());
		add_point(v,b[v].ub());
        }
}

//generate a new literal for a point at position it
void PtMatrix::add_point(Enode * v, double const p) {

	assert(v);

	set<double> * row_v = point_rows[v];
	unsigned old_size = row_v -> size();
	set<double>::iterator it = (row_v->emplace(p)).first;

	if (old_size == row_v -> size()) return;
	
	//I'm gonna move around on the sequence and it's risky to do on the set. 
	//should be improved. 
	vector<double> row_helper;	
	int id = 0;
	//copy the first part of the vector
	for (set<double>::iterator itt = row_v->begin(); itt != it; itt++ ) {
		row_helper.push_back(*itt);
		id++;
	}//exits when id is for the new point
	//copy the rest of the set to the vector
	for (set<double>::iterator itt = it; itt!=row_v->end(); itt++ ) {
		row_helper.push_back(*itt);
	}

	//get the lists of lb literals and ub literals for v
	vector<int> * lb_list = lb_lits[v];
	vector<int> * ub_list = ub_lits[v];

	//new point always adds a new lower bound literal and a new upper bound literal
	int new_lb_lit = ++top_lit;
 	lb_list -> push_back(new_lb_lit); //lb lit is odd
	lb_lit_map.emplace(std::make_pair(v,p),new_lb_lit);

	int new_ub_lit = ++top_lit;
	ub_list -> push_back(new_ub_lit); //ub lit is even
	assert(new_ub_lit == new_lb_lit + 1);//sanity check
	assert(top_lit%2 == 0); //should be even

	/*compare the new point with its neighbors*/	
	//if it's not the smallest, add constraints regarding the left side
	int pre = id-1;
	if (pre > 0) {
		//value of the previous point
		double pre_value = row_helper[pre];
		//index of the previous lb literal
		int pre_index = lb_lit_map[std::make_pair(v,pre_value)];
		//first, add a new linear clause for "x>it \implies x>it--"
		vector<int> * lc_left_low = new vector<int>;
		lc_left_low -> push_back( -new_lb_lit ); //negation of x>it
		lc_left_low -> push_back( pre_index ); //x>it--
		lc_left_low -> push_back( 0 ); //end of clause
		linear_clauses.push_back(lc_left_low); //put it in the global storage
		//next, add a new linear clause for "x<it-- \implies x<it"
		//WARNING: note that all ub literals requires +1 on the index
		vector<int> * lc_left_up = new vector<int>;
		lc_left_up -> push_back( - (pre_index+1) ); //negation of x<it--
		lc_left_up -> push_back( new_ub_lit ); //x>it--
		lc_left_up -> push_back( 0 ); //end of clause
		linear_clauses.push_back(lc_left_up); //put it in the global storage
	
		//now add a new clause saying "u_i \implies \not l_i-1"
		vector<int> * lu_left = new vector<int>;
		lu_left -> push_back( - new_ub_lit );
		lu_left -> push_back( - pre_index );
		lu_left -> push_back( 0 );
		lu_clauses.push_back( lu_left );

	}//else, do nothing

	unsigned succ = id + 1;
	if (succ < row_helper.size()) {
		//value of the next point
		double succ_value = row_helper[succ];
		//index of the next lb literal
		int succ_index = lb_lit_map[std::make_pair(v,succ_value)]; 
		//first, add a new linear clause for "x>it++ \implies x>it"
		vector<int> * lc_right_low = new vector<int>;
		lc_right_low -> push_back( -succ_index ); //negation of x>it
		lc_right_low -> push_back( new_lb_lit ); //x>it--
		lc_right_low -> push_back( 0 ); //end of clause
		linear_clauses.push_back(lc_right_low); //put it in the global storage
		//next, add a new linear clause for "x<it \implies x<it++"
		//WARNING: note that all ub literals requires +1 on the index
		vector<int> * lc_right_up = new vector<int>;
		lc_right_up -> push_back( - new_ub_lit ); //negation of x<it
		lc_right_up -> push_back( succ_index+1 ); //x>it++
		lc_right_up -> push_back( 0 ); //end of clause
		linear_clauses.push_back(lc_right_up); //put it in the global storage
	
		//now add a new clause saying "l_i \implies \not u_i+1"
		vector<int> * lu_right = new vector<int>;
		lu_right -> push_back( -new_lb_lit );
		lu_right -> push_back( - (succ_index+1) );
		lu_right -> push_back( 0 );
		lu_clauses.push_back( lu_right );

	}//else, do nothing
} 

}


