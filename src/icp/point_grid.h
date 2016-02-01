#pragma once

#include <functional>
#include <iostream>
#include <limits>
#include <map>
#include <random>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include "util/logging.h"
#include "util/box.h"
#include "util/stat.h"
#include "contractor/contractor.h"
#include "opensmt/smtsolvers/SMTConfig.h"
#include "picosat/picosat.h"

namespace dreal {

class Grid { //Point Matrix
private:
	//top_lit keeps the index of the last literal
	unsigned top_lit;

	//we need a map from a variable to the ordered points in its dimension
	std::unordered_map< Enode*, std::set<double>* > point_rows;

	//encoding: lower bound constraints x>c use odd numbers. upper bounds use even numbers.
	std::unordered_map< Enode*, std::vector<int>* >	lb_lits;
	std::unordered_map< Enode*, std::vector<int>* >	ub_lits;

	/* a global mapping from (var, point) to lit index. 
	lb and ub are automatically inferred from whether the index is odd or even
	by default it returns the lb constraint. Use +1 for the ub constraint. */
	std::map< std::pair<Enode*, double>, int > lb_lit_map;

	//clauses 
	std::vector<std::vector<int> *> linear_clauses; //linear ordering on bounds
	std::vector<std::vector<int> *> lu_clauses; //lower and upper bounds should be consistent
	std::vector<int>		full_lb_clauses;
	std::vector<int>		full_ub_clauses;	

	//new clauses that should be pushed
	std::vector<std::vector<int> *> push_linear_clauses;
	std::vector<std::vector<int> *> push_lu_clauses; 

	//prepare for sat
	std::vector<int> * current_formula;
	std::vector<int> * push_formula;
public:
	Grid(box const &);
	~Grid();
	void add_box(box const &);
	void add_point(Enode * v, double const p);
	inline void build_current_formula() {
		current_formula->clear();
		for (auto cl : linear_clauses)
			for (auto i : *cl)
				current_formula->push_back(i);
		for (auto cl : lu_clauses)  
			for (auto i : *cl) 
				current_formula->push_back(i);
		for (auto i : full_lb_clauses)	
			current_formula->push_back(i);
		for (auto i : full_ub_clauses)	
			current_formula->push_back(i);
	}
	inline void build_push_formula() {
		push_formula->clear();
		for (auto cl : push_linear_clauses)
			for (auto i : *cl)
				push_formula->push_back(i);
		for (auto cl : push_lu_clauses)  
			for (auto i : *cl) 
				push_formula->push_back(i);
		for (auto i : full_lb_clauses)	
			push_formula->push_back(i);
		for (auto i : full_ub_clauses)	
			push_formula->push_back(i);
	}
	inline void build_push_nobounds_formula() {
		push_formula->clear();
		for (auto cl : push_linear_clauses)
			for (auto i : *cl)
				push_formula->push_back(i);
		for (auto cl : push_lu_clauses)  
			for (auto i : *cl) 
				push_formula->push_back(i);
	}
	inline void build_push_bounds_only_formula() {
		push_formula->clear();
		for (auto i : full_lb_clauses)	
			push_formula->push_back(i);
		for (auto i : full_ub_clauses)	
			push_formula->push_back(i);
	}

	inline std::vector<int> * get_current_formula() { 
		build_current_formula();
		return current_formula; 
	}
	inline std::vector<int> * get_push_formula() { 
		build_push_formula();
		return push_formula; 
	}
	inline std::vector<int> * get_push_nobounds_formula() { 
		build_push_nobounds_formula();
		return push_formula; 
	}
	inline std::vector<int> * get_push_bounds_only_formula() { 
		build_push_bounds_only_formula();
		return push_formula; 
	}
};
}
