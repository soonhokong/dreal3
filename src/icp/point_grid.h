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

	//derived clauses 
	std::vector<std::vector<int> *> linear_clauses; //linear ordering on bounds
	std::vector<std::vector<int> *> lu_clauses; //lower and upper bounds should be consistent

public:
	Grid(box const &);
	~Grid();
	void add_box(box const &);
	void add_point(Enode * v, double const p);
};
}
