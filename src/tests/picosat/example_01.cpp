/*********************************************************************
Author: Soonho Kong <soonhok@cs.cmu.edu>

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

#include <iostream>
#include "picosat/picosat.h"

int main() {
    PicoSAT * ps = picosat_init();
    picosat_save_original_clauses(ps);
    picosat_add(ps, 99); picosat_add(ps,  0);

    // Add: -1 -2 3 0
    std::cerr << "Add: -1 -2 3 4 0" << std::endl;
    picosat_add(ps, -1);
    picosat_add(ps, -2);
    picosat_add(ps,  3);
    picosat_add(ps,  4);
    picosat_add(ps,  0);
    picosat_print(ps, stderr);
    std::cerr << std::endl << std::endl;

    // Add: -1 0
    // Add: -2 0
    std::cerr << "Add: -1 -2 0" << std::endl;
    picosat_add(ps, -1); picosat_add(ps, 0);
    picosat_add(ps, -2); picosat_add(ps, 0);
    picosat_print(ps, stderr);
    std::cerr << std::endl << std::endl;

    picosat_sat(ps, -1);
    picosat_print(ps, stderr);
    for (int i = 1; i <= 4; ++i) {
        std::cerr << " i = \t"
                  << picosat_deref(ps, i) << "\t"
                  << picosat_deref_partial(ps, i) << std::endl;
    }
    std::cerr << std::endl << std::endl;

    std::cerr << "simplify...";
    picosat_simplify(ps);
    std::cerr << "done\n";

    picosat_sat(ps, -1);
    picosat_print(ps, stderr);
    for (int i = 1; i <= 4; ++i) {
        std::cerr << " i = \t"
                  << picosat_deref(ps, i) << "\t"
                  << picosat_deref_partial(ps, i) << std::endl;
    }
    std::cerr << std::endl << std::endl;
    return 0;
}
