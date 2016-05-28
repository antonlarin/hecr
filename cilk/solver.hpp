#include <cmath>
#include <vector>

#include "algorithms.hpp"

class Solver
{
public:
	vec cyclic_reduction(const TridiagonalMatrix& m, const vec& b);
	double solve_problem(int nx, int nt);
        static vec construct_rhs(int nx, int nt, const vec& previous_layer,
                double t, function_2var f, function_2var u);
};

