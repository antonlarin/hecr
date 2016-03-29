#include "algorithms.hpp"
#include <vector>

class Solver
{
public:

	virtual vec cyclic_reduction(const TridiagonalMatrix& m, const vec& b);
	virtual double solve_problem(int nx, int nt);
};

class SolverOmp : Solver
{
public:

	vec cyclic_reduction(const TridiagonalMatrix& m, const vec& b);
	double solve_problem(int nx, int nt);
};
