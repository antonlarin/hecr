#define _USE_MATH_DEFINES 
#include <math.h>
#include <vector>
#include <chrono>
#include <fstream>
#include  "tbb.hpp"

vec cyclic_reduction(const TridiagonalMatrix& m, const vec& b);
double solve_problem(int nx, int nt);

