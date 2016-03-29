#ifndef __ALG_INCLUDED__   // if x.h hasn't been included yet...
#define __ALG_INCLUDED__ 
#include <vector>

typedef std::vector<double> vec;
typedef double (*function_2var)(double, double);

const double X = 0.25 * M_PI;
const double T = 1.0;

class TridiagonalMatrix
{
public:
    TridiagonalMatrix(const vec& as, const vec& bs, const vec& cs):
        as(as), bs(bs), cs(cs), size(as.size()) {}

    vec as;
    vec bs;
    vec cs;
    int size;
};

int compute_linear_system_size(int actual_size);
TridiagonalMatrix construct_matrix(int nx, int nt);
vec construct_rhs(int nx, int nt, const vec& previous_layer, double t,
        function_2var f, function_2var u);
double exact_solution(double x, double t);
double f(double x, double t);
double compute_error(const vec& last_layer, int nx);

#endif 