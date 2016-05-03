#include "algorithms.hpp"

int compute_linear_system_size(int actual_size)
{
    int padded_size = 2;
    while (padded_size - 1 < actual_size)
    {
        padded_size *= 2;
    }

    return padded_size - 1;
}


TridiagonalMatrix construct_matrix(int nx, int nt)
{
    int size = compute_linear_system_size(nx - 1);

    vec as(size, 1.0);
    vec bs(size, 0.0);
    vec cs(size, 0.0);

    double h = X / nx;
    double tau = T / nt;

    double a = - 1.0 / tau - 2.0 / (h * h);
    double b_and_c = 1.0 / (h * h);

    as[0] = a;
    cs[0] = b_and_c;
    for (int i = 1; i < nx - 2; i++)
    {
        as[i] = a;
        bs[i] = b_and_c;
        cs[i] = b_and_c;
    }
    as[nx - 2] = a;
    bs[nx - 2] = b_and_c;

    return TridiagonalMatrix(as, bs, cs);
}

vec construct_rhs(int nx, int nt, const vec& previous_layer, double t,
        function_2var f, function_2var u)
{
    static int size = compute_linear_system_size(nx - 1);

    vec rhs(size, 0.0);

    static double h = X / nx;
    static double tau = T / nt;

    for (int i = 0; i < nx - 1; i++)
    {
        rhs[i] = -previous_layer[i] / tau - f(h * (i + 1), t);
    }

    rhs[0] -= u(0.0, t) / h / h;
    rhs[nx - 2] -= u(X, t) / h / h;

    return rhs;
}


double exact_solution(double x, double t)
{
    return 5 * exp(-2 * t) * sin(3 * x) - 2 * exp(-6 * t) * cos(2 * x);
}

double f(double x, double t)
{
    return 35 * exp(-2 * t) * sin(3 * x) + 4 * exp(-6 * t) * cos(2 * x);
}

double compute_error(const vec& last_layer, int nx)
{
    double max_error = 0.0;
    double h = X / nx;

    for (int i = 0; i < nx - 1; i++)
    {
        double current_error =
            std::abs(last_layer[i] - exact_solution(h * (i + 1), T));
        if (current_error > max_error)
        {
            max_error = current_error;
        }
    }

    return max_error;
}

