#include <fstream>
#include <chrono>
#include <cmath>
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
private:
    vec as;
    vec bs;
    vec cs;
    int size;
};

// compute the size of the system suitable for the cyclic reduction method
int compute_linear_system_size(int nx)
{
    int size = 2;
    while (size - 1 <= nx - 1)
    {
        size = size * 2;
    }

    return size - 1;
}


TridiagonalMatrix construct_matrix(int nx, int nt)
{
    int size = compute_linear_system_size(nx);

    vec as(size, 1.0);
    vec bs_cs(size - 1, 0.0);

    double a = - nt / T - 2 * nx * nx / X / X;
    double b_and_c = nx * nx / X / X;

    for (int i = 0; i < nx - 2; i++)
    {
        as[i] = a;
        bs_cs[i] = b_and_c;
    }
    as[nx - 1] = a;

    return TridiagonalMatrix(as, bs_cs, bs_cs);
}

vec construct_rhs(int nx, int nt, const vec& previous_layer, double t,
        function_2var f, function_2var u)
{
    static int size = compute_linear_system_size(nx);

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

vec cyclic_reduction(const TridiagonalMatrix& /*a*/, const vec& /*b*/)
{
    // TODO replace with an algorithm
    return vec(1000, 0.0);
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

    for (int i = 0; i <= nx; i++)
    {
        double current_error =
            std::abs(last_layer[i] - exact_solution(h * i, 1.0));
        if (current_error > max_error)
        {
            max_error = current_error;
        }
    }

    return max_error;
}

double solve_problem(int nx, int nt)
{
    vec initial_conditions;
    double h = X / nx;
    double tau = T / nt;

    for (int i = 0; i < nx - 1; i++)
    {
        initial_conditions.push_back(exact_solution(h * (i + 1), 0.0));
    }

    TridiagonalMatrix a = construct_matrix(nx, nt);
    int size = compute_linear_system_size(nx);
    vec current_layer(size, 0.0);
    for (int i = 1; i <= nt; i++)
    {
        vec b = construct_rhs(nx, nt, current_layer, i * tau,
                f, exact_solution);
        current_layer = cyclic_reduction(a, b);
    }

    double error = compute_error(current_layer, nx);

    return error;
}

int main(int, char** argv)
{
    std::ios_base::sync_with_stdio(false);

    std::string input_file_name(argv[1]);
    std::string output_file_name(argv[2]);
    std::string timing_file_name(argv[3]);

    std::fstream input_file(input_file_name, std::ios::in);
    double nx, nt;
    input_file >> nx >> nt;

    auto start_time = std::chrono::high_resolution_clock::now();
    double error = solve_problem(nx, nt);
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_time = end_time - start_time;

    std::fstream output_file(output_file_name, std::ios::out);
    output_file << error << std::endl;

    std::fstream timing_file(timing_file_name, std::ios::out);
    timing_file << elapsed_time.count() << std::endl;

    return 0;
}
