#include <fstream>
#include <chrono>
#include <cmath>
#include <vector>

typedef std::vector<double> vec;

class TridiagonalMatrix
{
};

TridiagonalMatrix construct_matrix()
{
    return TridiagonalMatrix();
}

vec construct_rhs()
{
    return vec();
}

vec cyclic_reduction(const TridiagonalMatrix& /*a*/, const vec& /*b*/)
{
    return vec();
}

double exact_solution(double x, double t)
{
    return sin(x - t);
}

double compute_error(const vec& /*last_layer*/)
{
    return 1e-4;
}

double solve_problem(int nx, int nt)
{
    vec initial_conditions;
    double h = 1.0 / nx;

    for (int i = 0; i <= nx; i++)
    {
        initial_conditions.push_back(exact_solution(h * i, 0.0));
    }

    vec current_layer;
    for (int i = 1; i <= nt; i++)
    {
        current_layer.clear();
        TridiagonalMatrix a = construct_matrix();
        vec b = construct_rhs();
        current_layer = cyclic_reduction(a, b);
    }

    double error = compute_error(current_layer);

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
