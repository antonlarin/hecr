#include <iostream>
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
    vec bs(size, 0.0);
    vec cs(size, 0.0);

    double a = - nt / T - 2 * nx * nx / X / X;
    double b_and_c = nx * nx / X / X;

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

vec cyclic_reduction(const TridiagonalMatrix& m, const vec& b)
{
    vec old_a = m.as;
    vec old_b = m.bs;
    vec old_c = m.cs;
    vec old_rhs = b;
    int equation_count = m.size;

    while (equation_count > 1)
    {
        equation_count /= 2;
        vec new_a(equation_count, 0.0);
        vec new_b(equation_count, 0.0);
        vec new_c(equation_count, 0.0);
        vec new_rhs(equation_count, 0.0);

        for (int i = 0; i < equation_count; i++)
        {
            int old_i = 2 * i + 1;
            new_a[i] = old_a[old_i] -
                old_b[old_i] * old_c[old_i - 1] / old_a[old_i - 1] -
                old_c[old_i] * old_b[old_i + 1] / old_a[old_i + 1];
            new_b[i] = -old_b[old_i - 1] * old_b[old_i] / old_a[old_i - 1];
            new_c[i] = -old_c[old_i] * old_c[old_i + 1] / old_a[old_i + 1];
            new_rhs[i] = old_rhs[old_i] -
                old_rhs[old_i - 1] * old_b[old_i] / old_a[old_i - 1] -
                old_rhs[old_i + 1] * old_c[old_i] / old_a[old_i + 1];
        }
        std::cout << "eq_count=" << equation_count << std::endl;
        std::cout << "as = [";
        for (int i = 0; i < equation_count; i++)
        {
            std::cout << new_a[i] << ", ";
        }
        std::cout << "\b\b]" << std::endl;
        std::cout << "bs = [";
        for (int i = 0; i < equation_count; i++)
        {
            std::cout << new_b[i] << ", ";
        }
        std::cout << "\b\b]" << std::endl;
        std::cout << "cs = [";
        for (int i = 0; i < equation_count; i++)
        {
            std::cout << new_c[i] << ", ";
        }
        std::cout << "\b\b]" << std::endl;
        std::cout << "rhs = [";
        for (int i = 0; i < equation_count; i++)
        {
            std::cout << new_rhs[i] << ", ";
        }
        std::cout << "\b\b]" << std::endl;

        old_a = new_a;
        old_b = new_b;
        old_c = new_c;
        old_rhs = new_rhs;
    }

    // result is padded with one zero one both sides to emulate
    // implicit variables x_0 and x_N
    vec result(m.size + 2, 0.0);
    int fuv = m.size / 2; // fuv = first updated var
    int stride = (fuv + 1) * 2;
    result[fuv + 1] = old_rhs[0] / old_a[0];
    std::cout << "result[fuv]=" << result[fuv + 1] << std::endl;

    while (equation_count < m.size)
    {
        fuv /= 2;
        stride /= 2;
        int halfstride = stride / 2;

        for (unsigned int i = fuv; i < result.size(); i += stride)
        {
            result[i + 1] =
                (b[i] - m.bs[i] * result[i - halfstride + 1] -
                 m.cs[i] * result[i + halfstride + 1]) / m.as[i];
        }

        equation_count = equation_count * 2 + 1;
        std::cout << "eq_count=" << equation_count << std::endl;
        std::cout << "result = [";
        for (unsigned int i = 0; i < result.size(); i++)
        {
            std::cout << result[i] << ", ";
        }
        std::cout << "\b\b]" << std::endl;
    }

    // drop the padding
    result.pop_back();
    result.erase(result.begin());

    return result;
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

double solve_problem(int nx, int nt)
{
    int size = compute_linear_system_size(nx);
    std::cout << "N=" << size << std::endl;
    vec current_layer(size, 0.0);
    double h = X / nx;
    double tau = T / nt;

    std::cout << "u0 = [";
    for (int i = 0; i < nx - 1; i++)
    {
        current_layer[i] = exact_solution(h * (i + 1), 0.0);
        std::cout << current_layer[i] << ", ";
    }
    std::cout << "\b\b]" << std::endl;

    TridiagonalMatrix m = construct_matrix(nx, nt);
    std::cout << "as = [";
    for (int i = 0; i < m.size; i++)
    {
        std::cout << m.as[i] << ", ";
    }
    std::cout << "\b\b]" << std::endl;
    std::cout << "bs = [";
    for (int i = 0; i < m.size; i++)
    {
        std::cout << m.bs[i] << ", ";
    }
    std::cout << "\b\b]" << std::endl;
    std::cout << "cs = [";
    for (int i = 0; i < m.size; i++)
    {
        std::cout << m.cs[i] << ", ";
    }
    std::cout << "\b\b]" << std::endl;
    for (int i = 1; i <= nt; i++)
    {
        vec b = construct_rhs(nx, nt, current_layer, i * tau,
                f, exact_solution);
        std::cout << "rhs = [";
        for (unsigned int i = 0; i < b.size(); i++)
        {
            std::cout << b[i] << ", ";
        }
        std::cout << "\b\b]" << std::endl;
        current_layer = cyclic_reduction(m, b);
        std::cout << "last = [";
        for (unsigned int i = 0; i < current_layer.size(); i++)
        {
            std::cout << current_layer[i] << ", ";
        }
        std::cout << "\b\b]" << std::endl;
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
