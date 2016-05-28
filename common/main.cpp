#include <fstream>
#include <chrono>

#include "solver.hpp"

int main(int, char** argv)
{
    std::ios_base::sync_with_stdio(false);

    std::string input_file_name(argv[1]);
    std::string output_file_name(argv[2]);
    std::string timing_file_name(argv[3]);

    std::fstream input_file(input_file_name, std::ios::in);
    int nx, nt;
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
