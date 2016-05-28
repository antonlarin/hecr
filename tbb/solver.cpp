#include "solver.hpp"

vec cyclic_reduction(const TridiagonalMatrix& m, const vec& b)
{
	std::vector<vec > a_vecs;
	a_vecs.push_back(m.as);
	std::vector<vec > b_vecs;
	b_vecs.push_back(m.bs);
	std::vector<vec > c_vecs;
	c_vecs.push_back(m.cs);
	std::vector<vec > rhs_vecs;
	rhs_vecs.push_back(b);
	int equation_count = m.size;
	int grainsizeForCoef = 1;
	int grainsizeForRes = 1;

	while (equation_count > 1)
	{
		equation_count /= 2;
		vec new_a(equation_count, 0.0);
		vec new_b(equation_count, 0.0);
		vec new_c(equation_count, 0.0);
		vec new_rhs(equation_count, 0.0);
		vec& old_a = a_vecs.back();
		vec& old_b = b_vecs.back();
		vec& old_c = c_vecs.back();
		vec& old_rhs = rhs_vecs.back();

		//if (equation_count > 128)
		//	grainsize = 32;
		//else
		//	grainsize = 1;

		grainsizeForCoef = equation_count;
		/*if (equation_count > 256)
			grainsize = 64;*/
		if (equation_count > 1024)
			grainsizeForCoef = 256;

		fcompnewcoef(&new_a, &new_b, &new_c, &new_rhs, &old_a, &old_b, &old_c, &old_rhs,
			equation_count, grainsizeForCoef);

		a_vecs.push_back(new_a);
		b_vecs.push_back(new_b);
		c_vecs.push_back(new_c);
		rhs_vecs.push_back(new_rhs);
	}

	// result is padded with one zero one both sides to emulate
	// implicit variables x_0 and x_N
	vec result(m.size + 2, 0.0);
	// factor maps indices of equations in full system extended with x_0=0
	// and x_N=0 onto indices of equations in reduced extended system
	int factor = m.size / 2 + 1;
	int n = 1;
	for (unsigned int i = a_vecs.size(); i-- > 0;)
	{
		vec& as = a_vecs[i];
		vec& bs = b_vecs[i];
		vec& cs = c_vecs[i];
		vec& rhs = rhs_vecs[i];
		grainsizeForRes = equation_count + 1;
		if (n > 512)
			grainsizeForRes = 128;
		// j is an index of an equation in the reduced system extended with
		// x_0=0 and x_N=0
		
		fcompresult(&as, &bs, &cs, &rhs, &result, factor, n + 1, grainsizeForRes );

		n = equation_count + 1;
		equation_count = equation_count * 2 + 1;
		factor /= 2;
	}

	// drop the padding
	result.pop_back();
	result.erase(result.begin());

	return result;
}

double solve_problem(int nx, int nt)
{
	int size = compute_linear_system_size(nx - 1);
	vec current_layer(size, 0.0);
	double h = X / nx;
	double tau = T / nt;

	for (int i = 0; i < nx - 1; i++)
	{
		current_layer[i] = exact_solution(h * (i + 1), 0.0);
	}

	TridiagonalMatrix m = construct_matrix(nx, nt);

	tbb::task_scheduler_init init(2);

	std::chrono::duration<double> elapsed_time, elapsed_timeOnce;
	std::fstream timing_file("timeReduct.txt", std::ios::out);

	for (int i = 1; i <= nt; i++)
	{
		vec b = construct_rhs(nx, nt, current_layer, i * tau,
			f, exact_solution);

		auto start_time = std::chrono::high_resolution_clock::now();
		current_layer = cyclic_reduction(m, b);
		auto end_time = std::chrono::high_resolution_clock::now();
		elapsed_time += end_time - start_time;
	}

	timing_file << elapsed_time.count() << std::endl;
	init.terminate();
	double error = compute_error(current_layer, nx);

	return error;
}

vec construct_rhs(int nx, int nt, const vec& previous_layer, double t,
	function_2var f, function_2var u)
{
	static int size = compute_linear_system_size(nx - 1);

	vec rhs(size, 0.0);

	static double h = X / nx;
	static double tau = T / nt;

        fcomprhs(&rhs, &previous_layer, nx, h, tau, t, f);

	rhs[0] -= u(0.0, t) / h / h;
	rhs[nx - 2] -= u(X, t) / h / h;

	return rhs;
}

