#include "solver.hpp"

vec cyclic_reduction(const TridiagonalMatrix& m, const vec& b)
{
	int equation_count = m.size;
	int reductions_count = static_cast<int>(floor(log(equation_count) /
		log(2)));
	int coeff_array_size = 2 * (equation_count + 1) - reductions_count - 3;
	double* a_vecs = new double[coeff_array_size];
	double* b_vecs = new double[coeff_array_size];
	double* c_vecs = new double[coeff_array_size];
	double* rhs_vecs = new double[coeff_array_size];
	memcpy(a_vecs, m.as.data(), equation_count * sizeof(double));
	memcpy(b_vecs, m.bs.data(), equation_count * sizeof(double));
	memcpy(c_vecs, m.cs.data(), equation_count * sizeof(double));
	memcpy(rhs_vecs, b.data(), equation_count * sizeof(double));
	int base_idx = equation_count;
	int grainsizeForCoef = 1;
	int grainsizeForRes = 1;

	while (equation_count > 1)
	{
		int old_base_idx = base_idx - equation_count;
		equation_count /= 2;

		double* old_a = a_vecs + old_base_idx;
		double* old_b = b_vecs + old_base_idx;
		double* old_c = c_vecs + old_base_idx;
		double* old_rhs = rhs_vecs + old_base_idx;
		double* new_a = a_vecs + base_idx;
		double* new_b = b_vecs + base_idx;
		double* new_c = c_vecs + base_idx;
		double* new_rhs = rhs_vecs + base_idx;

		grainsizeForCoef = equation_count;
		/*if (equation_count > 256)
		grainsize = 64;*/
		if (equation_count > 128)
			grainsizeForCoef = 32;

		fcompnewcoef(new_a, new_b, new_c, new_rhs, old_a, old_b, old_c, old_rhs,
			equation_count, grainsizeForCoef);

		base_idx += equation_count;
	}

	// result is padded with one zero one both sides to emulate
	// implicit variables x_0 and x_N
	vec result(m.size + 2, 0.0);
	// factor maps indices of equations in full system extended with x_0=0
	// and x_N=0 onto indices of equations in reduced extended system
	int factor = m.size / 2 + 1;
	int n = 1;
	for (unsigned int i = reductions_count + 1; i-- > 0;)
	{
		base_idx -= equation_count;
		double* as = a_vecs + base_idx;
		double* bs = b_vecs + base_idx;
		double* cs = c_vecs + base_idx;
		double* rhs = rhs_vecs + base_idx;

		// j is an index of an equation in the reduced system extended with
		// x_0=0 and x_N=0

		grainsizeForRes = equation_count + 1;
		if (n > 256)
			grainsizeForRes = 64;

		fcompresult(as, bs, cs, rhs, &result, factor, n + 1, grainsizeForRes);

		n = equation_count + 1;
		equation_count = equation_count * 2 + 1;
		factor /= 2;
	}

	// drop the padding
	result.pop_back();
	result.erase(result.begin());

	delete[] a_vecs;
	delete[] b_vecs;
	delete[] c_vecs;
	delete[] rhs_vecs;

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

	tbb::task_scheduler_init init(4);

	for (int i = 1; i <= nt; i++)
	{
		vec b = construct_rhs(nx, nt, current_layer, i * tau,
			f, exact_solution);

		current_layer = cyclic_reduction(m, b);
	}

	double error = compute_error(current_layer, nx);
	init.terminate();

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

