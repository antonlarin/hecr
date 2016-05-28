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

		for (int i = 0; i < equation_count; i++)
		{
			int old_i = 2 * i + 1;
			double oldA_next = old_a[old_i + 1];
			double oldA_prev = old_a[old_i - 1];

			double CA = old_c[old_i] / oldA_next;
			double BA = old_b[old_i] / oldA_prev;

			new_a[i] = old_a[old_i] -
				BA * old_c[old_i - 1] -
				CA * old_b[old_i + 1];
			new_b[i] = -old_b[old_i - 1] * BA;
			new_c[i] = -CA * old_c[old_i + 1];
			new_rhs[i] = old_rhs[old_i] -
				old_rhs[old_i - 1] * BA -
				old_rhs[old_i + 1] * CA;
		}

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
		for (int j = 1; j < n + 1; j++)
		{
			int k = j * 2 - 1;
			// index in coefficents arrays is less by 1 since they don't
			// include x_0=0 and x_N=0 equations
			int coeff_j = k - 1;
			// mapping j to the index in the full extended system
			int result_j = k * factor;
			int result_j_minus1 = result_j - factor;
			int result_j_plus1 = result_j + factor;

			result[result_j] =
				(rhs[coeff_j] - bs[coeff_j] * result[result_j_minus1] -
				cs[coeff_j] * result[result_j_plus1]) / as[coeff_j];
		}

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

	for (int i = 1; i <= nt; i++)
	{
		vec b = construct_rhs(nx, nt, current_layer, i * tau,
			f, exact_solution);

		current_layer = cyclic_reduction(m, b);
	}

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

	for (int i = 0; i < nx - 1; i++)
	{
		rhs[i] = -previous_layer[i] / tau - f(h * (i + 1), t);
	}

	rhs[0] -= u(0.0, t) / h / h;
	rhs[nx - 2] -= u(X, t) / h / h;

	return rhs;
}