#include "tbb.hpp"

class FFCompNewCoefFunctor
{
private:
	double* new_a;
	double* new_b;
	double* new_c;
	double* new_rhs;
	double* old_a;
	double* old_b; 
	double* old_c;
	double* old_rhs;
public:
	FFCompNewCoefFunctor(double* _new_a, double* _new_b, double* _new_c, double* _new_rhs,
		double* _old_a, double* _old_b, double* _old_c, double* _old_rhs) 
	{
		new_a = _new_a;
		new_b = _new_b;
		new_c = _new_c;
		new_rhs = _new_rhs;
		old_a = _old_a;
		old_b = _old_b;
		old_c = _old_c;
		old_rhs = _old_rhs;
	}

	void operator()(const tbb::blocked_range<int> &r) const
	{
		int i, k;
		for (i = r.begin(); i < r.end(); i++)
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
	}
};

void fcompnewcoef(double* _new_a, double* _new_b, double* _new_c, double* _new_rhs,
	double* _old_a, double* _old_b, double* _old_c, double* _old_rhs,
	int equation_count, int grainsize)
{
	tbb::parallel_for<tbb::blocked_range<int>, FFCompNewCoefFunctor>(
		tbb::blocked_range<int>(0, equation_count, grainsize),
		FFCompNewCoefFunctor(_new_a, _new_b, _new_c, _new_rhs, _old_a, _old_b, _old_c, _old_rhs));
}

class FFCompResultFunctor
{
private:
	double* as;
	double* bs;
	double* cs;
	double* rhs;
	vec* result;
	int factor;
public:
	FFCompResultFunctor(double* _as, double* _bs, double* _cs, double* _rhs, vec* _result, int _factor)
	{
		as = _as;
		bs = _bs;
		cs = _cs;
		rhs = _rhs;
		result = _result;
		factor = _factor;
	}

	void operator()(const tbb::blocked_range<int> &r) const
	{
		for (int j = r.begin(); j < r.end(); j++)
		{
			int k = j * 2 - 1;
			// index in coefficents arrays is less by 1 since they don't
			// include x_0=0 and x_N=0 equations
			int coeff_j = k - 1;
			// mapping j to the index in the full extended system
			int result_j = k * factor;
			int result_j_minus1 = result_j - factor;
			int result_j_plus1 = result_j + factor;

			(*result)[result_j] =
				(rhs[coeff_j] - bs[coeff_j] * (*result)[result_j_minus1] -
				cs[coeff_j] * (*result)[result_j_plus1]) / as[coeff_j];
		}
	}
};

void fcompresult(double* as, double* bs, double* cs, double* rhs, vec* result, int factor, int equation_count, int grainsize)
{
	tbb::parallel_for<tbb::blocked_range<int>, FFCompResultFunctor>(
		tbb::blocked_range<int>(1, equation_count, grainsize),
		FFCompResultFunctor(as, bs, cs, rhs, result, factor));
}