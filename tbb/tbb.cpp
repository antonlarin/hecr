#include "tbb.hpp"

class FFCompNewCoefFunctor
{
private:
	vec* new_a;
	vec* new_b;
	vec* new_c;
	vec* new_rhs;
	vec* old_a;
	vec* old_b; 
	vec* old_c;
	vec* old_rhs;
public:
	FFCompNewCoefFunctor(vec* _new_a, vec* _new_b, vec* _new_c, vec* _new_rhs,
		vec* _old_a, vec* _old_b, vec* _old_c, vec* _old_rhs) 
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
			(*new_a)[i] = (*old_a)[old_i] -
				(*old_b)[old_i] * (*old_c)[old_i - 1] / (*old_a)[old_i - 1] -
				(*old_c)[old_i] * (*old_b)[old_i + 1] / (*old_a)[old_i + 1];
			(*new_b)[i] = -(*old_b)[old_i - 1] * (*old_b)[old_i] / (*old_a)[old_i - 1];
			(*new_c)[i] = -(*old_c)[old_i] * (*old_c)[old_i + 1] / (*old_a)[old_i + 1];
			(*new_rhs)[i] = (*old_rhs)[old_i] -
				(*old_rhs)[old_i - 1] * (*old_b)[old_i] / (*old_a)[old_i - 1] -
				(*old_rhs)[old_i + 1] * (*old_c)[old_i] / (*old_a)[old_i + 1];
		}
	}
};

void fcompnewcoef(vec* _new_a, vec* _new_b, vec* _new_c, vec* _new_rhs,
	vec* _old_a, vec* _old_b, vec* _old_c, vec* _old_rhs,
	int equation_count, int grainsize)
{
	tbb::parallel_for<tbb::blocked_range<int>, FFCompNewCoefFunctor>(
		tbb::blocked_range<int>(0, equation_count, grainsize),
		FFCompNewCoefFunctor(_new_a, _new_b, _new_c, _new_rhs, _old_a, _old_b, _old_c, _old_rhs));
}

class FFCompResultFunctor
{
private:
	vec* as;
	vec* bs;
	vec* cs;
	vec* rhs;
	vec* result;
	int factor;
public:
	FFCompResultFunctor(vec* _as, vec* _bs, vec* _cs, vec* _rhs, vec* _result, int _factor)
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
				((*rhs)[coeff_j] - (*bs)[coeff_j] * (*result)[result_j_minus1] -
				(*cs)[coeff_j] * (*result)[result_j_plus1]) / (*as)[coeff_j];
		}
	}
};

void fcompresult(vec* as, vec* bs, vec* cs, vec* rhs, vec* result, int factor, int equation_count, int grainsize)
{
	tbb::parallel_for<tbb::blocked_range<int>, FFCompResultFunctor>(
		tbb::blocked_range<int>(1, equation_count, grainsize),
		FFCompResultFunctor(as, bs, cs, rhs, result, factor));
}