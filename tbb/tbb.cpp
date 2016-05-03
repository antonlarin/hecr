#include "tbb.hpp"

class FFCompNewCoefFunctor
{
private:
	mutable vec* new_a;
	mutable vec* new_b;
	mutable vec* new_c;
	mutable vec* new_rhs;
	vec old_a;
	vec old_b; 
	vec old_c;
	vec old_rhs;
public:
	FFCompNewCoefFunctor(vec* _new_a, vec* _new_b, vec* _new_c, vec* _new_rhs,
		vec& _old_a, vec& _old_b, vec& _old_c, vec& _old_rhs)
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
			(*new_a)[i] = old_a[old_i] -
				old_b[old_i] * old_c[old_i - 1] / old_a[old_i - 1] -
				old_c[old_i] * old_b[old_i + 1] / old_a[old_i + 1];
			(*new_b)[i] = -old_b[old_i - 1] * old_b[old_i] / old_a[old_i - 1];
			(*new_c)[i] = -old_c[old_i] * old_c[old_i + 1] / old_a[old_i + 1];
			(*new_rhs)[i] = old_rhs[old_i] -
				old_rhs[old_i - 1] * old_b[old_i] / old_a[old_i - 1] -
				old_rhs[old_i + 1] * old_c[old_i] / old_a[old_i + 1];
		}
	}
};

void fcompnewcoef(vec* _new_a, vec* _new_b, vec* _new_c, vec* _new_rhs,
	vec& _old_a, vec& _old_b, vec& _old_c, vec& _old_rhs,
	int equation_count, int grainsize)
{
	tbb::parallel_for<tbb::blocked_range<int>, FFCompNewCoefFunctor>(
		tbb::blocked_range<int>(0, equation_count, grainsize),
		FFCompNewCoefFunctor(_new_a, _new_b, _new_c, _new_rhs, _old_a, _old_b, _old_c, _old_rhs));
}