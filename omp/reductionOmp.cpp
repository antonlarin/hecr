//#include "algorithms.hpp"
//#include <omp.h>
//#define OMP_SCHEDULE dynamic
//#define CHUNK 32
//
//vec cyclic_reduction(const TridiagonalMatrix& m, const vec& b)
//{
//	std::vector<vec > a_vecs;
//	a_vecs.push_back(m.as);
//	std::vector<vec > b_vecs;
//	b_vecs.push_back(m.bs);
//	std::vector<vec > c_vecs;
//	c_vecs.push_back(m.cs);
//	std::vector<vec > rhs_vecs;
//	rhs_vecs.push_back(b);
//	int equation_count = m.size;
//
//	while (equation_count > 1)
//	{
//		equation_count /= 2;
//		vec new_a(equation_count, 0.0);
//		vec new_b(equation_count, 0.0);
//		vec new_c(equation_count, 0.0);
//		vec new_rhs(equation_count, 0.0);
//		vec& old_a = a_vecs.back();
//		vec& old_b = b_vecs.back();
//		vec& old_c = c_vecs.back();
//		vec& old_rhs = rhs_vecs.back();
//
//#pragma omp parallel for if (equation_count > 128) num_threads(4)
//		for (int i = 0; i < equation_count; i++)
//		{
//			int old_i = 2 * i + 1;
//			double oldA_next = old_a[old_i + 1];
//			double oldA_prev = old_a[old_i - 1];
//
//			double CA = old_c[old_i] / oldA_next;
//			double BA = old_b[old_i] / oldA_prev;
//
//			new_a[i] = old_a[old_i] -
//				BA * old_c[old_i - 1] -
//				CA * old_b[old_i + 1];
//			new_b[i] = -old_b[old_i - 1] * BA;
//			new_c[i] = -CA * old_c[old_i + 1];
//			new_rhs[i] = old_rhs[old_i] -
//				old_rhs[old_i - 1] * BA -
//				old_rhs[old_i + 1] * CA;
//		}
//
//		a_vecs.push_back(new_a);
//		b_vecs.push_back(new_b);
//		c_vecs.push_back(new_c);
//		rhs_vecs.push_back(new_rhs);
//	}
//
//	// result is padded with one zero one both sides to emulate
//	// implicit variables x_0 and x_N
//	vec result(m.size + 2, 0.0);
//	// factor maps indices of equations in full system extended with x_0=0
//	// and x_N=0 onto indices of equations in reduced extended system
//	int factor = m.size / 2 + 1;
//	for (unsigned int i = a_vecs.size(); i-- > 0;)
//	{
//		vec& as = a_vecs[i];
//		vec& bs = b_vecs[i];
//		vec& cs = c_vecs[i];
//		vec& rhs = rhs_vecs[i];
//
//		// j is an index of an equation in the reduced system extended with
//		// x_0=0 and x_N=0
//#pragma omp parallel for if (equation_count > 128) num_threads(4)
//		for (int j = 1; j < equation_count + 1; j += 2)
//		{
//			// index in coefficents arrays is less by 1 since they don't
//			// include x_0=0 and x_N=0 equations
//			int coeff_j = j - 1;
//			// mapping j to the index in the full extended system
//			int result_j = j * factor;
//			int result_j_minus1 = result_j - factor;
//			int result_j_plus1 = result_j + factor;
//
//			result[result_j] =
//				(rhs[coeff_j] - bs[coeff_j] * result[result_j_minus1] -
//				cs[coeff_j] * result[result_j_plus1]) / as[coeff_j];
//		}
//
//		equation_count = equation_count * 2 + 1;
//		factor /= 2;
//	}
//
//	// drop the padding
//	result.pop_back();
//	result.erase(result.begin());
//
//	return result;
//}