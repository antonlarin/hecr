//#include  "tbb.hpp"
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
//	int grainsize = 1;
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
//		//if (equation_count > 128)
//		//	grainsize = 32;
//		//else
//		//	grainsize = 1;
//
//		grainsize = equation_count;
//		/*if (equation_count > 256)
//			grainsize = 64;*/
//		if (equation_count > 1024)
//			grainsize = 256;
//
//		fcompnewcoef(&new_a, &new_b, &new_c, &new_rhs, old_a, old_b,old_c, old_rhs,
//			equation_count, grainsize);
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