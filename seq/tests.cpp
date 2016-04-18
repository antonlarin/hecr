#include <cmath>

#include "gtest/gtest.h"

#include "solver.hpp"

const double EPS = 1e-15;

double max_diff(const vec& a, const vec& b)
{
    double max_diff = 0.0;
    for (unsigned int i = 0; i < a.size(); i++)
    {
        double diff = std::abs(a[i] - b[i]);
        if (diff > max_diff)
        {
            max_diff = diff;
        }
    }

    return max_diff;
}

vec make_vector(int N, const double* data)
{
    return vec(data, data + N);
}


TEST(HecrSequential, computeLinearSystemSize)
{
    EXPECT_EQ(1, compute_linear_system_size(1));
    EXPECT_EQ(7, compute_linear_system_size(4));
    EXPECT_EQ(7, compute_linear_system_size(7));
    EXPECT_EQ(15, compute_linear_system_size(8));
}

TEST(HecrSequential, constructMatrix)
{
    int nx = 10;
    int nt = 10;
    int size = 15;
    double a = -334.22778765548086862041428227;
    const double as_values[] = { a, a, a, a, a, a, a, a, a, 1, 1, 1, 1, 1, 1 };
    vec as = make_vector(size, as_values);
    double bc = 162.1138938277404343102071411;
    const double bs_values[] =
        { 0, bc, bc, bc, bc, bc, bc, bc, bc, 0, 0, 0, 0, 0, 0 };
    vec bs = make_vector(size, bs_values);
    const double cs_values[] =
        { bc, bc, bc, bc, bc, bc, bc, bc, 0, 0, 0, 0, 0, 0, 0 };
    vec cs = make_vector(size, cs_values);

    TridiagonalMatrix m = construct_matrix(nx, nt);

    EXPECT_EQ(size, m.size);
    EXPECT_LT(max_diff(m.as, as), EPS);
    EXPECT_LT(max_diff(m.bs, bs), EPS);
    EXPECT_LT(max_diff(m.cs, cs), EPS);
}

double f_test(double x, double t)
{
    return 1.0;
}

double u_test(double x, double t)
{
    return 4.0;
}

TEST(HecrSequential, constructRhs)
{
    int nx = 10;
    int nt = 10;
    int size = 15;
    vec prev_layer(size, -2.0);
    double t = 0.5;
    double end = -629.4555753109617372408285645;
    double mid = 19.0;
    const double expected_rhs_values[] =
        { end, mid, mid, mid, mid, mid, mid, mid, end, 0, 0, 0, 0, 0, 0 };
    vec expected_rhs = make_vector(size, expected_rhs_values);

    vec rhs = construct_rhs(nx, nt, prev_layer, t, f_test, u_test);

    EXPECT_EQ(size, rhs.size());
    EXPECT_LT(max_diff(expected_rhs, rhs), EPS);
}

TEST(HecrSequential, cyclicReductionIdentityMatrix)
{
	Solver solver;
    int size = 15;
    vec diag(size, 1.0);
    vec subdiag(size, 0.0);
    TridiagonalMatrix m(diag, subdiag, subdiag);
    const double rhs_values[] = { 1.0, 2.0, -0.5, 0.4, 33.3, -11.0, 2.3, 4.4,
        -15.7, -5.8, -6.58, 0.11, 4.04, 5.02, 8.0 };
    vec rhs = make_vector(size, rhs_values);

	vec solution = solver.cyclic_reduction(m, rhs);

    EXPECT_LT(max_diff(rhs, solution), EPS);
}

TEST(HecrSequential, cyclicRedictionFullMatrix)
{
	Solver solver;
    int size = 7;
    const double as_values[] = { 3, 4, 3, 7, 5.5, -6, 4 };
    vec as = make_vector(size, as_values);
    const double bs_values[] = { 0, -1, 2, 1.1, 3.2, 2, 1 };
    vec bs = make_vector(size, bs_values);
    const double cs_values[] = { 2, 0.5, -0.7, 1, 0.8, 1, 0 };
    vec cs = make_vector(size, cs_values);
    TridiagonalMatrix m(as, bs, cs);
    const double rhs_values[] = { 7, 6.5, 1, 0.9, 8.6, 26, 13 };
    vec rhs = make_vector(size, rhs_values);
    const double xs_values[] = { 1, 2, -1, 0, 2, -3, 4 };
    vec xs = make_vector(size, xs_values);

	vec solution = solver.cyclic_reduction(m, rhs);

    EXPECT_LT(max_diff(xs, solution), EPS);
}

TEST(HecrSequential, cyclicReductionExtendedMatrix)
{
	Solver solver;
    int size = 15;
    const double as_values[] = { 3, 4, 3, 7, 5.5, -6, 4, -1, 3,
        1, 1, 1, 1, 1, 1 };
    vec as = make_vector(size, as_values);
    const double bs_values[] = { 0, -1, 2, 1.1, 3.2, 2, 1, -0.3, -2,
        0, 0, 0, 0, 0, 0 };
    vec bs = make_vector(size, bs_values);
    const double cs_values[] = { 2, 0.5, -0.7, 1, 0.8, 1, 0.2, 0.5, 0,
        0, 0, 0, 0, 0, 0 };
    vec cs = make_vector(size, cs_values);
    TridiagonalMatrix m(as, bs, cs);
    const double rhs_values[] = { 7, 6.5, 1, 0.9, 8.6, 26, 12.6, 0.3, 1,
        0, 0, 0, 0, 0, 0 };
    vec rhs = make_vector(size, rhs_values);
    const double xs_values[] = { 1, 2, -1, 0, 2, -3, 4, -2, -1,
        0, 0, 0, 0, 0, 0 };
    vec xs = make_vector(size, xs_values);

    vec solution = solver.cyclic_reduction(m, rhs);

    EXPECT_LT(max_diff(xs, solution), EPS);
}

TEST(HecrSequential, cyclicReductionNumericalSchemeLikeSystem)
{
	Solver solver;
	int size = 15;
	const double as_values[] = { 3, 3, 3, 3, 3, 3, 3, 3, 3,
		3, 3, 1, 1, 1, 1 };
	vec as = make_vector(size, as_values);
	const double bs_values[] = { 0, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5,
		1.5, 1.5, 0, 0, 0, 0 };
	vec bs = make_vector(size, bs_values);
	const double cs_values[] = { 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5,
		1.5, 0, 0, 0, 0, 0 };
	vec cs = make_vector(size, cs_values);
	TridiagonalMatrix m(as, bs, cs);
	const double rhs_values[] = { 21, 13.5, -10.5, 9, 52.5, 51, 22.5, 9, -27,
		-43.5, -1.5, 0, 0, 0, 0 };
	vec rhs = make_vector(size, rhs_values);
	const double xs_values[] = { 4, 6, -7, 1, 11, 12, -1, 5, -3,
		-17, 8, 0, 0, 0, 0 };
	vec xs = make_vector(size, xs_values);

	vec solution = solver.cyclic_reduction(m, rhs);

	EXPECT_LT(max_diff(xs, solution), EPS);
}

TEST(HecrSequential, cyclicReductionLargeSystem)
{
	Solver solver;
	int size = 65535;
        int actual_size = 47000;
        double a = 3.21;
        double bc = 1.4;
        double x = -0.7;
        double rhs_middle = -4.207;
        double rhs_end = -3.227;
	vec as(size, 1.0);
	vec bs(size, 0.0);
	vec cs(size, 0.0);

        as[0] = a;
        cs[0] = bc;
        for (int i = 1; i < actual_size - 1; i++)
        {
            as[i] = a;
            bs[i] = bc;
            cs[i] = bc;
        }
        as[actual_size - 1] = a;
        bs[actual_size - 1] = bc;

	TridiagonalMatrix m(as, bs, cs);

	vec rhs(size, 0.0);
        rhs[0] = rhs_end;
        for (int i = 1; i < actual_size - 1; i++)
        {
            rhs[i] = rhs_middle;
        }
        rhs[actual_size - 1] = rhs_end;

	vec xs(size, 0.0);
        for (int i = 0; i < actual_size; i++)
        {
            xs[i] = x;
        }

	vec solution = solver.cyclic_reduction(m, rhs);

	EXPECT_LT(max_diff(xs, solution), EPS);
}

