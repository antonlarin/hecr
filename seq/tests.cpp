#include "gtest/gtest.h"
#include "algorithms.hpp"

TEST(HecrSequential, computeLinearSystemSize)
{
    EXPECT_EQ(1, compute_linear_system_size(1));
    EXPECT_EQ(7, compute_linear_system_size(4));
    EXPECT_EQ(7, compute_linear_system_size(7));
    EXPECT_EQ(15, compute_linear_system_size(8));
}

