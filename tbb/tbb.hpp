#ifndef _TBB_H
#define _TBB_H

#include "algorithms.hpp"
#include "tbb/parallel_for.h"
#include "tbb/blocked_range.h"
#include "tbb/partitioner.h"
#include "tbb/task_scheduler_init.h"

void fcompnewcoef(double* new_a, double* new_b, double* new_c, double* new_rhs, double* old_a, double* old_b, double* old_c, double* old_rhs,
	int equation_count, int grainsize);

void fcompresult(double* as, double* bs, double* cs, double* rhs, vec* result, int factor, int equation_count, int grainsize);

#endif