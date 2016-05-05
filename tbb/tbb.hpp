#ifndef _TBB_H
#define _TBB_H

#include "algorithms.hpp"
#include "tbb/parallel_for.h"
#include "tbb/blocked_range.h"
#include "tbb/partitioner.h"
#include "tbb/task_scheduler_init.h"

void fcompnewcoef(vec* new_a, vec* new_b, vec* new_c, vec* new_rhs, vec* old_a, vec* old_b, vec* old_c, vec* old_rhs,
	int equation_count, int grainsize);

void fcompresult(vec* as, vec* bs, vec* cs, vec* rhs, vec* result, int factor, int equation_count, int grainsize);

#endif