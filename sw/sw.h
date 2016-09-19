/*
Copyright (C) 2016 Dmitry Sigaev

Contact: Dmitry Sigaev <dima.sigaev@gmail.com>
*/


#ifndef _SW_H_
#define _SW_H_

#include "lal_typedefs.h"

double sw_constant_gap_double(const search_swcg_profile_t * sp, const sequence_t * dseq, const sequence_t * qseq);
int64_t sw_constant_gap_int(const search_swcg_profile_int_t * sp, const sequence_t * dseq, const sequence_t * qseq);

#endif /* _SW_H_ */