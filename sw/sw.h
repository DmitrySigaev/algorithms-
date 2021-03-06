/*
Copyright (C) 2016 Dmitry Sigaev

Contact: Dmitry Sigaev <dima.sigaev@gmail.com>
*/


#ifndef _SW_H_
#define _SW_H_

#include "lal_typedefs.h"

double sw_constant_gap_double(const search_swcg_profile_t * sp, const sequence_t * dseq, const sequence_t * qseq);
int64_t sw_constant_gap_int(const search_swcg_profile_int_t * sp, const sequence_t * dseq, const sequence_t * qseq);

double sw_affine_gap(const search_swag_profile_t * sp, const sequence_t * dseq, const sequence_t * qseq);
region_t sw_alignment_swipe(const search_swag_profile_t * sp, const sequence_t *xseq, const sequence_t *yseq);
score_matrix_t sw_directions(const search_swag_profile_t * sp, const sequence_t *xseq, const sequence_t *yseq);
double sw_thr(search_thr_profile_t * sp, const sequence_t *dseq, const sequence_t *qseq);
search_thr_profile_t * search_thr_init(search_swag_profile_t *s, size_t thr);
void search_thr_deinit(search_thr_profile_t * sthr, size_t thr);
#endif /* _SW_H_ */