/*
Copyright (C) 2016 Dmitry Sigaev

Contact: Dmitry Sigaev <dima.sigaev@gmail.com>
*/


#ifndef _FP_H_
#define _FP_H_

#include "lal_typedefs.h"

double fp_thr(search_fp_thr_profile_t * sp, const sequence_t *dseq, const sequence_t *qseq);
search_fp_thr_profile_t * search_fp_thr_init(search_fp_profile_t *s, size_t thr);
void search_fp_thr_deinit(search_fp_thr_profile_t * sthr, size_t thr);
#endif /* _FP_H_ */