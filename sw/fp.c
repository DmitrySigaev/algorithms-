/*
Copyright (C) 2016

Contact: Dmitry Sigaev <dima.sigaev@gmail.com>
*/

#include <stdio.h>
#include <malloc.h>
#include <string.h>
#include "lal_matrix.h"
#include "lal_report.h"
#include "sw.h"

#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define SCORE(a, b, c, d) (((a) == (b)) ? (c) : (d))
#define VDTABLE(a, b) (sp->mtx->sc_double_matrix.ddata[(a)][(b)])
#define VITABLE(a, b) (sp->mtx->sc_int_matrix.idata[(a)][(b)])
#if 0
static double fp_ms(search_fp_thr_profile_t *s, const sequence_t * dseq, const sequence_t * qseq)
{
	const search_fp_profile_t *sp = s->sp;
	double *ph = s->h;
	double *pe = s->e;
	double v;
	double F, global_max = 0, h11, h01;
	double M10, M11, M12, M13, M14, M15;
	double M00, M01, M02, M03, M04, M05;

	double M01 = M02 = M03 = M04 = M05 = 0;
	M14 = M4[j];
	M13 = M3[j];
	M12 = M2[j];
	M11 = M1[j];

	for (size_t i = 0; i < 3; i++)
		for (size_t j = 0; j < qseq->len + 1; j++)
			mattr[i][j][s_fgap] = mattr[i][j][s_xgap] = mattr[i][j][s_ygap] = mattr[i][j][s_match] = 0;

	sse2_t M15 = sse2_shift_left(M[sse2_len - 1], 1);
	sse2_t M14 = sse2_shift_left(M4[sse2_len - 1], 1);
	sse2_t M13 = sse2_shift_left(M3[sse2_len - 1], 1);
	sse2_t M12 = sse2_shift_left(M2[sse2_len - 1], 1);
	sse2_t M11 = sse2_shift_left(M1[sse2_len - 1], 1);



	for (size_t i = 0; i < dseq->len; i++) {
		F = 0;
		h11 = 0;
		// [i][0].F = mattr[i][0].X = [i][0].Y = [i][0].M = 0;
		for (size_t j = 0; j < qseq->len; j++) {
			const double M05 = M[j];
			const double M4 = M[j-1];
			M13 = M3[j];
			M12 = M2[j];
			M11 = M1[j];

			if (!sp->mtx)
				v = SCORE(dseq->seq[i], qseq->seq[j], 1.0, -1.0);
			else
				v = VDTABLE(dseq->seq[i], qseq->seq[j]);

			M[j] = MAX(MAX(MAX(0,MAX(MAX(MAX(M13, F[j]), X[j]), pY)+ v), M12 + g->matchmax5), M11 + g->matchmax6);
			/*tmp.M = match + MAX(MAX(MAX(MAX([i - 3][j - 1].F, [i - 3][j - 1].X), [i - 3][j - 1].Y), [i - 3][j - 1].M), 0); */
			/*mattr[i][j].M = MAX(MAX(MAX(0, tmp.M),[i - 2][j - 1].M + matchmax5), [i - 1][j - 1].M + matchmax6);*/
			global_max = MAX(global_max, M[j]);

			X[j] = MAX(MAX(M13 + g->xopen1,X[j] + g->xext1), F[j]+ g->fx1);
			/*[i][j].X = MAX(MAX([i - 3][j].X + xext1, [i - 3][j].M + xopen1),[i - 3][j].F + fx1);*/

			F[j] = MAX(M11+ g->fopen, F1[j]+ g->fext);
			/*[i][j].F = MAX([i - 1][j].M + f_open,[i - 1][j].F + fext);*/

			/* set YGap */
			Y[j] = pY = MAX(MAX(MAX(M13/*i0*/ + g->yopen1, M14/*i1*/ + g->yopen2), M15/*i2*/+ g->yopen3),pY + g->yext);
			/*	[i][j].Y = MAX(MAX(MAX(MAX([i - 1][j - 1].M + yopen2, [i][j - 1].M + yopen1),[i - 2][j - 1].M + yopen3),[i][j - 1].Y + yext),
			*/

			M15 = M05;
			M14 = M4[j];
			M13 = M3[j];
			M12 = M2[j];
			M11 = M1[j];
		}
	}

	return global_max;
}
#endif

static double fp_ms_old(search_fp_thr_profile_t *s, const sequence_t * dseq, const sequence_t * qseq)
{
	const search_fp_profile_t *sp = s->sp;
	double *ph = s->h;
	double *pe = s->e;
	double v;
	double F, global_max = 0, h11, h01;

	for (size_t i = 0; i < dseq->len; i++) {
		F = 0;
		h11 = 0;
		for (size_t j = 0; j < qseq->len; j++) {
			if (!sp->mtx)
				v = SCORE(dseq->seq[i], qseq->seq[j], 1.0, -1.0);
			else
				v = VDTABLE(dseq->seq[i], qseq->seq[j]);

			h01 = ph[j];
			ph[j] = MAX(MAX(MAX(h11, pe[j]), F) + v, 0);
			global_max = MAX(global_max, ph[j]);
			pe[j] = MAX(h11 + sp->gapOpen, pe[j] + sp->gapExt);
			F = MAX(h11 + sp->gapOpen, F + sp->gapExt);
			h11 = h01;
		}
	}

	return global_max;
}

typedef struct tag_align {
	double score;
	int x;
	int y;
} tag_align;


/* optimization of matix Frame plus algorithm */
double fp_ms_release(const search_fp_thr_profile_t * s, const sequence_t * dseq, const sequence_t * qseq)
{
	const search_fp_profile_t *sp = s->sp;

	enum state_t { s_start = 0, s_fgap, s_xgap, s_ygap, s_match, s_tmpmatch, s_end, s_state_size };
	/* Datastructure: */

	double  ***mattr;

	double score = 0.0;

	mattr = malloc((dseq->len + 1) * sizeof(double *));

	for (size_t i = 0; i < dseq->len + 1; i++) {
		mattr[i] = malloc((qseq->len + 1) * sizeof(mattr[0][0]));
		for (size_t j = 0; j < (qseq->len + 1); j++)
			mattr[i][j] = malloc(s_state_size * sizeof(mattr[0][0][0]));
	}

	for (size_t i = 0; i < dseq->len + 1; i++)
		for (size_t j = 0; j < qseq->len + 1; j++)
			mattr[i][j][s_fgap] = mattr[i][j][s_xgap] = mattr[i][j][s_ygap] = mattr[i][j][s_match] = 0;

	for (size_t i = 3; i < dseq->len + 1; i++) {
		mattr[i][0][s_fgap] = mattr[i][0][s_xgap] = mattr[i][0][s_ygap] = mattr[i][0][s_match] = 0;
		for (size_t j = 1; j < qseq->len + 1; j++) {
			/* gap penalties for line n are stored in line n-1 */

			/* Note: the code below was modified to read profile at [j][xopen1] instead of [j-1][xopen1]
			the reasons are to be found out, but this seems to be the right place to read it to be consistent with om/ms runmodes */

			const double xopen1 = sp->gapOpen;		//43
			const double xext1 = sp->gapExt;		//44
			const double fx1 = sp->gapFrame;		//45 latest
			const double matchmax5 = sp->matchMax5;	//32
			const double matchmax6 = sp->matchMax6;	//33
			const double yopen1 = sp->gapOpen;		//39
			const double yopen2 = sp->gapOpen2;		//40
			const double yopen3 = sp->gapOpen3;		//41
			const double yext = sp->gapExt;			//42
													/* constant */
			const double fext = sp->gapFrameExt;		//38
			const double f_open = sp->gapFrameOpen;	//37
			double v;
			if (!sp->mtx)
				v = SCORE(qseq->seq[j - 1], dseq->seq[i - 1], 1.0, -1.0);
			else
				v = VDTABLE(qseq->seq[j - 1], dseq->seq[i - 1]);

			const double match = v;

			/* depends on [i-1][j] */
			mattr[i][j][s_fgap] = MAX(mattr[i - 1][j][s_match] + f_open, mattr[i - 1][j][s_fgap] + fext);

			/* depends on [i-3][j] */
			mattr[i][j][s_xgap] = MAX(MAX(mattr[i - 3][j][s_xgap] + xext1, mattr[i - 3][j][s_match] + xopen1),
				/*	    xopen1,  -- match above is non-negative*/ /* from s_start */
				mattr[i - 3][j][s_fgap] + fx1);
			/* depends on [i-2..i][j-1] */
			mattr[i][j][s_ygap] = MAX(MAX(MAX(MAX(mattr[i - 1][j - 1][s_match] + yopen2,
				mattr[i][j - 1][s_match] + yopen1),
				mattr[i - 2][j - 1][s_match] + yopen3),
				mattr[i][j - 1][s_ygap] + yext),
				0 /*from s_start */);
			/* depends on [i-3][j-1] */
			mattr[i][j][s_match] = match + MAX(MAX(MAX(MAX(mattr[i - 3][j - 1][s_fgap],
				mattr[i - 3][j - 1][s_xgap]),
				mattr[i - 3][j - 1][s_ygap]),
				mattr[i - 3][j - 1][s_match]),
				0);
			mattr[i][j][s_match] = MAX(MAX(MAX(0,
				mattr[i][j][s_match]),
				mattr[i - 2][j - 1][s_match] + matchmax5),
				mattr[i - 1][j - 1][s_match] + matchmax6);

			if (score < mattr[i][j][s_match])
				score = mattr[i][j][s_match];
		}
	}
	for (size_t i = 0; i < dseq->len + 1; i++)
	{
		for (size_t j = 0; j < qseq->len + 1; j++)
			free(mattr[i][j]);
		free(mattr[i]);
	}
	return score;
}

/* optimization of matix Frame plus algorithm */
tag_align score_frameplus_p2n_opt2(const search_fp_profile_t * sp, const sequence_t * dseq, const sequence_t * qseq)
{
	/* What I know about the sequence:
	1. before we do anything, nucleic seq. is prepened by a blank symbol.
	2. Next, each position [i] is subst. with AA corresponding to [i-2,i-1, i]. It follows that the first 3 positions are always 'X' (no corr. AA).
	3. Stop codon has code 27 (corr to '*' in blosum & oth matrices)
	*/
	enum state_t { s_start = 0, s_fgap, s_xgap, s_ygap, s_match, s_tmpmatch, s_end, s_state_size };
	/* Datastructure: */

	double  ***mattr;

	tag_align max_v = (tag_align) { 0.0, -1, -1 }; /* best seen alignment score in a match state */

	mattr = malloc((dseq->len + 1) * sizeof(double *));

	for (size_t i = 0; i < dseq->len + 1; i++) {
		mattr[i] = malloc((qseq->len + 1) * sizeof(mattr[0][0]));
		for (size_t j = 0; j < (qseq->len + 1); j++)
			mattr[i][j] = malloc(s_state_size * sizeof(mattr[0][0][0]));
	}

	for (size_t i = 0; i < dseq->len + 1; i++)
		for (size_t j = 0; j < qseq->len + 1; j++)
			mattr[i][j][s_fgap] = mattr[i][j][s_xgap] = mattr[i][j][s_ygap] = mattr[i][j][s_match] = 0;

	for (size_t i = 3; i < dseq->len + 1; i++) {
		mattr[i][0][s_fgap] = mattr[i][0][s_xgap] = mattr[i][0][s_ygap] = mattr[i][0][s_match] = 0;
		for (size_t j = 1; j < qseq->len + 1; j++) {
			/* gap penalties for line n are stored in line n-1 */

			/* Note: the code below was modified to read profile at [j][xopen1] instead of [j-1][xopen1]
			the reasons are to be found out, but this seems to be the right place to read it to be consistent with om/ms runmodes */

			const double xopen1 = sp->gapOpen;		//43
			const double xext1 = sp->gapExt;		//44
			const double fx1 = sp->gapFrame;		//45 latest
			const double matchmax5 = sp->matchMax5;	//32
			const double matchmax6 = sp->matchMax6;	//33
			const double yopen1 = sp->gapOpen;		//39
			const double yopen2 = sp->gapOpen2;		//40
			const double yopen3 = sp->gapOpen3;		//41
			const double yext = sp->gapExt;			//42
													/* constant */
			const double fext = sp->gapFrameExt;		//38
			const double f_open = sp->gapFrameOpen;	//37
			double v;
			if (!sp->mtx)
				v = SCORE(qseq->seq[j - 1], dseq->seq[i - 1], 1.0, -1.0);
			else
				v = VDTABLE(qseq->seq[j - 1], dseq->seq[i - 1]);

			const double match = v;

			/* depends on [i-1][j] */
			mattr[i][j][s_fgap] = MAX(mattr[i - 1][j][s_match] + f_open, mattr[i - 1][j][s_fgap] + fext);

			/* depends on [i-3][j] */
			mattr[i][j][s_xgap] = MAX(MAX(mattr[i - 3][j][s_xgap] + xext1, mattr[i - 3][j][s_match] + xopen1),
				/*	    xopen1,  -- match above is non-negative*/ /* from s_start */
				mattr[i - 3][j][s_fgap] + fx1);
			/* depends on [i-2..i][j-1] */
			mattr[i][j][s_ygap] = MAX(MAX(MAX(MAX(mattr[i - 1][j - 1][s_match] + yopen2,
				mattr[i][j - 1][s_match] + yopen1),
				mattr[i - 2][j - 1][s_match] + yopen3),
				mattr[i][j - 1][s_ygap] + yext),
				0 /*from s_start */);
			/* depends on [i-3][j-1] */
			mattr[i][j][s_match] = match + MAX(MAX(MAX(MAX(mattr[i - 3][j - 1][s_fgap],
				mattr[i - 3][j - 1][s_xgap]),
				mattr[i - 3][j - 1][s_ygap]),
				mattr[i - 3][j - 1][s_match]),
				0);
			mattr[i][j][s_match] = MAX(MAX(MAX(0,
				mattr[i][j][s_match]),
				mattr[i - 2][j - 1][s_match] + matchmax5),
				mattr[i - 1][j - 1][s_match] + matchmax6);

			if (max_v.score < mattr[i][j][s_match])
				max_v = (tag_align) { mattr[i][j][s_match], i, j };
		}
	}
	for (size_t i = 0; i < dseq->len + 1; i++)
	{
		for (size_t j = 0; j < qseq->len + 1; j++)
			free(mattr[i][j]);
		free(mattr[i]);
	}
	return max_v;
}


double fp_thr(search_fp_thr_profile_t * sp, const sequence_t *dseq, const sequence_t *qseq) {
	search_fp_profile_t *s = sp->sp;
	memset(sp->h, 0, sizeof(double) * (qseq->len));
	memset(sp->e, 0, sizeof(double) * (qseq->len));
	return fp_ms_release(sp, dseq, qseq);
}


search_fp_thr_profile_t * search_fp_thr_init(search_fp_profile_t *s, size_t thr)
{
	search_fp_thr_profile_t * sthr = malloc(sizeof(search_fp_thr_profile_t) * thr);
	for (size_t i = 0; i < thr; i++) {
			sthr[i].h = malloc(sizeof(double) * s->max_query_len);
			sthr[i].e = malloc(sizeof(double) * s->max_query_len);
			sthr[i].sp = s;
	}
	return sthr;
}

void search_fp_thr_deinit(search_fp_thr_profile_t * sthr, size_t thr)
{
	for (size_t i = 0; i < thr; i++) {
		free(sthr[i].h);
		free(sthr[i].e);
	}
	free(sthr);
}
