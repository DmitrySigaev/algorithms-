/*
Copyright (C) 2016

Contact: Dmitry Sigaev <dima.sigaev@gmail.com>
*/

#include <stdio.h>
#include <malloc.h>
#include "lal_matrix.h"
#include "lal_report.h"
#include "sw.h"

#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define SCORE(a, b, c, d) (((a) == (b)) ? (c) : (d))
#define VDTABLE(a, b) (sp->mtx->sc_double_matrix.ddata[(a)][(b)])
#define VITABLE(a, b) (sp->mtx->sc_int_matrix.idata[(a)][(b)])
/*
 * Sequence alignments are used in different area of computer science.
 * Main feature of alignment is a gap function because complexity of
 * computation depends on this function. The most common and simple
 * case is using a constant gap penalties.
 */
double sw_constant_gap_double(const search_swcg_profile_t * sp, const sequence_t * dseq, const sequence_t * qseq)
{
	double d_last, u_last, l_last;
	double d_new, u_new, l_new;
	double v;
	matrix_t score_mat = matrix(dseq->len + 1, qseq->len + 1, DOUBLETYPE);
	for (size_t i = 0; i <= dseq->len; i++) {
		for (size_t j = 0; j <= qseq->len; j++) {
			// This is the first row / column which is all zeros
			if (i == 0 || j == 0) {
				score_mat.ddata[i][j] = 0;
				continue;
			}
			else {
				d_last = score_mat.ddata[i - 1][j - 1];
				u_last = score_mat.ddata[i - 1][j];
				l_last = score_mat.ddata[i][j - 1];
			}
			if (!sp->mtx)
				v = SCORE(dseq->seq[i - 1], qseq->seq[j - 1], 1.0, -1.0);
			else
				v = VDTABLE(dseq->seq[i - 1], qseq->seq[j - 1]);
			d_new = d_last + v;
			u_new = u_last + sp->gap;
			l_new = l_last + sp->gap;
			score_mat.ddata[i][j] = MAX(MAX(d_new, u_new), MAX(l_new, 0));
		}
	}

	element_t score = find_max(&score_mat);
//	print_matrix(&score_mat);
	free_matrix(&score_mat);
	return score.d;
}

/*
* Sequence alignments are used in different area of computer science.
* Main feature of alignment is a gap function because complexity of
* computation depends on this function. The most common and simple
* case is using a constant gap penalties.
*/
int64_t sw_constant_gap_int(const search_swcg_profile_int_t * sp, const sequence_t * dseq, const sequence_t * qseq)
{
	int64_t d_last, u_last, l_last;
	int64_t d_new, u_new, l_new;
	int64_t v;
	matrix_t score_mat = matrix(dseq->len + 1, qseq->len + 1, INTTYPE);

	for (size_t i = 0; i <= dseq->len; i++) {
		for (size_t j = 0; j <= qseq->len; j++) {
			// This is the first row / column which is all zeros
			if (i == 0 || j == 0) {
				score_mat.idata[i][j] = 0;
				continue;
			}
			else {
				d_last = score_mat.idata[i - 1][j - 1];
				u_last = score_mat.idata[i - 1][j];
				l_last = score_mat.idata[i][j - 1];
			}
			if (!sp->mtx)
				v = SCORE(dseq->seq[i - 1], qseq->seq[j - 1], 1, -1);
			else
				v = VITABLE(dseq->seq[i - 1], qseq->seq[j - 1]);
			d_new = d_last + v;
			u_new = u_last + sp->gap;
			l_new = l_last + sp->gap;
			score_mat.idata[i][j] = MAX(MAX(d_new, u_new), MAX(l_new, 0));
		}
	}

	element_t score = find_max(&score_mat);
//	print_matrix(&score_mat);
	free_matrix(&score_mat);
	return score.i;
}
/*
 * Sequence alignments with constant and linear gap penalties can be computed in time O(n*m) for two
 * sequences of lingth m and n. With affine gap penalties the time increases to O(n*m*(n+m)),
 * since for each cell the algorithm has to checkif a gap is extended or a new one is opened.
 * In 1982 Gotoh described a method to compute optimal sequence alignments, with affine gap
 * penalties, in time O(n*m). His version uses two additional matrices(E and F) to keep track of open gap.
 * E keeps track of gaps in the query sequence and F if gaps in the database sequence.
 * The values in E are calculated as the maximum of the previous value in H plus the costs Q for opening
 * a gap an the previous balue in E plus the costs T for extending a gap.  The values in F are
 * computed the same way, except for gaps in the other sequence. The values in H are computed like
 * in the Smith-Waterman with linear gap costs, except the value in E and F are used ad the gap costs.
 */
double sw_affine_gap(const search_swag_profile_t * sp, const sequence_t * dseq, const sequence_t * qseq)
{
	double d_last, u_last, l_last, e_last, f_last;
	double d_new, u_new, l_new, e_new, f_new;
	double v;
	matrix_t score_mat = matrix(dseq->len + 1, qseq->len + 1, DOUBLETYPE);
	matrix_t ee = matrix(dseq->len + 1, qseq->len + 1, DOUBLETYPE);
	matrix_t ff = matrix(dseq->len + 1, qseq->len + 1, DOUBLETYPE);

	for (size_t i = 0; i <= dseq->len; i++) {
		for (size_t j = 0; j <= qseq->len; j++) {
			// This is the first row / column which is all zeros
			if (i == 0 || j == 0) {
				if (i == 0)
					ee.ddata[i][j] = 0;
				if (j == 0)
					ff.ddata[i][j] = 0;
				score_mat.ddata[i][j] = 0;
				continue;
			}
			else {
				e_last = ee.ddata[i][j - 1];
				f_last = ff.ddata[i - 1][j];
				d_last = score_mat.ddata[i - 1][j - 1];
				u_last = score_mat.ddata[i - 1][j];
				l_last = score_mat.ddata[i][j - 1];
			}
			if (!sp->mtx)
				v = SCORE(dseq->seq[i - 1], qseq->seq[j - 1], 1.0, -1.0);
			else
				v = VDTABLE(dseq->seq[i - 1], qseq->seq[j - 1]);
			d_new = d_last + v;
			u_new = u_last + sp->gapOpen;
			l_new = l_last + sp->gapOpen;
			e_new = e_last + sp->gapExt;
			f_new = f_last + sp->gapExt;
			ee.ddata[i][j] = MAX(e_new, l_new);
			ff.ddata[i][j] = MAX(f_new, u_new);
			score_mat.ddata[i][j] = MAX(MAX(d_new, ee.ddata[i][j]), MAX(ff.ddata[i][j], 0));
		}
	}

	element_t score = find_max(&score_mat);
//	print_matrix(&score_mat);
//	print_matrix(&ee);
//	print_matrix(&ff);
	free_matrix(&score_mat);
	free_matrix(&ee);
	free_matrix(&ff);
	return score.d;
}

region_t sw_alignment_swipe(const search_swag_profile_t * sp, const sequence_t *xseq, const sequence_t *yseq) {
	region_t region;
	double score = 0; /* score of direct path*/
	double nick = 0; /* score of backward path*/
	size_t size = MAX(xseq->len, yseq->len);
	double *EE = malloc(size * sizeof(double));
	double *HH = malloc(size * sizeof(double));
	double v;

	/* direct aisle */
	for (size_t j = 0; j < xseq->len; j++) {
		HH[j] = 0; /* memset */
		EE[j] = sp->gapOpen;
	}

	for (size_t i = 0; i < yseq->len; i++) {
		double h = 0;
		double p = 0;
		double f = sp->gapOpen;

		for (size_t j = 0; j < xseq->len; j++) {
			f = MAX(f, h + sp->gapOpen) + sp->gapExt;
			EE[j] = MAX(EE[j], HH[j] + sp->gapOpen) + sp->gapExt;

			if (!sp->mtx)
				v = SCORE(xseq->seq[j], yseq->seq[i], 1.0, -1.0);
			else
				v = VDTABLE(xseq->seq[j], yseq->seq[i]);
			h = p + v;
			if (h < 0) h = 0;
			if (f > h) h = f;
			if (EE[j] > h) h = EE[j];
			p = HH[j];
			HH[j] = h;
			if (h > score) {
				score = h;
				region.fdscore = score;
				region.x_f = j;
				region.y_f = i;
			}
		}
	}
	/* backward pass */
	for (size_t j = 0; j < region.y_f + 1; j++) {
		HH[j] = -1;
		EE[j] = -1;
	}

	for (size_t i = 0; i < region.y_f + 1; i++) {
		size_t bi = region.y_f - i;
		double h = -1;
		double f = -1;
		double p = (i != 0) ? -1 : 0;

		for (size_t j = 0; j < region.x_f + 1; j++) {
			size_t bj = region.x_f - j;
			f = MAX(f, h + sp->gapOpen) + sp->gapExt;
			EE[bj] = MAX(EE[bj], HH[bj] + sp->gapOpen) + sp->gapExt;

			if (!sp->mtx)
				v = SCORE(xseq->seq[bj], yseq->seq[bi], 1.0, -1.0);
			else
				v = VDTABLE(xseq->seq[bj], yseq->seq[bi]);

			h = p + v;

			if (f > h) h = f;
			if (EE[bj] > h) h = EE[bj];
			p = HH[bj];
			HH[bj] = h;

			if (h > nick) {
				nick = h;
				region.x_s = bj;
				region.y_s = bi;
				if (nick >= score) {
					region.bdscore = nick;
					free(EE);
					free(HH);
					return region;
				}
			}
		}
	}
	report_warning("check in align function.");
	free(EE);
	free(HH);
	return (region_t) { 0, 0, 0, 0, 0, 0 };

}
