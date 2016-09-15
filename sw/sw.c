#include <stdio.h>
#include "lal_matrix.h"
#include "sw.h"

#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define SCORE(a, b, c, d) (((a) == (b)) ? (c) : (d))

/*
 * Sequence alignments are used in different area of computer science.
 * Main feature of alignment is a gap function because complexity of
 * computation depends on this function. The most common and simple
 * case is using a constant gap penalties.
 */
double sw_constant_gap(const search_swcg_profile_t * sp, const sequence_t * dseq, const sequence_t * qseq)
{
	double d_last, u_last, l_last;
	double d_new, u_new, l_new;
	matrix_t score_mat = matrix(dseq->len + 1, qseq->len + 1);

	for (size_t i = 0; i <= dseq->len; i++) {
		for (size_t j = 0; j <= qseq->len; j++) {
			// This is the first row / column which is all zeros
			if (i == 0 || j == 0) {
				score_mat.data[i][j] = 0;
				continue;
			}
			else {
				d_last = score_mat.data[i - 1][j - 1];
				u_last = score_mat.data[i - 1][j];
				l_last = score_mat.data[i][j - 1];
			}
			d_new = d_last + SCORE(dseq->seq[i-1], dseq->seq[j-1], 1.0, -1.0);
			u_new = u_last + sp->gap;
			l_new = l_last + sp->gap;
			score_mat.data[i][j] = MAX(MAX(d_new, u_new), MAX(l_new, 0));
		}
	}

	double score = find_max(&score_mat);
	print_matrix(&score_mat);
	free_matrix(&score_mat);
	return score;
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

	matrix_t score_mat = matrix(dseq->len + 1, qseq->len + 1);
	matrix_t ee = matrix(dseq->len + 1, qseq->len + 1);
	matrix_t ff = matrix(dseq->len + 1, qseq->len + 1);

	for (size_t i = 0; i <= dseq->len; i++) {
		for (size_t j = 0; j <= qseq->len; j++) {
			// This is the first row / column which is all zeros
			if (i == 0 || j == 0) {
				if (i == 0)
					ee.data[i][j] = 0;
				if (j == 0)
					ff.data[i][j] = 0;
				score_mat.data[i][j] = 0;
				continue;
			}
			else {
				e_last = ee.data[i][j - 1];
				f_last = ff.data[i - 1][j];
				d_last = score_mat.data[i - 1][j - 1];
				u_last = score_mat.data[i - 1][j];
				l_last = score_mat.data[i][j - 1];
			}
			d_new = d_last + SCORE(dseq->seq[i - 1], dseq->seq[j - 1], 1.0, -1.0);
			u_new = u_last + sp->gapOpen;
			l_new = l_last + sp->gapOpen;
			e_new = e_last + sp->gapExt;
			f_new = f_last + sp->gapExt;
			ee.data[i][j] = MAX(e_new, l_new);
			ff.data[i][j] = MAX(f_new, u_new);
			score_mat.data[i][j] = MAX(MAX(d_new, ee.data[i][j]), MAX(ff.data[i][j], 0));
		}
	}

	double score = find_max(&score_mat);
	print_matrix(&score_mat);
	print_matrix(&ee);
	print_matrix(&ff);
	free_matrix(&score_mat);
	free_matrix(&ee);
	free_matrix(&ff);
	return score;
}