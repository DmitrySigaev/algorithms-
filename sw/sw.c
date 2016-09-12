#include <stdio.h>
#include <malloc.h>
#include "sw.h"

#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define SCORE(a, b, c, d) (((a) == (b)) ? (c) : (d))
/*
 * Sequence alignments are used in different area of computer science.
 * Main feature of alignment is a gap function because complexity of
 * computation depends on this function. The most common and simple
 * case is using a constant gap penalties.
 */

double sw_constant_gap(const search_swcd_profile_t * sp, const sequence_t * dseq, const sequence_t * qseq)
{
	double d_last, u_last, l_last;
	double d_new, u_new, l_new;
	double ** score_mat = (double **)malloc((dseq->len + 1)* sizeof(double *));
	if (NULL == score_mat) {
		return -1.0;
	}
	for (size_t i = 0; i <= dseq->len; i++) {
		score_mat[i] = (double *)malloc((qseq->len + 1) * sizeof(double));
		if (NULL == score_mat[i]) {
			for (size_t j = 0; j < i; j++) {
				free(score_mat[j]);
			}
			free(score_mat);
			return -1.0;
		}
	}

	for (size_t i = 0; i <= dseq->len; i++) {
		for (size_t j = 0; j <= qseq->len; j++) {
			// This is the first row / column which is all zeros
			if (i == 0 || j == 0) {
				score_mat[i][j] = 0;
				continue;
			}
			else {
				d_last = score_mat[i - 1][j - 1];
				u_last = score_mat[i - 1][j];
				l_last = score_mat[i][j - 1];
			}
			d_new = d_last + SCORE(dseq->seq[i-1], dseq->seq[j-1], 1.0, -1.0);
			u_new = u_last + sp->gap;
			l_new = l_last + sp->gap;
			score_mat[i][j] = MAX(MAX(d_new, u_new), MAX(l_new, 0));
		}
	}

	double score = score_mat[0][0];
	for (size_t i = 0; i <= dseq->len; i++) {
		for (size_t j = 0; j <= qseq->len; j++) {
			printf(" %f", score_mat[i][j]);
			if (score < score_mat[i][j]) {
				score = score_mat[i][j];
			}
		}
		printf(" \n");
	}

	for (size_t i = 0; i < dseq->len; i++) {
		free(score_mat[i]);
	}
	free(score_mat);
	return score;
}


