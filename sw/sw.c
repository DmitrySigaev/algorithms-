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
double sw_constant_gap(const search_swcd_profile_t * sp, const sequence_t * dseq, const sequence_t * qseq)
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
