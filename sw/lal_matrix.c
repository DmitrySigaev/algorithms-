#include <stdio.h>
#include <float.h>
#include <malloc.h>
#include "lal_matrix.h"

#define MAX(a, b) ((a) > (b) ? (a) : (b))

matrix_t matrix(const size_t nrows, const size_t ncols)
{
	double ** score_mat = (double **)malloc(nrows * sizeof(double *));
	if (NULL == score_mat) {
		return (matrix_t){ NULL, 0, 0 };
	}

	for (size_t i = 0; i < nrows; i++) {
		score_mat[i] = (double *)malloc(ncols * sizeof(double));
		if (NULL == score_mat[i]) {
			for (size_t j = 0; j < i; j++) {
				free(score_mat[j]);
			}
			free(score_mat);
			return (matrix_t){ NULL, 0, 0 };
		}
	}
	return (matrix_t) { score_mat, nrows, ncols };
}

double find_max(const matrix_t *matrix)
{
	if (!matrix) return DBL_MAX;
	if (!matrix->data) return DBL_MAX;
	double score = matrix->data[0][0];
	for (size_t i = 0; i < matrix->nrows; i++) {
		for (size_t j = 0; j < matrix->ncols; j++) {
			if (score < matrix->data[i][j]) {
				score = matrix->data[i][j];
			}
		}
	}
	return score;
}

void print_matrix(const matrix_t *matrix)
{
	if (!matrix) return;
	if (!matrix->data) return;
	for (size_t i = 0; i < matrix->nrows; i++) {
		for (size_t j = 0; j < matrix->ncols; j++) {
			printf(" %f", matrix->data[i][j]);
		}
		printf(" \n");
	}
	printf(" \n");
}

void free_matrix(matrix_t *matrix)
{
	if (!matrix)
		return;
	if (!matrix->data)
		return;
	for (size_t i = 0; i < matrix->nrows; i++) {
		free(matrix->data[i]);
	}
	free(matrix->data);
	matrix->data = NULL;
	matrix->nrows = 0;
	matrix->nrows = 0;
}

