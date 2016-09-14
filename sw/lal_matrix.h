/*
Copyright (C) 2016 Dmitry Sigaev

Contact: Dmitry Sigaev <dima.sigaev@gmail.com>
*/


#ifndef _LAL_MARTIX_H_
#define _LAL_MARTIX_H_

/**
 *  Structure of a matrix[nrows][ncols]
 *  @typedef {struct} matrix_t
 *  @field data  matrix data
 *  @field nrows the number of rows 
 *  @field ncols the number of columns 
 */
typedef struct tag_matrix {
	double **data;
	size_t nrows;
	size_t ncols;
} matrix_t;

matrix_t matrix(const size_t nrows, const size_t ncols);
double find_max(const matrix_t *matrix);
void print_matrix(const matrix_t *matrix);
void free_matrix(matrix_t *matrix);

#endif /* _LAL_MARTIX_H_ */