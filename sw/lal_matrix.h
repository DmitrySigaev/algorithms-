/*
Copyright (C) 2016 Dmitry Sigaev

Contact: Dmitry Sigaev <dima.sigaev@gmail.com>
*/


#ifndef _LAL_MARTIX_H_
#define _LAL_MARTIX_H_

#include <stdint.h>

enum MATTYPE { VOIDTYPE, DOUBLETYPE, INTTYPE, FLOATTYPE };
/**
 *  Structure of a matrix[nrows][ncols]
 *  @typedef {struct} matrix_t
 *  @field data  matrix data
 *  @field nrows the number of rows
 *  @field ncols the number of columns
 */
typedef struct tag_matrix {
	union {
		double **ddata;
		int64_t **idata;
	};
	size_t nrows;
	size_t ncols;
	enum MATTYPE type;
} matrix_t;

typedef struct tag_element {
	union {
		double d;
		int64_t i;
	};
	enum MATTYPE type;
} element_t;


matrix_t matrix(const size_t nrows, const size_t ncols, enum MATTYPE type);
element_t find_max(const matrix_t *matrix);
int matrix_set_value(matrix_t *matrix, element_t value);
void print_matrix(const matrix_t *matrix);
void free_matrix(matrix_t *matrix);

#endif /* _LAL_MARTIX_H_ */