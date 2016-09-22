/*
Copyright (C) 2016

Contact: Dmitry Sigaev <dima.sigaev@gmail.com>
*/

#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "lal_matrix.h"
#include "lal_tables.h"
#include "lal_scoring_matrix.h"


int8_t * score_matrix_8 = NULL; // char
int16_t * score_matrix_16 = NULL; // short
int64_t * score_matrix_64 = NULL; // long
double * score_matrix_d = NULL; // long
#define SCORE_MATRIX_64(x, y) (score_matrix_64[(x << 5) + y])
#define SCORE_MATRIX_8(x, y) (score_matrix_8[(x << 5) + y])
#define SCORE_MATRIX_16(x, y) (score_matrix_16[(x << 5) + y])
#define SCORE_MATRIX_D(x, y) (score_matrix_d[(x << 5) + y])


 /**
  * Allocates memory for storing the scoring matrices.
  *
  * score_matrix_64 is initialised with -1.
  */
static void prepare_matrices() {

	score_matrix_8 = malloc(SCORE_MATRIX_DIM * SCORE_MATRIX_DIM * sizeof(int8_t));
	score_matrix_16 = malloc(SCORE_MATRIX_DIM * SCORE_MATRIX_DIM * sizeof(int16_t));
	score_matrix_64 = malloc(SCORE_MATRIX_DIM * SCORE_MATRIX_DIM * sizeof(int64_t));
	score_matrix_d = malloc(SCORE_MATRIX_DIM * SCORE_MATRIX_DIM * sizeof(double));
	memset(score_matrix_64, -1, SCORE_MATRIX_DIM * SCORE_MATRIX_DIM * sizeof(int64_t));
}

/**
* Copies the values read into score_matrix_64 into the other representations
* of the scoring matrices.
*/
static void finalize_matrices() {
	int64_t sc, lo, hi;

	hi = -100;
	lo = 100;

	for (size_t a = 0; a < SCORE_MATRIX_DIM; a++) {
		for (size_t b = 0; b < SCORE_MATRIX_DIM; b++) {
			sc = SCORE_MATRIX_64(a, b);
			if (sc < lo)
				lo = sc;
			if (sc > hi)
				hi = sc;
		}
	}

	for (size_t a = 0; a < SCORE_MATRIX_DIM; a++) {
		for (size_t b = 0; b < SCORE_MATRIX_DIM; b++) {
			sc = SCORE_MATRIX_64(a, b);

			SCORE_MATRIX_8(a, b) = (int8_t)sc;
			SCORE_MATRIX_16(a, b) = (int16_t)sc;
		}
	}
}

void gencore_mat_init_constant_scoring(scoring_matrix_t *scmat) {
	prepare_matrices();
	for (size_t a = 0; a < SCORE_MATRIX_DIM; a++) { // 0
		for (size_t b = 0; b < SCORE_MATRIX_DIM; b++) { // 0
			SCORE_MATRIX_64(a, b) = (int64_t)scmat->sc_int_matrix.ddata[a][b];
			SCORE_MATRIX_D(a, b) = (double)scmat->sc_double_matrix.ddata[a][b];
		}
	}

	finalize_matrices();
}

static double min_double(double a, double b) {
	return (a < b) ? a : b;
}

static double max_double(double a, double b) {
	return (a > b) ? a : b;
}

static double get_second_max_neg(double  a1, double a2, double x) {
	if (x >= 0) return a2;
	if (x >= a1) /* a1 <= x < 0 */ return a1;
	if (x >= a2) /* a2 <= x < a1 */ return x;
	return a2; /* x < a2 <= a1 < 0 */
}

static void update(double * second_max_negative, double * max_negative, double * min_positive, double val) {
	const double tolerance = 0.000001;
	if (val < -tolerance) {
		*second_max_negative = get_second_max_neg(*max_negative, *second_max_negative, val);
		*max_negative = max_double(*max_negative, val);
	}
	else if (val > tolerance)
		*min_positive = min_double(*min_positive, val);
}

typedef struct tag_descritor
{
	char *data;
	char **list;
	size_t len;
}descritor_t;

descritor_t create_descritor(const char *matrixstring, size_t len)
{
	descritor_t dsc;
	dsc.data = malloc(len + 1);
	memcpy(dsc.data, matrixstring, len + 1); /*+1 \null simbol for string data*/
	char *Letters = strtok(dsc.data, "\t\n");
	size_t rows = 0;
	while (Letters) {
		Letters = strtok(NULL, "\t\n");
		rows++;
	}
	dsc.len = rows;
	dsc.list = (char **)malloc(rows * sizeof(char*));
	dsc.list[0] = dsc.data;
	for (size_t i = 0; i < dsc.len - 1 /*dispose of heap corruption */; i++) {
		dsc.list[i + 1] = dsc.list[i] + strlen(dsc.list[i]) + 1;
	}
	return dsc;
}

void free_descritor(descritor_t * desc) {
	desc->list[0] = NULL;
	if (desc->data)
		free(desc->data); // should be first;
	desc->data = NULL;
	if (desc->list)
		free(desc->list);
	*desc = (descritor_t) { NULL, NULL, 0 };
}

void read_scoring_matrix(scoring_matrix_t *mtx, const char *matrixstring, size_t len) {
	descritor_t  desc = create_descritor(matrixstring, len);
	const char *s;
	const char *pdest;
	size_t cnt, next_line;
	char *Letter, *Number;
	int Columns[32], Row, Count;
	float Val;
	int rows_count, columns_count;
	int i_row, i_col;
	/* scaling */
	const double big_num = 100000000.0;
	const float data4vec_infty = INT_MIN;

	double min_positive = big_num; /* Minimal positive element  */
	double second_max_negative = -big_num;
	double max_negative = -big_num; /* negative closest to 0 */

	double result = -1;
	/* return scale based on analysis of the profile
	The scale is not optimal but works around the following ideas:
	the smallest (in abs. numbers) penalty is scaled at least to '-1', so we still penalize it.
	No we look at the smallest positive score (after the scale) and make sure it's at least 1 as well (so we don't increase it's score too much by rounding).
	*/
	/*end*/

	mtx->sc_double_matrix = matrix(SCORE_MATRIX_DIM, SCORE_MATRIX_DIM, DOUBLETYPE);
	mtx->sc_int_matrix = matrix(SCORE_MATRIX_DIM, SCORE_MATRIX_DIM, INTTYPE);
	mtx->scale = 1.0;
	mtx->scaleback = 1.0;

	matrix_set_value(&mtx->sc_double_matrix, (element_t) { 0, DOUBLETYPE });
	matrix_set_value(&mtx->sc_int_matrix, (element_t) { 0, INTTYPE });

	/* Skip the comment */
	mtx->Doc[0] = '\0';
	cnt = 0;
	for (size_t i = 0; i < desc.len; i++) {
		if (s = desc.list[i]) {
			if (!(pdest = strstr(s, "#"))) {
				next_line = i;
				break;
			}
			if ((cnt += strlen(pdest)) + 1 < MAX_DOC_LEN)
				strcat(mtx->Doc, pdest);
		}

	}

	/* Parse the column description line */
	Letter = strtok(desc.list[next_line], " \t\n");
	columns_count = 0;
	while (Letter) {
		Columns[columns_count] = lal_encode31[*Letter];
		Letter = strtok(NULL, " \t\n");
		columns_count++;
	}

	/* Read the table - next step */
	rows_count = 0;
	for (size_t i = next_line + 1; i < desc.len; i++)
		if (s = desc.list[i]) {
			/* Read the row description letter */
			if (!(Letter = strtok((char *)s, " \t\n"))) continue; /* ignore empty lines */
			Row = lal_encode31[*Letter];

			Count = 0;
			while (Number = strtok(NULL, " \t\n")) {
				if (sscanf(Number, "%f", &Val) != 1) {
					printf("Bad type of comparison matrix \n");
					break;
				}
				/* ncbi enable rectangular (not square) matrix. We are make it simmetric. */
				mtx->sc_double_matrix.ddata[Row][Columns[Count]] = Val;
				if (Val != data4vec_infty)
					update(&second_max_negative, &max_negative, &min_positive, Val);
				Count++;
			}
			rows_count++;

		}
	free_descritor(&desc);
	/* calucate scaling factor and integer matrix*/
	mtx->man2mip[0] = second_max_negative;
	mtx->man2mip[1] = max_negative;
	mtx->man2mip[2] = min_positive;
	assert(max_negative < 0);
	if (min_positive > -max_negative) {
		result = rint(-second_max_negative / max_negative) / second_max_negative; /* |max_negative| < |min_positive|  */
	}
	else {
		result = rint(max_negative / min_positive) / max_negative;
	}

	/* The following part is written on order to take care of rectangular matrices that are
	   not square. In such a case the matrix is forced to be symmetric. This is not
	   necessarily the right thing to be done but this is the convention here */

	if (rows_count != columns_count) { /* TODO: unit tests */
		if (rows_count > columns_count) {
			for (i_row = 0; i_row < mtx->sc_double_matrix.nrows; i_row++) {
				for (i_col = 0; i_col < i_row; i_col++) {
					mtx->sc_double_matrix.ddata[i_col][i_row] = mtx->sc_double_matrix.ddata[i_row][i_col];
					mtx->sc_int_matrix.idata[i_row][i_col] = lrint(mtx->sc_double_matrix.ddata[i_row][i_col] * result);
				}
			}
		}
		else {
			for (i_col = 0; i_col < mtx->sc_double_matrix.ncols; i_col++) {
				for (i_row = i_col + 1; i_row < mtx->sc_double_matrix.nrows; i_row++) {
					mtx->sc_double_matrix.ddata[i_row][i_col] = mtx->sc_double_matrix.ddata[i_col][i_row];
					mtx->sc_int_matrix.idata[i_row][i_col] = lrint(mtx->sc_double_matrix.ddata[i_col][i_row] * result);
				}
			}
		}
	}
	else {
		for (i_row = 0; i_row < mtx->sc_double_matrix.nrows; i_row++)
			for (i_col = 0; i_col < mtx->sc_double_matrix.ncols; i_col++)
				mtx->sc_int_matrix.idata[i_col][i_row] = lrint(mtx->sc_double_matrix.ddata[i_row][i_col] * result);
	}

	mtx->scale = result;
	mtx->scaleback = 1.0 / result;
}


void update_gap_scoring_matrix(scoring_matrix_t *mtx, double gapOpen, double gapExtend)
{
	double result = -1;
	update(mtx->man2mip, mtx->man2mip + 1, mtx->man2mip + 2, gapOpen);
	update(mtx->man2mip, mtx->man2mip + 1, mtx->man2mip + 2, gapExtend);
	update(mtx->man2mip, mtx->man2mip + 1, mtx->man2mip + 2, gapOpen);
	update(mtx->man2mip, mtx->man2mip + 1, mtx->man2mip + 2, gapExtend);
	assert(mtx->man2mip[1] < 0);
	if (mtx->man2mip[2] > -mtx->man2mip[1]) {
		result = rint(-mtx->man2mip[0] / mtx->man2mip[1]) / mtx->man2mip[0]; /* |max_negative| < |min_positive|  */
	}
	else {
		result = rint(mtx->man2mip[1] / mtx->man2mip[2]) / mtx->man2mip[1];
	}
	for (size_t i_row = 0; i_row < mtx->sc_double_matrix.nrows; i_row++)
		for (size_t i_col = 0; i_col < mtx->sc_double_matrix.ncols; i_col++)
			mtx->sc_int_matrix.idata[i_col][i_row] = lrint(mtx->sc_double_matrix.ddata[i_row][i_col] * result);
	mtx->scale = result;
	mtx->scaleback = 1.0 / result;
	mtx->gapOpen = gapOpen;
	mtx->gapExtend = gapExtend;
	printf("s=%f, go=%f, ge=%f \n", result, gapOpen, gapExtend);
}

void free_scoring_matrix(scoring_matrix_t *mtx)
{
	free_matrix(&mtx->sc_double_matrix);
	free_matrix(&mtx->sc_int_matrix);
}
