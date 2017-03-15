/*
Copyright (C) 2016

Contact: Dmitry Sigaev <dima.sigaev@gmail.com>
*/

#include "../../utests/tests.h"

#include <stdint.h>
#include <stdio.h>
#include "../lal_tables.h"
#include "../lal_scoring_matrix.h"
#include "../lal_report.h"


char badmatirx[] = { "#  Matrix contains errors \n \
#  number of columns mach more then 32    \n \
   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  * A  R  N  D  C  Q  E  G  H  I  L    \n \
A  8 -3 -4 -5 -2 -2 -3 -1 -4 -4 -4 -2 -3 -5 -2  1 -1 -6 -5 -2 -4 -2 -2 -10   \n \
" };
START_TEST(test_bad_matrix)
{
	scoring_matrix_t mtx;
	int status = read_scoring_matrix(&mtx, badmatirx, strlen(badmatirx));
	ck_assert_int_eq(status, 0); /* bad format*/
	report("bad format reported");
	ck_assert_int_eq((int)mtx.scale, -1);
	free_scoring_matrix(&mtx);

}END_TEST

START_TEST(test_gaptest1_matrix)
{
	scoring_matrix_t mtx;
	int status = read_scoring_matrix(&mtx, gaptest1, strlen(gaptest1));
	//	print_matrix(&mtx.sc_int_matrix);
	ck_assert_int_eq(status, 1); /* ok */
	ck_assert_int_eq((int)mtx.scale, 1);
	ck_assert_int_eq(mtx.sc_int_matrix.idata[0][0], 2); /* in the upper-left corner */
	free_scoring_matrix(&mtx);

}END_TEST

START_TEST(test_read_scoring_matrix)
{
	scoring_matrix_t mtx;
	int status = read_scoring_matrix(&mtx, identity_nuc, strlen(identity_nuc));
	ck_assert_int_eq(status, 1); /* ok */
	//	print_matrix(&mtx.sc_double_matrix);
	//	print_matrix(&mtx.sc_int_matrix);
	ck_assert_int_eq(mtx.sc_int_matrix.idata[0][0], 5); /* in the upper-left corner */
	free_scoring_matrix(&mtx);
	status = read_scoring_matrix(&mtx, blosum100, strlen(blosum100));
	ck_assert_int_eq(status, 1); /* ok */
	ck_assert_int_eq(mtx.sc_int_matrix.idata[27][27], 1); /* *x* in the bottom-rignt corner */
	ck_assert_int_eq(mtx.sc_int_matrix.idata[25][27], -10); /* zx*  */
	free_scoring_matrix(&mtx);
}END_TEST

void addScoringMatrixTC(Suite *s) {
	TCase *tc_core = tcase_create("Scoring Matrix");
	tcase_add_test(tc_core, test_gaptest1_matrix);
	tcase_add_test(tc_core, test_bad_matrix);
	tcase_add_test(tc_core, test_read_scoring_matrix);

	suite_add_tcase(s, tc_core);
}

