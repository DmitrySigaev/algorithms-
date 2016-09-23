/*
Copyright (C) 2016

Contact: Dmitry Sigaev <dima.sigaev@gmail.com>
*/

#include "../../utests/tests.h"

#include <stdint.h>
#include <stdio.h>
#include "../lal_tables.h"
#include "../lal_scoring_matrix.h"


START_TEST(test_read_scoring_matrix)
{
	scoring_matrix_t mtx;
	read_scoring_matrix(&mtx, identity_nuc, strlen(identity_nuc));
//	print_matrix(&mtx.sc_double_matrix);
//	print_matrix(&mtx.sc_int_matrix);
	ck_assert_int_eq(mtx.sc_int_matrix.idata[0][0], 5); /* in the upper-left corner */
	free_scoring_matrix(&mtx);
}END_TEST

void addScoringMatrixTC(Suite *s) {
	TCase *tc_core = tcase_create("Scoring Matrix");
	tcase_add_test(tc_core, test_read_scoring_matrix);

	suite_add_tcase(s, tc_core);
}

