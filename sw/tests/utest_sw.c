/*
Copyright (C) 2016

Contact: Dmitry Sigaev <dima.sigaev@gmail.com>
*/

#include "../../utests/tests.h"

#include <stdint.h>
#include <stdio.h>
#include "../sw.h"


START_TEST(test_sw)
{
	char seq1[] = { 1 , 2 , 1 , 5  };
	size_t len1 = sizeof(seq1) / sizeof(seq1[0]);
	char seq2[] = { 1 , 3 , 1 , 5, 5, 6 };
	size_t len2 = sizeof(seq2) / sizeof(seq2[0]);
	sequence_t sseq1 = { 1, (char *)seq1, len1 };
	sequence_t sseq2 = { 2, (char *)seq2, len2 };
	search_swcd_profile_t sp = { -1, NULL };
	double score = sw_constant_gap(&sp, &sseq1, &sseq2);
	ck_assert_int_eq((int)score, 4); /* Max score */

}END_TEST

void addSWTC(Suite *s) {
	TCase *tc_core = tcase_create("SW");
	tcase_add_test(tc_core, test_sw);

	suite_add_tcase(s, tc_core);
}

