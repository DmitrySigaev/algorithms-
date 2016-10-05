/*
Copyright (C) 2016

Contact: Dmitry Sigaev <dima.sigaev@gmail.com>
*/

#include "../../utests/tests.h"

#include <stdint.h>
#include <stdio.h>
#include <malloc.h>
#include "../sw.h"
#include "../lal_encoding.h"
#include "../lal_tables.h"

START_TEST(test_sw_double_symbols)
{
	char seq1[] = { "ABAC" };
	size_t len1 = strlen(seq1);
	char seq2[] = { "ADACCG" };
	size_t len2 = strlen(seq2);
	sequence_t sseq1 = { 1, (char *)seq1, len1 };
	sequence_t sseq2 = { 2, (char *)seq2, len2 };
	search_swcg_profile_t sp = { -1, NULL };
	double score = sw_constant_gap_double(&sp, &sseq1, &sseq2);
	ck_assert_int_eq((int)score, 4); /* Max score */
}END_TEST


START_TEST(test_sw_double_encoded)
{
	char seq1[] = { "ABAC" };
	size_t len1 = strlen(seq1);
	char seq2[] = { "ADACCG" };
	size_t len2 = strlen(seq2);
	sequence_t inseq1 = { 1, (char *)seq1, len1 };
	sequence_t inseq2 = { 2, (char *)seq2, len2 };
	sequence_t enseq1 = { 1, malloc(len1 + 1), len1 };
	sequence_t enseq2 = { 2, malloc(len2 + 1), len2 };
	seq2encodedseq(inseq1, enseq1, lal_encode31);
	seq2encodedseq(inseq2, enseq2, lal_encode31);
	search_swcg_profile_t sp = { -1, NULL };
	double score = sw_constant_gap_double(&sp, &enseq1, &enseq2);
	ck_assert_int_eq((int)score, 4); /* Max score */
}END_TEST

START_TEST(test_sw_double_encoded_vtable)
{
	scoring_matrix_t mtx;
	int status = read_scoring_matrix(&mtx, gaptest1, strlen(gaptest1));
	char seq1[] = { "ABAC" };
	size_t len1 = strlen(seq1);
	char seq2[] = { "ADACCG" };
	size_t len2 = strlen(seq2);
	sequence_t inseq1 = { 1, (char *)seq1, len1 };
	sequence_t inseq2 = { 2, (char *)seq2, len2 };
	sequence_t enseq1 = { 1, malloc(len1 + 1), len1 };
	sequence_t enseq2 = { 2, malloc(len2 + 1), len2 };
	seq2encodedseq(inseq1, enseq1, lal_encode31);
	seq2encodedseq(inseq2, enseq2, lal_encode31);
	search_swcg_profile_t sp = { -1, (!status) ? (NULL) : (&mtx) };
	double score = sw_constant_gap_double(&sp, &enseq1, &enseq2);
	ck_assert_int_eq((int)score, 8); /* Max score */
}END_TEST


START_TEST(test_sw_double)
{
	char seq1[] = { 1 , 2 , 1 , 5 };
	size_t len1 = sizeof(seq1) / sizeof(seq1[0]);
	char seq2[] = { 1 , 3 , 1 , 5, 5, 6 };
	size_t len2 = sizeof(seq2) / sizeof(seq2[0]);
	sequence_t sseq1 = { 1, (char *)seq1, len1 };
	sequence_t sseq2 = { 2, (char *)seq2, len2 };
	search_swcg_profile_t sp = { -1, NULL };
	double score = sw_constant_gap_double(&sp, &sseq1, &sseq2);
	ck_assert_int_eq((int)score, 4); /* Max score */

}END_TEST

START_TEST(test_sw_int)
{
	char seq1[] = { 1 , 2 , 1 , 5 };
	size_t len1 = sizeof(seq1) / sizeof(seq1[0]);
	char seq2[] = { 1 , 3 , 1 , 5, 5, 6 };
	size_t len2 = sizeof(seq2) / sizeof(seq2[0]);
	sequence_t sseq1 = { 1, (char *)seq1, len1 };
	sequence_t sseq2 = { 2, (char *)seq2, len2 };
	search_swcg_profile_int_t sp = { -1, NULL };
	int64_t score = sw_constant_gap_int(&sp, &sseq1, &sseq2);
	ck_assert_int_eq(score, 4); /* Max score */

}END_TEST

void addSWTC(Suite *s) {
	TCase *tc_core = tcase_create("SW");
	tcase_add_test(tc_core, test_sw_double_encoded_vtable);
	tcase_add_test(tc_core, test_sw_double_encoded);
	tcase_add_test(tc_core, test_sw_double_symbols);
	tcase_add_test(tc_core, test_sw_double);
	tcase_add_test(tc_core, test_sw_int);

	suite_add_tcase(s, tc_core);
}

