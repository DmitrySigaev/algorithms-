/*
Copyright (C) 2017

Contact: Dmitry Sigaev <dima.sigaev@gmail.com>
*/

#include "../../utests/tests.h"

#include <stdint.h>
#include <stdio.h>
#include <malloc.h>
#include "../sw.h"
#include "../fp.h"
#include "../gc_sw.h"
#include "../gc_fp.h"
#include "../lal_encoding.h"
#include "../lal_tables.h"
#include "../lal_translate_table.h"

START_TEST(test_fp_first_refference_test)
{
	scoring_matrix_t mtx;
	int status = read_scoring_matrix(&mtx, blosum62, strlen(blosum62));

	// >ACHA_ELEEL P09688 electrophorus electricus (electric eel). acetylcholine receptor protein, alpha chain (fragment). 2/94
	char seq1[] = { "SEDETRLVKNLFSGYNKVVRPVNH" };
	size_t len1 = strlen(seq1);
	// >HSBGL2
	char seq2[] = { "ATGTCATACCTCTTATCTCCTCCCACAGCTCCTGGGCAACGTGCTGGTCTGTGTGCTGGCCCATCACTTTGGCAAAGAATTC" };
	char any_symbol = 'x';
	size_t len2 = strlen(seq2);
	sequence_t inseq1 = { 1, (char *)seq1, len1 };
	sequence_t inseq2 = { 2, (char *)seq2, len2 };
	sequence_t enseq1 = { 1, malloc(len1 + 1), len1 };
	sequence_t enseq2 = { 2, malloc(len2 + 1), len2 };
	sequence_t any = { 3, malloc(1 + 1), 1 };
	lal_seq2encodedseq(inseq1, enseq1, lal_encode31);
	lal_seq2encodedseq(inseq2, enseq2, lal_encode31);
	lal_seq2encodedseq((sequence_t){ 3, &any_symbol, 1 }, any, lal_encode31);
	mtx.scale = 10.0;
	search_fp_profile_t sp = { -10.5,/*gapOpen */\
		-0.5,/*	gapExt */\
		-10.5,  /*gapFrame*/ \
		-13.0, /* matchMax5 */ \
		-20.0, /* matchMax6 */ \
		-30.0, /*gapOpen2 */\
		-23.0, /*gapOpen3 */ \
		-7.0, /* gapFrameOpen */ \
		-13.0, /* gapFrameExt */ \
		(!status) ? (NULL) : (&mtx), any.seq[0] };
	double score = fp_gencore(&sp, &enseq2, &enseq1);
	ck_assert_int_eq((int)score, 7);
	free_scoring_matrix(&mtx);
}END_TEST

START_TEST(test_fp_first_refference_test_trans)
{
	scoring_matrix_t mtx;
	translate_table_t tt;
	int status = read_translate_table(&tt, human40, strlen(human40));
	status = read_scoring_matrix(&mtx, blosum62, strlen(blosum62));

	// >ACHA_ELEEL P09688 electrophorus electricus (electric eel). acetylcholine receptor protein, alpha chain (fragment). 2/94
	char seq1[] = { "SEDETRLVKNLFSGYNKVVRPVNH" };
	size_t len1 = strlen(seq1);
	// >HSBGL2
	char seq2[] = { "ATGTCATACCTCTTATCTCCTCCCACAGCTCCTGGGCAACGTGCTGGTCTGTGTGCTGGCCCATCACTTTGGCAAAGAATTC" };
	char any_symbol = 'x';
	size_t len2 = strlen(seq2);
	sequence_t inseq1 = { 1, (char *)seq1, len1 };
	sequence_t inseq2 = { 2, (char *)seq2, len2 };
	sequence_t inseqrev = { 3, malloc(len2), len2 };
	sequence_t enseq1 = { 1, malloc(len1 + 1), len1 };
	sequence_t enseq2 = { 2, malloc(len2 + 1), len2 };
	sequence_t enseqrev = { 3, malloc(len2 + 1), len2 };
	sequence_t any = { 3, malloc(1 + 1), 1 };
	lal_seq2encodedseq(inseq1, enseq1, lal_encode31);
	lal_seq2encodedseq_trans(inseq2, enseq2, lal_na2indx, &tt);
	lal_seq2encodedseq((sequence_t) { 3, &any_symbol, 1 }, any, lal_encode31);
	mtx.scale = 10.0;
	search_fp_profile_t sp = { -10.5,/*gapOpen */\
		- 0.5, /* gapExt */\
		- 10.5,  /*gapFrame*/ \
		- 13.0, /* matchMax5 */ \
		- 20.0, /* matchMax6 */ \
		- 30.0, /*gapOpen2 */\
		- 23.0, /*gapOpen3 */ \
		- 7.0, /* gapFrameOpen */ \
		- 13.0, /* gapFrameExt */ \
		(!status) ? (NULL) : (&mtx), any.seq[0] };
	double score = fp_gencore(&sp, &enseq2, &enseq1);
	ck_assert_int_eq((int)score, 15);
	lal_reverse(inseq2.seq, inseq2.len, inseqrev.seq, lal_revers31);
	lal_seq2encodedseq_trans(inseqrev, enseqrev, lal_na2indx, &tt);
	score = fp_gencore(&sp, &enseqrev, &enseq1);
	ck_assert_int_eq((int)score, 12);
	free_scoring_matrix(&mtx);

}END_TEST

START_TEST(test_fp_first_refference_test_trans_matrix)
{
	scoring_matrix_t mtx;
	translate_table_t tt;
	int status = read_translate_table(&tt, human40, strlen(human40));
	status = read_scoring_matrix(&mtx, blosum62, strlen(blosum62));

	// >ACHA_ELEEL P09688 electrophorus electricus (electric eel). acetylcholine receptor protein, alpha chain (fragment). 2/94
	char seq1[] = { "SEDETRLVKNLFSGYNKVVRPVNH" };
	size_t len1 = strlen(seq1);
	// >HSBGL2
	char seq2[] = { "ATGTCATACCTCTTATCTCCTCCCACAGCTCCTGGGCAACGTGCTGGTCTGTGTGCTGGCCCATCACTTTGGCAAAGAATTC" };
	char any_symbol = 'x';
	size_t len2 = strlen(seq2);
	sequence_t inseq1 = { 1, (char *)seq1, len1 };
	sequence_t inseq2 = { 2, (char *)seq2, len2 };
	sequence_t inseqrev = { 3, malloc(len2), len2 };
	sequence_t enseq1 = { 1, malloc(len1 + 1), len1 };
	sequence_t enseq2 = { 2, malloc(len2 + 1), len2 };
	sequence_t enseqrev = { 3, malloc(len2 + 1), len2 };
	sequence_t any = { 3, malloc(1 + 1), 1 };
	lal_seq2encodedseq(inseq1, enseq1, lal_encode31);
	lal_seq2encodedseq_trans(inseq2, enseq2, lal_na2indx, &tt);
	lal_seq2encodedseq((sequence_t) { 3, &any_symbol, 1 }, any, lal_encode31);
	mtx.scale = 10.0;
	search_fp_profile_t sp = { -10.5,/*gapOpen */\
		- 0.5, /* gapExt */\
		- 10.5,  /*gapFrame*/ \
		- 13.0, /* matchMax5 */ \
		- 20.0, /* matchMax6 */ \
		- 30.0, /*gapOpen2 */\
		- 23.0, /*gapOpen3 */ \
		- 7.0, /* gapFrameOpen */ \
		- 13.0, /* gapFrameExt */ \
		(!status) ? (NULL) : (&mtx), any.seq[0] };
	double score = fp_gencore_matrix(&sp, &enseq2, &enseq1);
	ck_assert_int_eq((int)score, 15);
	lal_reverse(inseq2.seq, inseq2.len, inseqrev.seq, lal_revers31);
	lal_seq2encodedseq_trans(inseqrev, enseqrev, lal_na2indx, &tt);
	score = fp_gencore_matrix(&sp, &enseqrev, &enseq1);
	ck_assert_int_eq((int)score, 12);
	free_scoring_matrix(&mtx);
}END_TEST

START_TEST(test_fp_first_refference_test_trans_seq)
{
	scoring_matrix_t mtx;
	translate_table_t tt;
	int status = read_translate_table(&tt, human40, strlen(human40));
	status = read_scoring_matrix(&mtx, blosum62, strlen(blosum62));

	// >ACHA_ELEEL P09688 electrophorus electricus (electric eel). acetylcholine receptor protein, alpha chain (fragment). 2/94
	char seq1[] = { "SEDETRLVKNLFSGYNKVVRPVNH" };
	size_t len1 = strlen(seq1);
	// >HSBGL2
	char seq2[] = { "ATGTCATACCTCTTATCTCCTCCCACAGCTCCTGGGCAACGTGCTGGTCTGTGTGCTGGCCCATCACTTTGGCAAAGAATTC" };
	char any_symbol = 'x';
	size_t len2 = strlen(seq2);
	sequence_t inseq1 = { 1, (char *)seq1, len1 };
	sequence_t inseq2 = { 2, (char *)seq2, len2 };
	sequence_t inseqrev = { 3, malloc(len2), len2 };
	sequence_t enseq1 = { 1, malloc(len1 + 1), len1 };
	sequence_t enseq2 = { 2, malloc(len2 + 1), len2 };
	sequence_t enseqrev = { 3, malloc(len2 + 1), len2 };
	sequence_t any = { 3, malloc(1 + 1), 1 };
	lal_seq2encodedseq(inseq1, enseq1, lal_encode31);
	lal_seq2encodedseq_trans(inseq2, enseq2, lal_na2indx, &tt);
	lal_seq2encodedseq((sequence_t) { 3, &any_symbol, 1 }, any, lal_encode31);
	mtx.scale = 10.0;
	search_fp_profile_t sp = { -10.5,/*gapOpen */\
		- 0.5, /* gapExt */\
		- 10.5,  /*gapFrame*/ \
		- 13.0, /* matchMax5 */ \
		- 20.0, /* matchMax6 */ \
		- 30.0, /*gapOpen2 */\
		- 23.0, /*gapOpen3 */ \
		- 7.0, /* gapFrameOpen */ \
		- 13.0, /* gapFrameExt */ \
		(!status) ? (NULL) : (&mtx), any.seq[0] };
	double score = fp_gencore_seq(&sp, &enseq2, &enseq1);
	ck_assert_int_eq((int)score, 15);
	lal_reverse(inseq2.seq, inseq2.len, inseqrev.seq, lal_revers31);
	lal_seq2encodedseq_trans(inseqrev, enseqrev, lal_na2indx, &tt);
	score = fp_gencore_seq(&sp, &enseqrev, &enseq1);
	ck_assert_int_eq((int)score, 12);
	free_scoring_matrix(&mtx);

}END_TEST

START_TEST(test_fp_ACHA_ELEEL_test_model_specific_double)
{
	scoring_matrix_t mtx;
	translate_table_t tt;
	int status = read_translate_table(&tt, human40, strlen(human40));
	status = read_scoring_matrix(&mtx, blosum62, strlen(blosum62));

	// >ACHA_ELEEL P09688 electrophorus electricus (electric eel). acetylcholine receptor protein, alpha chain (fragment). 2/94
	char seq1[] = { "SEDETRLVKNLFSGYNKVVRPVNH" };
	size_t len1 = strlen(seq1);
	// >HSBGL2
	char seq2[] = { "ATGTCATACCTCTTATCTCCTCCCACAGCTCCTGGGCAACGTGCTGGTCTGTGTGCTGGCCCATCACTTTGGCAAAGAATTC" };
	char any_symbol = 'x';
	size_t len2 = strlen(seq2);
	sequence_t inseq1 = { 1, (char *)seq1, len1 };
	sequence_t inseq2 = { 2, (char *)seq2, len2 };
	sequence_t inseqrev = { 3, malloc(len2), len2 };
	sequence_t enseq1 = { 1, malloc(len1 + 1), len1 };
	sequence_t enseq2 = { 2, malloc(len2 + 1), len2 };
	sequence_t enseqrev = { 3, malloc(len2 + 1), len2 };
	sequence_t any = { 3, malloc(1 + 1), 1 };
	lal_seq2encodedseq(inseq1, enseq1, lal_encode31);
	lal_seq2encodedseq_trans(inseq2, enseq2, lal_na2indx, &tt);
	lal_seq2encodedseq((sequence_t) { 3, &any_symbol, 1 }, any, lal_encode31);
	mtx.scale = 10.0;
	search_fp_profile_t sp = { -10.5,/*gapOpen */\
		- 0.5, /* gapExt */\
		- 10.5,  /*gapFrame*/ \
		- 13.0, /* matchMax5 */ \
		- 20.0, /* matchMax6 */ \
		- 30.0, /*gapOpen2 */\
		- 23.0, /*gapOpen3 */ \
		- 7.0, /* gapFrameOpen */ \
		- 13.0, /* gapFrameExt */ \
		(!status) ? (NULL) : (&mtx), any.seq[0], enseq1.len /*very impotant params:  max of sequences length of slice date */};
	search_fp_thr_profile_t *sp_thr = search_fp_thr_init(&sp, 2);
	double score = fp_thr(sp_thr, &enseq2, &enseq1);
	ck_assert_int_eq((int)score, 15);
	lal_reverse(inseq2.seq, inseq2.len, inseqrev.seq, lal_revers31);
	lal_seq2encodedseq_trans(inseqrev, enseqrev, lal_na2indx, &tt);
	score = fp_thr(sp_thr+1, &enseqrev, &enseq1);
	ck_assert_int_eq((int)score, 12);
	search_fp_thr_deinit(sp_thr, 2);
	free_scoring_matrix(&mtx);
}END_TEST

void addFPTC(Suite *s) {
	TCase *tc_core = tcase_create("FP");
	
	tcase_add_test(tc_core, test_fp_ACHA_ELEEL_test_model_specific_double);
	tcase_add_test(tc_core, test_fp_first_refference_test);
	tcase_add_test(tc_core, test_fp_first_refference_test_trans);
	tcase_add_test(tc_core, test_fp_first_refference_test_trans_matrix);
	tcase_add_test(tc_core, test_fp_first_refference_test_trans_seq);
	suite_add_tcase(s, tc_core);
}

