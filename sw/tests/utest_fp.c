/*
Copyright (C) 2017

Contact: Dmitry Sigaev <dima.sigaev@gmail.com>
*/

#include "../../utests/tests.h"

#include <stdint.h>
#include <stdio.h>
#include <malloc.h>
#include "../sw.h"
#include "../gc_sw.h"
#include "../gc_fp.h"
#include "../lal_encoding.h"
#include "../lal_tables.h"


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
	seq2encodedseq(inseq1, enseq1, lal_encode31);
	seq2encodedseq(inseq2, enseq2, lal_encode31);
	seq2encodedseq((sequence_t){ 3, &any_symbol, 1 }, any, lal_encode31);
	mtx.scale = 10.0;
	search_swag_profile_t sp = { -10.5, -0.5, (!status) ? (NULL) : (&mtx), any.seq[0] };
	double score = sw_gencore(&sp, &enseq2, &enseq1);
	ck_assert_int_eq((int)score, 7);
}END_TEST

START_TEST(test_fp_double_symbols)
{
	char seq1[] = { "ABAC" };
	size_t len1 = strlen(seq1);
	char seq2[] = { "ADACCG" };
	size_t len2 = strlen(seq2);
	sequence_t sseq1 = { 1, (char *)seq1, len1 };
	sequence_t sseq2 = { 2, (char *)seq2, len2 };
	search_swcg_profile_t sp = { -1, NULL };
	double score = sw_constant_gap_double(&sp, &sseq1, &sseq2);
	ck_assert_int_eq((int)score, 2); /* Max score */
}END_TEST


START_TEST(test_fp_double_encoded)
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
	ck_assert_int_eq((int)score, 2); /* Max score */
}END_TEST

START_TEST(test_fp_double_encoded_vtable)
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
	ck_assert_int_eq((int)score, 5); /* Max score */
}END_TEST


START_TEST(test_fp_double)
{
	char seq1[] = { 1 , 2 , 1 , 5 };
	size_t len1 = sizeof(seq1) / sizeof(seq1[0]);
	char seq2[] = { 1 , 3 , 1 , 5, 5, 6 };
	size_t len2 = sizeof(seq2) / sizeof(seq2[0]);
	sequence_t sseq1 = { 1, (char *)seq1, len1 };
	sequence_t sseq2 = { 2, (char *)seq2, len2 };
	search_swcg_profile_t sp = { -1, NULL };
	double score = sw_constant_gap_double(&sp, &sseq1, &sseq2);
	ck_assert_int_eq((int)score, 2); /* Max score */

}END_TEST

START_TEST(test_fp_int)
{
	char seq1[] = { 1 , 2 , 1 , 5 };
	size_t len1 = sizeof(seq1) / sizeof(seq1[0]);
	char seq2[] = { 1 , 3 , 1 , 5, 5, 6 };
	size_t len2 = sizeof(seq2) / sizeof(seq2[0]);
	sequence_t sseq1 = { 1, (char *)seq1, len1 };
	sequence_t sseq2 = { 2, (char *)seq2, len2 };
	search_swcg_profile_int_t sp = { -1, NULL };
	int64_t score = sw_constant_gap_int(&sp, &sseq1, &sseq2);
	ck_assert_int_eq(score, 2); /* Max score */

}END_TEST

START_TEST(test_fp_int_encoded_vtable)
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
	search_swcg_profile_int_t sp = { -1, (!status) ? (NULL) : (&mtx) };
	int64_t score = sw_constant_gap_int(&sp, &enseq1, &enseq2);
	ck_assert_int_eq((int)score, 5); /* Max score */
}END_TEST

START_TEST(test_fp_affine_double)
{
	char seq1[] = { 1 , 2 , 1 , 5 };
	size_t len1 = sizeof(seq1) / sizeof(seq1[0]);
	char seq2[] = { 1 , 3 , 1 , 5, 5, 6 };
	size_t len2 = sizeof(seq2) / sizeof(seq2[0]);
	sequence_t sseq1 = { 1, (char *)seq1, len1 };
	sequence_t sseq2 = { 2, (char *)seq2, len2 };
	search_swag_profile_t sp = { 0, -1, NULL };
	double score = sw_affine_gap(&sp, &sseq1, &sseq2);
	ck_assert_int_eq((int)score, 3); /* Max score */
}END_TEST



START_TEST(test_fp_affine_double_encoded_vtable)
{
	scoring_matrix_t mtx;
	int status = read_scoring_matrix(&mtx, gaptest1, strlen(gaptest1));
	char seq1[] = { "CCC" };
	size_t len1 = strlen(seq1);
	char seq2[] = { "ACACCTT" };
	size_t len2 = strlen(seq2);
	sequence_t inseq1 = { 1, (char *)seq1, len1 };
	sequence_t inseq2 = { 2, (char *)seq2, len2 };
	sequence_t enseq1 = { 1, malloc(len1 + 1), len1 };
	sequence_t enseq2 = { 2, malloc(len2 + 1), len2 };
	seq2encodedseq(inseq1, enseq1, lal_encode31);
	seq2encodedseq(inseq2, enseq2, lal_encode31);
	search_swag_profile_t sp = { -1, -1, (!status) ? (NULL) : (&mtx) };
	double score = sw_affine_gap(&sp, &enseq1, &enseq2);
	ck_assert_int_eq((int)score, 5); /* Max score */
}END_TEST

START_TEST(test_fp_affine_double_encoded_vtable_195)
{
	scoring_matrix_t mtx;
	int status = read_scoring_matrix(&mtx, gaptest1, strlen(gaptest1));
//	char seq1[] = { "TCGTACGCTGCAGACGATGGTAGAAGTGATAGCGCCAGTTGCTCCACCCCTCCGTAGGCATTGCACGCCGCACTACTATGACCCAACGTAGGAAGTTGG" };
	char seq1[] = { "TCGTACGCTGCAGACGATGGTAGAAGTGATAGCGCCAGTTGCTCCACCCCTCCGTAGGCATTGCACGCCGCACTACTATGACCCAACGTAGGAAGTTGG" };
	size_t len1 = strlen(seq1);
	char seq2[] = { "TCGTACGCTGCAGACGATGGTAGAAGTGATAGCGCCAGTTGCTCCACCCCTCCGTAGGCATTGCCCACGCCGCACTACTATGACCCAACGTAGGAAGTTG" };
	size_t len2 = strlen(seq2);
	sequence_t inseq1 = { 1, (char *)seq1, len1 };
	sequence_t inseq2 = { 2, (char *)seq2, len2 };
	sequence_t enseq1 = { 1, malloc(len1 + 1), len1 };
	sequence_t enseq2 = { 2, malloc(len2 + 1), len2 };
	seq2encodedseq(inseq1, enseq1, lal_encode31);
	seq2encodedseq(inseq2, enseq2, lal_encode31);
	search_swag_profile_t sp = { -1, 0, (!status) ? (NULL) : (&mtx) };
	double score = sw_affine_gap(&sp, &enseq1, &enseq2);
	ck_assert_int_eq((int)score, 195); /* Max score */
}END_TEST

START_TEST(test_fp_affine_double_encoded_vtable_195swipe)
{
	scoring_matrix_t mtx;
	int status = read_scoring_matrix(&mtx, gaptest1, strlen(gaptest1));
	//	char seq1[] = { "TCGTACGCTGCAGACGATGGTAGAAGTGATAGCGCCAGTTGCTCCACCCCTCCGTAGGCATTGCACGCCGCACTACTATGACCCAACGTAGGAAGTTGG" };
	char seq1[] = { "TCGTACGCTGCAGACGATGGTAGAAGTGATAGCGCCAGTTGCTCCACCCCTCCGTAGGCATTGCACGCCGCACTACTATGACCCAACGTAGGAAGTTGG" };
	size_t len1 = strlen(seq1);
	char seq2[] = { "TCGTACGCTGCAGACGATGGTAGAAGTGATAGCGCCAGTTGCTCCACCCCTCCGTAGGCATTGCCCACGCCGCACTACTATGACCCAACGTAGGAAGTTG" };
	size_t len2 = strlen(seq2);
	sequence_t inseq1 = { 1, (char *)seq1, len1 };
	sequence_t inseq2 = { 2, (char *)seq2, len2 };
	sequence_t enseq1 = { 1, malloc(len1 + 1), len1 };
	sequence_t enseq2 = { 2, malloc(len2 + 1), len2 };
	seq2encodedseq(inseq1, enseq1, lal_encode31);
	seq2encodedseq(inseq2, enseq2, lal_encode31);
	search_swag_profile_t sp = { -1, 0, (!status) ? (NULL) : (&mtx) };
	region_t score = sw_alignment_swipe(&sp, &enseq1, &enseq2);
	ck_assert_int_eq((int)score.fdscore, 195); /* Max forward score */
	ck_assert_int_eq((int)score.bdscore, 195); /* Max backward score */
}END_TEST

START_TEST(test_fp_affine_double_encoded_vtable_195gencore)
{
	scoring_matrix_t mtx;
	int status = read_scoring_matrix(&mtx, gaptest1, strlen(gaptest1));
	//	char seq1[] = { "TCGTACGCTGCAGACGATGGTAGAAGTGATAGCGCCAGTTGCTCCACCCCTCCGTAGGCATTGCACGCCGCACTACTATGACCCAACGTAGGAAGTTGG" };
	char seq1[] = { "TCGTACGCTGCAGACGATGGTAGAAGTGATAGCGCCAGTTGCTCCACCCCTCCGTAGGCATTGCACGCCGCACTACTATGACCCAACGTAGGAAGTTGG" };

	size_t len1 = strlen(seq1);
	char seq2[] = { "TCGTACGCTGCAGACGATGGTAGAAGTGATAGCGCCAGTTGCTCCACCCCTCCGTAGGCATTGCCCACGCCGCACTACTATGACCCAACGTAGGAAGTTG" };
//	0x01c79c10      "TCGTACGCTGCAGACGATGGTAGAAGTGATAGCGCCAGTTGCTCCACCCCTCCGTAGGCATTGCCCACGCCGCACTACTATGACCCAACGTAGGAAGTTG"

	size_t len2 = strlen(seq2);
	sequence_t inseq1 = { 1, (char *)seq1, len1 };
	sequence_t inseq2 = { 2, (char *)seq2, len2 };
	sequence_t enseq1 = { 1, malloc(len1 + 1), len1 };
	sequence_t enseq2 = { 2, malloc(len2 + 1), len2 };
	seq2encodedseq(inseq1, enseq1, lal_encode31);
	seq2encodedseq(inseq2, enseq2, lal_encode31);
	search_swag_profile_t sp = { -1, 0, (!status) ? (NULL) : (&mtx) };
	double score = sw_gencore(&sp, &enseq1, &enseq2);
	ck_assert_int_eq((int)score, 195); /* Max score */
}END_TEST

START_TEST(test_fp_affine_double_encoded_vtable_88)
{
	scoring_matrix_t mtx;
	int status = read_scoring_matrix(&mtx, identity_nuc, strlen(identity_nuc));
	char seq1[] = { "TCGTACGCTGCAACGATGGTAGAAGTGATAGCGCCAGTTGCTCCACCCCTCCGTAGGCATTGCCCACGCCGCACTACTATGACCCAACGTAGGAAGTTG" };
	size_t len1 = strlen(seq1);
	char seq2[] = { "TCGTACGCTGCAGACGATGGTAGAAGTGATAGCGCCAGTTGCTCCACCCCTCCGTAGGCATTGCCCACGCCGCACTACTATGACCCAACGTAGGAAGTTG" };
	size_t len2 = strlen(seq2);
	sequence_t inseq1 = { 1, (char *)seq1, len1 };
	sequence_t inseq2 = { 2, (char *)seq2, len2 };
	sequence_t enseq1 = { 1, malloc(len1 + 1), len1 };
	sequence_t enseq2 = { 2, malloc(len2 + 1), len2 };
	seq2encodedseq(inseq1, enseq1, lal_encode31);
	seq2encodedseq(inseq2, enseq2, lal_encode31);
	search_swag_profile_t sp = { -11, -1, (!status) ? (NULL) : (&mtx) };
	double score = sw_affine_gap(&sp, &enseq1, &enseq2);
	ck_assert_int_eq((int)score, 88); /* Max score */
}END_TEST

START_TEST(test_fp_affine_double_encoded_vtable_87swipe)
{
	scoring_matrix_t mtx;
	int status = read_scoring_matrix(&mtx, identity_nuc, strlen(identity_nuc));
	char seq1[] = { "TCGTACGCTGCAACGATGGTAGAAGTGATAGCGCCAGTTGCTCCACCCCTCCGTAGGCATTGCCCACGCCGCACTACTATGACCCAACGTAGGAAGTTG" };
	size_t len1 = strlen(seq1);
	char seq2[] = { "TCGTACGCTGCAGACGATGGTAGAAGTGATAGCGCCAGTTGCTCCACCCCTCCGTAGGCATTGCCCACGCCGCACTACTATGACCCAACGTAGGAAGTTG" };
	size_t len2 = strlen(seq2);
	sequence_t inseq1 = { 1, (char *)seq1, len1 };
	sequence_t inseq2 = { 2, (char *)seq2, len2 };
	sequence_t enseq1 = { 1, malloc(len1 + 1), len1 };
	sequence_t enseq2 = { 2, malloc(len2 + 1), len2 };
	seq2encodedseq(inseq1, enseq1, lal_encode31);
	seq2encodedseq(inseq2, enseq2, lal_encode31);
	search_swag_profile_t sp = { -11, -1, (!status) ? (NULL) : (&mtx) };
	region_t score = sw_alignment_swipe(&sp, &enseq1, &enseq2);
	ck_assert_int_eq((int)score.fdscore, 87); /* Max forward score */
	ck_assert_int_eq((int)score.bdscore, 87); /* Max backward score */
}END_TEST


START_TEST(test_fp_affine_double_encoded_vtable_88gencore)
{
	scoring_matrix_t mtx;
	int status = read_scoring_matrix(&mtx, identity_nuc, strlen(identity_nuc));
	char seq1[] = { "TCGTACGCTGCAACGATGGTAGAAGTGATAGCGCCAGTTGCTCCACCCCTCCGTAGGCATTGCCCACGCCGCACTACTATGACCCAACGTAGGAAGTTG" };
	size_t len1 = strlen(seq1);
	char seq2[] = { "TCGTACGCTGCAGACGATGGTAGAAGTGATAGCGCCAGTTGCTCCACCCCTCCGTAGGCATTGCCCACGCCGCACTACTATGACCCAACGTAGGAAGTTG" };
	size_t len2 = strlen(seq2);
	sequence_t inseq1 = { 1, (char *)seq1, len1 };
	sequence_t inseq2 = { 2, (char *)seq2, len2 };
	sequence_t enseq1 = { 1, malloc(len1 + 1), len1 };
	sequence_t enseq2 = { 2, malloc(len2 + 1), len2 };
	seq2encodedseq(inseq1, enseq1, lal_encode31);
	seq2encodedseq(inseq2, enseq2, lal_encode31);
	search_swag_profile_t sp = { -11, -1, (!status) ? (NULL) : (&mtx) };
	double score = sw_gencore(&sp, &enseq1, &enseq2);
	ck_assert_int_eq((int)score, 88); /* Max score */
}END_TEST

//0x040fa27c "CAACTTCCTACGTTGGGTCATAGTAGTGCGTGGGCAATGCCTACGGAGGGGTGGAGCAACTGGCGCTATCACTTCTACCATCGTCTGCAGCGTACGA"   //_261

START_TEST(test_fp_gaptest1_261_88gencore)
{
	scoring_matrix_t mtx;
	int status = read_scoring_matrix(&mtx, gaptest1, strlen(gaptest1));
	//0x040fa27c "CAACTTCCTACGTTGGGTCATAGTAGTGCGTGGGCAATGCCTACGGAGGGGTGGAGCAACTGGCGCTATCACTTCTACCATCGTCTGCAGCGTACGA"   //_261
	char seq1[] = { "CAACTTCCTACGTTGGGTCATAGTAGTGCGTGGGCAATGCCTACGGAGGGGTGGAGCAACTGGCGCTATCACTTCTACCATCGTCTGCAGCGTACGA" };

	size_t len1 = strlen(seq1);
//	char seq2[] = { "TCGTACGCTGCAGACGATGGTAGAAGTGATAGCGCCAGTTGCTCCACCCCTCCGTAGGCATTGCCCACGCCGCACTACTATGACCCAACGTAGGAAGTTG" };
	char seq2[] = { "tcgtacgctgcagacgatggtagaagtgatagcgccagttgctccacccctccgtaggcattgcccacgccgcactactatgacccaacgtaggaagttg" };
	//	0x01c79c10      "TCGTACGCTGCAGACGATGGTAGAAGTGATAGCGCCAGTTGCTCCACCCCTCCGTAGGCATTGCCCACGCCGCACTACTATGACCCAACGTAGGAAGTTG"

	size_t len2 = strlen(seq2);
	sequence_t inseq1 = { 1, (char *)seq1, len1 };
	sequence_t inseq2 = { 2, (char *)seq2, len2 };
	sequence_t enseq1 = { 1, malloc(len1 + 1), len1 };
	sequence_t enseq2 = { 2, malloc(len2 + 1), len2 };
	seq2encodedseq(inseq1, enseq1, lal_encode31);
	seq2encodedseq(inseq2, enseq2, lal_encode31);
	search_swag_profile_t sp = { -1, 0, (!status) ? (NULL) : (&mtx) };
	double score = sw_gencore(&sp, &enseq1, &enseq2);
	ck_assert_int_eq((int)score, 88);

	sequence_t reverse1 = { 3, malloc(len1 + 1), len1 };
	for (size_t i = 0; i < enseq1.len; i++) {
		// computes the reverse complement of the input sequence.
		reverse1.seq[i] = cns[(int)(enseq1.seq[enseq1.len - 1 - i])];
	}
	score = sw_gencore(&sp, &reverse1, &enseq2);
	ck_assert_int_eq((int)score, 193);
}END_TEST

START_TEST(test_fp_gaptest1_261_89gotoh)
{
	scoring_matrix_t mtx;
	int status = read_scoring_matrix(&mtx, gaptest1, strlen(gaptest1));
	//_261
	//0x040fa27c "CAACTTCCTACGTTGGGTCATAGTAGTGCGTGGGCAATGCCTACGGAGGGGTGGAGCAACTGGCGCTATCACTTCTACCATCGTCTGCAGCGTACGA"
	char seq1[] = { "CAACTTCCTACGTTGGGTCATAGTAGTGCGTGGGCAATGCCTACGGAGGGGTGGAGCAACTGGCGCTATCACTTCTACCATCGTCTGCAGCGTACGA" };

	size_t len1 = strlen(seq1);
	//	char seq2[] = { "TCGTACGCTGCAGACGATGGTAGAAGTGATAGCGCCAGTTGCTCCACCCCTCCGTAGGCATTGCCCACGCCGCACTACTATGACCCAACGTAGGAAGTTG" };
	char seq2[] = { "tcgtacgctgcagacgatggtagaagtgatagcgccagttgctccacccctccgtaggcattgcccacgccgcactactatgacccaacgtaggaagttg" };
	//	0x01c79c10      "TCGTACGCTGCAGACGATGGTAGAAGTGATAGCGCCAGTTGCTCCACCCCTCCGTAGGCATTGCCCACGCCGCACTACTATGACCCAACGTAGGAAGTTG"

	size_t len2 = strlen(seq2);
	sequence_t inseq1 = { 1, (char *)seq1, len1 };
	sequence_t inseq2 = { 2, (char *)seq2, len2 };
	sequence_t enseq1 = { 1, malloc(len1 + 1), len1 };
	sequence_t enseq2 = { 2, malloc(len2 + 1), len2 };
	seq2encodedseq(inseq1, enseq1, lal_encode31);
	seq2encodedseq(inseq2, enseq2, lal_encode31);
	search_swag_profile_t sp = { -1, 0, (!status) ? (NULL) : (&mtx) };
	double score = sw_affine_gap(&sp, &enseq1, &enseq2);
	ck_assert_int_eq((int)score, 89);

	sequence_t reverse1 = { 3, malloc(len1 + 1), len1 };
	for (size_t i = 0; i < enseq1.len; i++) {
		// computes the reverse complement of the input sequence.
		reverse1.seq[i] = cns[(int)(enseq1.seq[enseq1.len - 1 - i])];
	}
	score = sw_affine_gap(&sp, &reverse1, &enseq2);
	ck_assert_int_eq((int)score, 193);
}END_TEST

START_TEST(test_fp_gaptest1_261_90swipe)
{
	scoring_matrix_t mtx;
	int status = read_scoring_matrix(&mtx, gaptest1, strlen(gaptest1));
	//_261
	//0x040fa27c "CAACTTCCTACGTTGGGTCATAGTAGTGCGTGGGCAATGCCTACGGAGGGGTGGAGCAACTGGCGCTATCACTTCTACCATCGTCTGCAGCGTACGA"
	char seq1[] = { "CAACTTCCTACGTTGGGTCATAGTAGTGCGTGGGCAATGCCTACGGAGGGGTGGAGCAACTGGCGCTATCACTTCTACCATCGTCTGCAGCGTACGA" };

	size_t len1 = strlen(seq1);
	//	char seq2[] = { "TCGTACGCTGCAGACGATGGTAGAAGTGATAGCGCCAGTTGCTCCACCCCTCCGTAGGCATTGCCCACGCCGCACTACTATGACCCAACGTAGGAAGTTG" };
	char seq2[] = { "tcgtacgctgcagacgatggtagaagtgatagcgccagttgctccacccctccgtaggcattgcccacgccgcactactatgacccaacgtaggaagttg" };
	//	0x01c79c10      "TCGTACGCTGCAGACGATGGTAGAAGTGATAGCGCCAGTTGCTCCACCCCTCCGTAGGCATTGCCCACGCCGCACTACTATGACCCAACGTAGGAAGTTG"

	size_t len2 = strlen(seq2);
	sequence_t inseq1 = { 1, (char *)seq1, len1 };
	sequence_t inseq2 = { 2, (char *)seq2, len2 };
	sequence_t enseq1 = { 1, malloc(len1 + 1), len1 };
	sequence_t enseq2 = { 2, malloc(len2 + 1), len2 };
	seq2encodedseq(inseq1, enseq1, lal_encode31);
	seq2encodedseq(inseq2, enseq2, lal_encode31);
	search_swag_profile_t sp = { -1, 0, (!status) ? (NULL) : (&mtx) };
	region_t score = sw_alignment_swipe(&sp, &enseq1, &enseq2);
	ck_assert_int_eq((int)score.fdscore, 89); /* Max forward score */ 
	ck_assert_int_eq((int)score.bdscore, 90); /* Max backward score */ // המכהום בûעü 89??

	sequence_t reverse1 = { 3, malloc(len1 + 1), len1 };
	for (size_t i = 0; i < enseq1.len; i++) {
		// computes the reverse complement of the input sequence.
		reverse1.seq[i] = cns[(int)(enseq1.seq[enseq1.len - 1 - i])];
	}
	score = sw_alignment_swipe(&sp, &reverse1, &enseq2);
	ck_assert_int_eq((int)score.fdscore, 193); /* Max forward score */
	ck_assert_int_eq((int)score.bdscore, 193); /* Max backward score */
}END_TEST

// 0x04c9ef2c "CTTCCTACGTTGGGTCATAGTAGTGCGGCGTGGGCAATGCCTACGGAGGGGTGGAGCAACTGGCGCTATCACTTCTACCATCGTCTGCAGCGTACGA" //290
START_TEST(test_fp_gaptest1_290_90_194gencore)
{
	scoring_matrix_t mtx;
	int status = read_scoring_matrix(&mtx, gaptest1, strlen(gaptest1));
	char seq1[] = { "CTTCCTACGTTGGGTCATAGTAGTGCGGCGTGGGCAATGCCTACGGAGGGGTGGAGCAACTGGCGCTATCACTTCTACCATCGTCTGCAGCGTACGA" };
	size_t len1 = strlen(seq1);
	char seq2[] = { "tcgtacgctgcagacgatggtagaagtgatagcgccagttgctccacccctccgtaggcattgcccacgccgcactactatgacccaacgtaggaagttg" };
	size_t len2 = strlen(seq2);
	sequence_t inseq1 = { 1, (char *)seq1, len1 };
	sequence_t inseq2 = { 2, (char *)seq2, len2 };
	sequence_t enseq1 = { 1, malloc(len1 + 1), len1 };
	sequence_t enseq2 = { 2, malloc(len2 + 1), len2 };
	seq2encodedseq(inseq1, enseq1, lal_encode31);
	seq2encodedseq(inseq2, enseq2, lal_encode31);
	search_swag_profile_t sp = { -1, 0, (!status) ? (NULL) : (&mtx) };
	double score = sw_gencore(&sp, &enseq1, &enseq2);
	ck_assert_int_eq((int)score, 90);

	sequence_t reverse1 = { 3, malloc(len1 + 1), len1 };
	for (size_t i = 0; i < enseq1.len; i++) {
		// computes the reverse complement of the input sequence.
		reverse1.seq[i] = cns[(int)(enseq1.seq[enseq1.len - 1 - i])];
	}
	score = sw_gencore(&sp, &reverse1, &enseq2);
	ck_assert_int_eq((int)score, 194);
}END_TEST

START_TEST(test_fp_gaptest1_290_90_194gotoh)
{
	scoring_matrix_t mtx;
	int status = read_scoring_matrix(&mtx, gaptest1, strlen(gaptest1));
	char seq1[] = { "CTTCCTACGTTGGGTCATAGTAGTGCGGCGTGGGCAATGCCTACGGAGGGGTGGAGCAACTGGCGCTATCACTTCTACCATCGTCTGCAGCGTACGA" };
	size_t len1 = strlen(seq1);
	char seq2[] = { "tcgtacgctgcagacgatggtagaagtgatagcgccagttgctccacccctccgtaggcattgcccacgccgcactactatgacccaacgtaggaagttg" };
	size_t len2 = strlen(seq2);
	sequence_t inseq1 = { 1, (char *)seq1, len1 };
	sequence_t inseq2 = { 2, (char *)seq2, len2 };
	sequence_t enseq1 = { 1, malloc(len1 + 1), len1 };
	sequence_t enseq2 = { 2, malloc(len2 + 1), len2 };
	seq2encodedseq(inseq1, enseq1, lal_encode31);
	seq2encodedseq(inseq2, enseq2, lal_encode31);
	search_swag_profile_t sp = { -1, 0, (!status) ? (NULL) : (&mtx) };
	double score = sw_affine_gap(&sp, &enseq1, &enseq2);
	ck_assert_int_eq((int)score, 90);

	sequence_t reverse1 = { 3, malloc(len1 + 1), len1 };
	for (size_t i = 0; i < enseq1.len; i++) {
		// computes the reverse complement of the input sequence.
		reverse1.seq[i] = cns[(int)(enseq1.seq[enseq1.len - 1 - i])];
	}
	score = sw_affine_gap(&sp, &reverse1, &enseq2);
	ck_assert_int_eq((int)score, 194);
}END_TEST

START_TEST(test_fp_gaptest1_290_91_194swipe)
{
	scoring_matrix_t mtx;
	int status = read_scoring_matrix(&mtx, gaptest1, strlen(gaptest1));
	char seq1[] = { "CTTCCTACGTTGGGTCATAGTAGTGCGGCGTGGGCAATGCCTACGGAGGGGTGGAGCAACTGGCGCTATCACTTCTACCATCGTCTGCAGCGTACGA" };
	size_t len1 = strlen(seq1);
	char seq2[] = { "tcgtacgctgcagacgatggtagaagtgatagcgccagttgctccacccctccgtaggcattgcccacgccgcactactatgacccaacgtaggaagttg" };
	//char seq2[]={ "TCGTACGCTGCAGACGATGGTAGAAGTGATAGCGCCAGTTGCTCCACCCCTCCGTAGGCATTGCCCACGCCGCACTACTATGACCCAACGTAGGAAGTTG" };
	size_t len2 = strlen(seq2);
	sequence_t inseq1 = { 1, (char *)seq1, len1 };
	sequence_t inseq2 = { 2, (char *)seq2, len2 };
	sequence_t enseq1 = { 1, malloc(len1 + 1), len1 };
	sequence_t enseq2 = { 2, malloc(len2 + 1), len2 };
	seq2encodedseq(inseq1, enseq1, lal_encode31);
	seq2encodedseq(inseq2, enseq2, lal_encode31);
	search_swag_profile_t sp = { -1, 0, (!status) ? (NULL) : (&mtx) };
	region_t score = sw_alignment_swipe(&sp, &enseq1, &enseq2);
	ck_assert_int_eq((int)score.fdscore, 90); /* Max forward score */ // ok
	ck_assert_int_eq((int)score.bdscore, 91); /* Max backward score */ // should it be equil to 90 ??

	sequence_t reverse1 = { 3, malloc(len1 + 1), len1 };
	for (size_t i = 0; i < enseq1.len; i++) {
		// computes the reverse complement of the input sequence.
		reverse1.seq[i] = cns[(int)(enseq1.seq[enseq1.len - 1 - i])];
	}
	score = sw_alignment_swipe(&sp, &reverse1, &enseq2);
	ck_assert_int_eq((int)score.fdscore, 194); /* Max forward score */
	ck_assert_int_eq((int)score.bdscore, 194); /* Max backward score */
}END_TEST

START_TEST(test_fp_gaptest1_290_89_193swdirection)
{
	scoring_matrix_t mtx;
	int status = read_scoring_matrix(&mtx, gaptest1, strlen(gaptest1));
	char seq1[] = { "CTTCCTACGTTGGGTCATAGTAGTGCGGCGTGGGCAATGCCTACGGAGGGGTGGAGCAACTGGCGCTATCACTTCTACCATCGTCTGCAGCGTACGA" };
	size_t len1 = strlen(seq1);
	char seq2[] = { "tcgtacgctgcagacgatggtagaagtgatagcgccagttgctccacccctccgtaggcattgcccacgccgcactactatgacccaacgtaggaagttg" };
	//char seq2[]={ "TCGTACGCTGCAGACGATGGTAGAAGTGATAGCGCCAGTTGCTCCACCCCTCCGTAGGCATTGCCCACGCCGCACTACTATGACCCAACGTAGGAAGTTG" };
	size_t len2 = strlen(seq2);
	sequence_t inseq1 = { 1, (char *)seq1, len1 };
	sequence_t inseq2 = { 2, (char *)seq2, len2 };
	sequence_t enseq1 = { 1, malloc(len1 + 1), len1 };
	sequence_t enseq2 = { 2, malloc(len2 + 1), len2 };
	seq2encodedseq(inseq1, enseq1, lal_encode31);
	seq2encodedseq(inseq2, enseq2, lal_encode31);
	search_swag_profile_t sp = { -1, 0, (!status) ? (NULL) : (&mtx) };
	score_matrix_t sd = sw_directions(&sp, &enseq1, &enseq2);
	element_t score = find_max(&sd.score);
	ck_assert_int_eq((int)score.d, 89); /* Max score */ // ok
	free_matrix(&sd.score);
	free_matrix(&sd.directions);
	sequence_t reverse1 = { 3, malloc(len1 + 1), len1 };
	for (size_t i = 0; i < enseq1.len; i++) {
		// computes the reverse complement of the input sequence.
		reverse1.seq[i] = cns[(int)(enseq1.seq[enseq1.len - 1 - i])];
	}
	sd = sw_directions(&sp, &reverse1, &enseq2);
	score = find_max(&sd.score);
	ck_assert_int_eq((int)score.d, 193); /* Max  score */
	free_matrix(&sd.score);
	free_matrix(&sd.directions);
}END_TEST

START_TEST(test_fp_gaptest1_290_76_194constant)
{
	scoring_matrix_t mtx;
	int status = read_scoring_matrix(&mtx, gaptest1, strlen(gaptest1));
	char seq1[] = { "CTTCCTACGTTGGGTCATAGTAGTGCGGCGTGGGCAATGCCTACGGAGGGGTGGAGCAACTGGCGCTATCACTTCTACCATCGTCTGCAGCGTACGA" };
	size_t len1 = strlen(seq1);
	char seq2[] = { "tcgtacgctgcagacgatggtagaagtgatagcgccagttgctccacccctccgtaggcattgcccacgccgcactactatgacccaacgtaggaagttg" };
	//char seq2[]={ "TCGTACGCTGCAGACGATGGTAGAAGTGATAGCGCCAGTTGCTCCACCCCTCCGTAGGCATTGCCCACGCCGCACTACTATGACCCAACGTAGGAAGTTG" };
	size_t len2 = strlen(seq2);
	sequence_t inseq1 = { 1, (char *)seq1, len1 };
	sequence_t inseq2 = { 2, (char *)seq2, len2 };
	sequence_t enseq1 = { 1, malloc(len1 + 1), len1 };
	sequence_t enseq2 = { 2, malloc(len2 + 1), len2 };
	seq2encodedseq(inseq1, enseq1, lal_encode31);
	seq2encodedseq(inseq2, enseq2, lal_encode31);
	search_swcg_profile_t sp = { -1, (!status) ? (NULL) : (&mtx) };
	double score = sw_constant_gap_double(&sp, &enseq1, &enseq2);
	ck_assert_int_eq((int)score, 76); /* Max score */ // המכהום בûעü 89

	sequence_t reverse1 = { 3, malloc(len1 + 1), len1 };
	for (size_t i = 0; i < enseq1.len; i++) {
		// computes the reverse complement of the input sequence.
		reverse1.seq[i] = cns[(int)(enseq1.seq[enseq1.len - 1 - i])];
	}
	score = sw_constant_gap_double(&sp, &reverse1, &enseq2);
	ck_assert_int_eq((int)score, 194); /* Max score */
}END_TEST

START_TEST(test_fp_gaptest1_290_89_193swdirection2)
{
	scoring_matrix_t mtx;
	int status = read_scoring_matrix(&mtx, gaptest1, strlen(gaptest1));
	char seq1[] = { "CTTCCTACGTTGGGTCATAGTAGTGCGGCGTGGGCAATGCCTACGGAGGGGTGGAGCAACTGGCGCTATCACTTCTACCATCGTCTGCAGCGTACGA" };
	size_t len1 = strlen(seq1);
	char seq2[] = { "tcgtacgctgcagacgatggtagaagtgatagcgccagttgctccacccctccgtaggcattgcccacgccgcactactatgacccaacgtaggaagttg" };
	//char seq2[]={ "TCGTACGCTGCAGACGATGGTAGAAGTGATAGCGCCAGTTGCTCCACCCCTCCGTAGGCATTGCCCACGCCGCACTACTATGACCCAACGTAGGAAGTTG" };
	size_t len2 = strlen(seq2);
	sequence_t inseq1 = { 1, (char *)seq1, len1 };
	sequence_t inseq2 = { 2, (char *)seq2, len2 };
	sequence_t enseq1 = { 1, malloc(len1 + 1), len1 };
	sequence_t enseq2 = { 2, malloc(len2 + 1), len2 };
	seq2encodedseq(inseq1, enseq1, lal_encode31);
	seq2encodedseq(inseq2, enseq2, lal_encode31);
	search_swag_profile_t sp = { -1, 0, (!status) ? (NULL) : (&mtx) };
	score_matrix_t sd = sw_directions(&sp, &enseq1, &enseq2);
	element_t score = find_max(&sd.score);
	ck_assert_int_eq((int)score.d, 89); /* Max score */ // ok
	free_matrix(&sd.score);
	free_matrix(&sd.directions);
	sequence_t reverse1 = { 3, malloc(len1 + 1), len1 };
	for (size_t i = 0; i < enseq1.len; i++) {
		// computes the reverse complement of the input sequence.
		reverse1.seq[i] = cns[(int)(enseq1.seq[enseq1.len - 1 - i])];
	}
	sd = sw_directions(&sp, &reverse1, &enseq2);
	score = find_max(&sd.score);
	ck_assert_int_eq((int)score.d, 193); /* Max  score */
	free_matrix(&sd.score);
	free_matrix(&sd.directions);
}END_TEST


void addFPTC(Suite *s) {
	TCase *tc_core = tcase_create("FP");
	
	tcase_add_test(tc_core, test_fp_first_refference_test);

	suite_add_tcase(s, tc_core);
}

