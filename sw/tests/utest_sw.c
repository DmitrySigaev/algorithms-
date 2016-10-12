/*
Copyright (C) 2016

Contact: Dmitry Sigaev <dima.sigaev@gmail.com>
*/

#include "../../utests/tests.h"

#include <stdint.h>
#include <stdio.h>
#include <malloc.h>
#include "../sw.h"
#include "../gc_sw.h"
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
	ck_assert_int_eq((int)score, 2); /* Max score */
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
	ck_assert_int_eq((int)score, 2); /* Max score */
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
	ck_assert_int_eq((int)score, 5); /* Max score */
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
	ck_assert_int_eq((int)score, 2); /* Max score */

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
	ck_assert_int_eq(score, 2); /* Max score */

}END_TEST

START_TEST(test_sw_int_encoded_vtable)
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

START_TEST(test_sw_affine_double)
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



START_TEST(test_sw_affine_double_encoded_vtable)
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

START_TEST(test_sw_affine_double_encoded_vtable_195)
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

START_TEST(test_sw_affine_double_encoded_vtable_195swipe)
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
	ck_assert_int_eq((int)score.dscore, 195); /* Max score */
}END_TEST

START_TEST(test_sw_affine_double_encoded_vtable_195gencore)
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

START_TEST(test_sw_affine_double_encoded_vtable_88)
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

START_TEST(test_sw_affine_double_encoded_vtable_87swipe)
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
	ck_assert_int_eq((int)score.dscore, 87); /* Max score */
}END_TEST


START_TEST(test_sw_affine_double_encoded_vtable_88gencore)
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

START_TEST(test_sw_gaptest1_261_88gencore)
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

START_TEST(test_sw_gaptest1_261_89gotoh)
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

START_TEST(test_sw_gaptest1_261_90swipe)
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
	ck_assert_int_eq((int)score.dscore, 90); /* Max score */ // המכהום בûעü 89

	sequence_t reverse1 = { 3, malloc(len1 + 1), len1 };
	for (size_t i = 0; i < enseq1.len; i++) {
		// computes the reverse complement of the input sequence.
		reverse1.seq[i] = cns[(int)(enseq1.seq[enseq1.len - 1 - i])];
	}
	score = sw_alignment_swipe(&sp, &reverse1, &enseq2);
	ck_assert_int_eq((int)score.dscore, 193); /* Max score */
}END_TEST

// 0x04c9ef2c "CTTCCTACGTTGGGTCATAGTAGTGCGGCGTGGGCAATGCCTACGGAGGGGTGGAGCAACTGGCGCTATCACTTCTACCATCGTCTGCAGCGTACGA" //290
START_TEST(test_sw_gaptest1_290_90_194gencore)
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

START_TEST(test_sw_gaptest1_290_90_194gotoh)
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

START_TEST(test_sw_gaptest1_290_91_194swipe)
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
	ck_assert_int_eq((int)score.dscore, 91); /* Max score */ // המכהום בûעü 89

	sequence_t reverse1 = { 3, malloc(len1 + 1), len1 };
	for (size_t i = 0; i < enseq1.len; i++) {
		// computes the reverse complement of the input sequence.
		reverse1.seq[i] = cns[(int)(enseq1.seq[enseq1.len - 1 - i])];
	}
	score = sw_alignment_swipe(&sp, &reverse1, &enseq2);
	ck_assert_int_eq((int)score.dscore, 194); /* Max score */
}END_TEST

START_TEST(test_sw_gaptest1_290_76_194constant)
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

void addSWTC(Suite *s) {
	TCase *tc_core = tcase_create("SW");
	tcase_add_test(tc_core, test_sw_gaptest1_290_76_194constant); 
	tcase_add_test(tc_core, test_sw_gaptest1_290_91_194swipe); // 90 check
	tcase_add_test(tc_core, test_sw_gaptest1_290_90_194gotoh);
	tcase_add_test(tc_core, test_sw_gaptest1_290_90_194gencore);
	tcase_add_test(tc_core, test_sw_gaptest1_261_90swipe); //89
	tcase_add_test(tc_core, test_sw_gaptest1_261_89gotoh);
	tcase_add_test(tc_core, test_sw_gaptest1_261_88gencore);
	tcase_add_test(tc_core, test_sw_affine_double_encoded_vtable_195gencore);
	tcase_add_test(tc_core, test_sw_affine_double_encoded_vtable_88gencore);
	tcase_add_test(tc_core, test_sw_affine_double_encoded_vtable_87swipe);
	tcase_add_test(tc_core, test_sw_affine_double_encoded_vtable_88);
	tcase_add_test(tc_core, test_sw_affine_double_encoded_vtable_195swipe);
	tcase_add_test(tc_core, test_sw_affine_double_encoded_vtable_195);
	tcase_add_test(tc_core, test_sw_affine_double_encoded_vtable);
	tcase_add_test(tc_core, test_sw_affine_double);
	tcase_add_test(tc_core, test_sw_int_encoded_vtable);
	tcase_add_test(tc_core, test_sw_double_encoded_vtable);
	tcase_add_test(tc_core, test_sw_double_encoded);
	tcase_add_test(tc_core, test_sw_double_symbols);
	tcase_add_test(tc_core, test_sw_double);
	tcase_add_test(tc_core, test_sw_int);

	suite_add_tcase(s, tc_core);
}

