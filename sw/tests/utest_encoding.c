/*
Copyright (C) 2016

Contact: Dmitry Sigaev <dima.sigaev@gmail.com>
*/

#include "../../utests/tests.h"

#include <stdint.h>
#include <stdio.h>
#include <malloc.h>
#include "../lal_typedefs.h"
#include "../lal_encoding.h"
#include "../lal_tables.h"


// 0x04c9ef2c "CTTCCTACGTTGGGTCATAGTAGTGCGGCGTGGGCAATGCCTACGGAGGGGTGGAGCAACTGGCGCTATCACTTCTACCATCGTCTGCAGCGTACGA" //290
START_TEST(test_encoding_290)
{
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


	sequence_t reverse1 = { 3, malloc(len1 + 1), len1 };
	for (size_t i = 0; i < enseq1.len; i++) {
		// computes the reverse complement of the input sequence.
		reverse1.seq[i] = cns[(int)(enseq1.seq[enseq1.len - 1 - i])];
	}
	char rev1[] = { "CTTCCTACGTTGGGTCATAGTAGTGCGGCGTGGGCAATGCCTACGGAGGGGTGGAGCAACTGGCGCTATCACTTCTACCATCGTCTGCAGCGTACGA" };
	size_t revlen1 = strlen(rev1);
	lal_reverse(seq1, len1, rev1, lal_revers31_s);
	//0x0104(rev1)  "TCGTACGCTGCAGACGATGGTAGAAGTGATAGCGCCAGTTGCTCCACCCCTCCGTAGGCATTGCCCACGCCGCACTACTATGACCCAACGTAGGAAG"
	char rev2[] = { "TCGTACGCTGCAGACGATGGTAGAAGTGATAGCGCCAGTTGCTCCACCCCTCCGTAGGCATTGCCCACGCCGCACTACTATGACCCAACGTAGGAAG" };
	
	for (size_t i = 0; i < reverse1.len; i++) {
		ck_assert_int_eq((int)rev1[i], (int)rev2[i]);
	}

	sequence_t inrev1 = { 3, (char *)rev1, len1 };
	sequence_t enrev1 = { 3, malloc(len1 + 1), len1 };
	seq2encodedseq(inrev1, enrev1, lal_encode31);

	for (size_t i = 0; i < reverse1.len; i++) {
		ck_assert_int_eq((int)reverse1.seq[i], (int)enrev1.seq[i]); 
	}
}END_TEST


// 0x04c9ef2c "CAACTTCCTACGTTGGGTCATAGTAGTGCGTGGGCAATGCCTACGGAGGGGTGGAGCAACTGGCGCTATCACTTCTACCATCGTCTGCAGCGTACGA" //261
START_TEST(test_encoding_261)
{
	char seq1[] = { "CAACTTCCTACGTTGGGTCATAGTAGTGCGTGGGCAATGCCTACGGAGGGGTGGAGCAACTGGCGCTATCACTTCTACCATCGTCTGCAGCGTACGA" };
	size_t len1 = strlen(seq1);
	char seq2[] = { "tcgtacgctgcagacgatggtagaagtgatagcgccagttgctccacccctccgtaggcattgcccacgccgcactactatgacccaacgtaggaagttg" };
	size_t len2 = strlen(seq2);
	sequence_t inseq1 = { 1, (char *)seq1, len1 };
	sequence_t inseq2 = { 2, (char *)seq2, len2 };
	sequence_t enseq1 = { 1, malloc(len1 + 1), len1 };
	sequence_t enseq2 = { 2, malloc(len2 + 1), len2 };
	seq2encodedseq(inseq1, enseq1, lal_encode31);
	seq2encodedseq(inseq2, enseq2, lal_encode31);


	sequence_t reverse1 = { 3, malloc(len1 + 1), len1 };
	for (size_t i = 0; i < enseq1.len; i++) {
		// computes the reverse complement of the input sequence.
		reverse1.seq[i] = cns[(int)(enseq1.seq[enseq1.len - 1 - i])];
	}
	char rev1[] = { "CAACTTCCTACGTTGGGTCATAGTAGTGCGTGGGCAATGCCTACGGAGGGGTGGAGCAACTGGCGCTATCACTTCTACCATCGTCTGCAGCGTACGA" };
	size_t revlen1 = strlen(rev1);
	lal_reverse(seq1, len1, rev1, lal_revers31_s);
	//0x003e(rev1)  "TCGTACGCTGCAGACGATGGTAGAAGTGATAGCGCCAGTTGCTCCACCCCTCCGTAGGCATTGCCCACGCACTACTATGACCCAACGTAGGAAGTTG"
	char rev2[] = { "TCGTACGCTGCAGACGATGGTAGAAGTGATAGCGCCAGTTGCTCCACCCCTCCGTAGGCATTGCCCACGCACTACTATGACCCAACGTAGGAAGTTG" };

	for (size_t i = 0; i < reverse1.len; i++) {
		ck_assert_int_eq((int)rev1[i], (int)rev2[i]);
	}

	sequence_t inrev1 = { 3, (char *)rev1, len1 };
	sequence_t enrev1 = { 3, malloc(len1 + 1), len1 };
	seq2encodedseq(inrev1, enrev1, lal_encode31);

	for (size_t i = 0; i < reverse1.len; i++) {
		ck_assert_int_eq((int)reverse1.seq[i], (int)enrev1.seq[i]);
	}
}END_TEST


void addEncodingTC(Suite *s) {
	TCase *tc_core = tcase_create("Encoding");
	tcase_add_test(tc_core, test_encoding_261);
	tcase_add_test(tc_core, test_encoding_290);
	suite_add_tcase(s, tc_core);
}

