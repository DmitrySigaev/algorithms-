/*
Copyright (C) 2016

Contact: Dmitry Sigaev <dima.sigaev@gmail.com>
*/
#include <stdio.h>
#include "lal_typedefs.h"
#include "lal_report.h"
#include "lal_encoding.h"
#include "lal_translate_table.h"

/**
* Returns a new created sequence, with the mapped values of the input sequence.
*
* @param inseq   the input sequence that should be coded
* @param enseq   the new coded sequence
* @param encode  the encode table
*/
void lal_seq2encodedseq(sequence_t inseq, sequence_t enseq, const char* encode) {
	size_t illicit_symbol = 0;
	for (size_t i = 0; i < inseq.len; i++) {
		char cm;
		if ((cm = encode[(size_t)inseq.seq[i]]) >= 0) {
			enseq.seq[i] = cm;
		}
		else {
			enseq.seq[i] = LAL_NONEXISTENT;  /*if encode table has -1 value*/
			illicit_symbol++;
		}
	}
	if (illicit_symbol > 0) {
		report_warning("illicit symbols found and set to LAL_NONEXISTENT");
	}
	enseq.seq[enseq.len] = LAL_NONEXISTENT;
}

/**
* Returns a new created sequence, with the mapped values of the input sequence.
* Can insert amino-acid translation of nucleic sequence into seq structure
*
* @param inseq   the input sequence that should be coded
* @param enseq   the new coded sequence
* @param encode  the encode table (nucleic triplets)
* @param tt      the translate table (amino-acid translation)
*/
void lal_seq2encodedseq_trans(sequence_t inseq, sequence_t enseq, const char* encode, translate_table_t *tt) {
	size_t illicit_symbol = 0;
	/* code the amino - acids from nucleic triplets */
	for (int i = 0; i < inseq.len + 1; i++) {
		int i2 = ((i - 1 - 2) >= 0) ? encode[inseq.seq[i - 1 - 2]] : 0;
		int i1 = ((i - 1 - 1) >= 0) ? encode[inseq.seq[i - 1 - 1]] : 0;
		int i0 = ((i - 1 - 0) >= 0) ? encode[inseq.seq[i - 1 - 0]] : 0;
		enseq.seq[i] = tt->TheIntTable[i2][i1][i0];
	}
//	enseq.seq[enseq.len] = LAL_NONEXISTENT;
}

/**
* Create reversed strands
*
* @param source  the input sequence that should be reversed
* @param len     the length of sequnce
* @param enseq   the new revers sequence
* @param encode  the reverse table 
*/
void lal_reverse(const char * source, size_t len, char *dest, const char *reverse_tab)
{
	char *ptmp;
	if (!len) return;
	ptmp = dest + (len - 1);
	while (len--)
		*ptmp-- = reverse_tab[*source++];
}
