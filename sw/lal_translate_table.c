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
#include <ctype.h>

#include "lal_report.h"
#include "lal_matrix.h"
#include "lal_tables.h"
#include "lal_translate_table.h"

typedef struct {
	char na_symb;  /* the nucleic acid symbol */
	int  na_len;   /* number of unambiguous na in the na */
	char *na_list; /* list of unambiguous na in the na */
} NCBI_NA_t;

#define lal_na2na_len(na) lal_na_s[(int)(lal_na2indx[(int)na])].na_len
#define lal_na2na_list(na, indx) lal_na_s[(int)(lal_na2indx[(int)na])].na_list[indx]

static NCBI_NA_t lal_na_s[NUCL_ACID_NUM] = {
	{ 'N',4,"CTAG" },
	{ 'C',1,"C" },
	{ 'T',1,"T" },
	{ 'Y',2,"CT" },
	{ 'A',1,"A" },
	{ 'M',2,"CA" },
	{ 'W',2,"TA" },
	{ 'H',3,"CTA" },
	{ 'G',1,"G" },
	{ 'S',2,"CG" },
	{ 'K',2,"TG" },
	{ 'B',3,"CTG" },
	{ 'R',2,"AG" },
	{ 'V',3,"CAG" },
	{ 'D',3,"TAG" },
	{ '.',0,NULL } };

typedef struct tag_descritor
{
	char *data;
	char **list;
	size_t current_position;
	size_t length;

}descritor_t;

static descritor_t create_descritor(const char *tablestring, size_t len)
{
	descritor_t dsc;
	dsc.data = malloc(len + 1);
	memcpy(dsc.data, tablestring, len + 1); /*+1 \null simbol for string data*/
	char *Letters = strtok(dsc.data, "\t\n");
	size_t rows = 0;
	while (Letters) {
		Letters = strtok(NULL, "\t\n");
		rows++;
	}
	dsc.length = rows;
	dsc.list = (char **)malloc(rows * sizeof(char*));
	dsc.list[0] = dsc.data;
	for (size_t i = 0; i < dsc.length - 1 /*dispose of heap corruption */; i++) {
		dsc.list[i + 1] = dsc.list[i] + strlen(dsc.list[i]) + 1;
	}
	return dsc;
}

static int read_docs(translate_table_t *tt, descritor_t * desc)
{ /* fills docs  */
	size_t cnt = 0;
	const char *s;
	const char *pshift;
	tt->Doc[0] = '\0';
	desc->current_position = 0;
	for (size_t i = 0; i < desc->length; i++) {
		if (s = desc->list[i]) {
			if (!(pshift = strstr(s, "#"))) {
				desc->current_position = i; /* next line */
				return 1; /* successful docs reading */
			}
			if ((cnt += strlen(pshift)) + 1 < MAX_DOC_LEN)
				strcat(tt->Doc, pshift);
		}
	}
	return 0;
}


static void set_tt_val(char *na3, char symb, TransTable_t TransTable)
{
	TransTable[lal_na2indx[na3[0]]][lal_na2indx[na3[1]]][lal_na2indx[na3[2]]] = toupper(symb);
}

static void set_tt_mval(char *na3, char symb, TransTable_t TransTable)
{
	int i, j, k;
	char new_na3[4];

	new_na3[3] = '\0';
	for (i = 0; i < lal_na2na_len(na3[0]); i++) {
		new_na3[0] = lal_na2na_list(na3[0], i);
		for (j = 0; j < lal_na2na_len(na3[1]); j++) {
			new_na3[1] = lal_na2na_list(na3[1], j);
			for (k = 0; k < lal_na2na_len(na3[2]); k++) {
				new_na3[2] = lal_na2na_list(na3[2], k);
				set_tt_val(new_na3, symb, TransTable);
			}
		}
	}
}

#define lal_indx2na(indx) lal_na_s[indx].na_symb

/*-------------------------------------------------------*
Check whether the nucleic sequence *na contains
an ambiguous nucleic acid.
Note: the sequence with '.' is un ambiguous.
Return:
1 if does, 0 if doesn't.
*-------------------------------------------------------*/
int lal_is_ambiguous_na(char *na, int len)
{
	int i;
	int retval = 0;

	for (i = 0; i < len; i++) {
		if (na[i] == '.') return 0;
		if (lal_na2na_len(na[i]) > 1)
			retval = 1;
	}
	return retval;
}

static char get_tt_val(char *na3, TransTable_t TransTable)
{
	int n1 = lal_na2indx[na3[0]];
	int n2 = lal_na2indx[na3[1]];
	int n3 = lal_na2indx[na3[2]];

	return TransTable[n1][n2][n3];
}
static void UpdateTransTable(char *na3, TransTable_t TransTable)
{
	int i, j, k, check = -1;
	char new_na3[4], symb;

	new_na3[3] = '\0';
	for (i = 0; i < lal_na2na_len(na3[0]); i++) {
		new_na3[0] = lal_na2na_list(na3[0], i);
		for (j = 0; j < lal_na2na_len(na3[1]); j++) {
			new_na3[1] = lal_na2na_list(na3[1], j);
			for (k = 0; k < lal_na2na_len(na3[2]); k++) {
				new_na3[2] = lal_na2na_list(na3[2], k);
				symb = get_tt_val(new_na3, TransTable);
				if (check == -1)
					check = (int)symb;
				else
					if (check != (int)symb)
						return;
			}
		}
	}
	set_tt_val(na3, (char)check, TransTable);
	return;
}

static void CompleteTransTable(TransTable_t TransTable)
{
	int i, j, k;
	char na3[4];
	na3[3] = '\0';
	for (i = 0; i < NUCL_ACID_NUM; i++) {
		na3[0] = lal_indx2na(i);
		for (j = 0; j < NUCL_ACID_NUM; j++) {
			na3[1] = lal_indx2na(j);
			for (k = 0; k < NUCL_ACID_NUM; k++) {
				na3[2] = lal_indx2na(k);
				if (lal_is_ambiguous_na(na3, 3))
					UpdateTransTable(na3, TransTable);
			}
		}
	}
}

static void GetIntTranslateTable(TransTable_t TransTable, IntTransTable_t IntTransTable)
{
	int i, j, k;

	for (i = 0; i<NUCL_ACID_NUM; i++)
		for (j = 0; j<NUCL_ACID_NUM; j++)
			for (k = 0; k<NUCL_ACID_NUM; k++)
				IntTransTable[i][j][k] = lal_encode31[TransTable[i][j][k]];
}

static int deflate_table(translate_table_t *tt, descritor_t * desc)
{
	char na3[4], *token, symb, *s;

	memset(tt->TheTable, 'X', sizeof(TransTable_t));
	for (size_t i = desc->current_position; i < desc->length; i++) {
		if (s = desc->list[i]) {
			if (*s == '#') continue; /* comments */
			if (*s == '\n') continue; /* empty string */
			if (!(token = strtok(s, " \t\n"))) continue;
			strncpy(na3, token, 3);
			if (!(token = strtok(NULL, " \t\n"))) continue;
			symb = toupper(token[0]);
			set_tt_mval(na3, symb, tt->TheTable);
		}
	}
	CompleteTransTable(tt->TheTable);
	GetIntTranslateTable(tt->TheTable, tt->TheIntTable);
	return 1;
}


int read_translate_table(translate_table_t *tt, const char *tablestring, size_t len) {
	descritor_t  desc = create_descritor(tablestring, len);
	double result = -1;
	/* return scale based on analysis of the profile
	The scale is not optimal but works around the following ideas:
	the smallest (in abs. numbers) penalty is scaled at least to '-1', so we still penalize it.
	No we look at the smallest positive score (after the scale) and make sure it's at least 1 as well (so we don't increase it's score too much by rounding).
	*/
	/*end*/

	if (!read_docs(tt, &desc)) {
		report_warning("docs not found");
		return 0;
	}
	deflate_table(tt, &desc);
	return 1;
}

