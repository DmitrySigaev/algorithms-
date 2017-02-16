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

#include "lal_report.h"
#include "lal_matrix.h"
#include "lal_tables.h"
#include "lal_translate_table.h"


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

	return 1;
}

