/*
Copyright (C) 2017

Contact: Dmitry Sigaev <dima.sigaev@gmail.com>
*/

#ifndef _LAL_TRANSLATE_TABLE_H_
#define _LAL_TRANSLATE_TABLE_H_

#include <stdint.h>
#include "lal_matrix.h"

#define MAX_LINE_LEN  1024
#define MAX_DOC_LEN   (MAX_LINE_LEN*4)

typedef struct tag_translate_table {
	char Doc[MAX_DOC_LEN];
} translate_table_t;


int read_translate_table(translate_table_t *tt, const char *tablestring, size_t len);

#endif /* _LAL_TRANSLATE_TABLE_H_ */