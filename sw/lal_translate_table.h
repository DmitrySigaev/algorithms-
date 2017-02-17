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
#define NUCL_ACID_NUM 16

typedef unsigned char TransTable_t[NUCL_ACID_NUM][NUCL_ACID_NUM][NUCL_ACID_NUM];
typedef unsigned char IntTransTable_t[NUCL_ACID_NUM][NUCL_ACID_NUM][NUCL_ACID_NUM];

typedef struct tag_translate_table {
	TransTable_t TheTable;
	IntTransTable_t TheIntTable;
	char Doc[MAX_DOC_LEN];
} translate_table_t;


int read_translate_table(translate_table_t *tt, const char *tablestring, size_t len);

#endif /* _LAL_TRANSLATE_TABLE_H_ */