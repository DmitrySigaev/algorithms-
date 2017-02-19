/*
Copyright (C) 2016 Dmitry Sigaev

Contact: Dmitry Sigaev <dima.sigaev@gmail.com>
*/


#ifndef _LAL_ENCODING_H_
#define _LAL_ENCODING_H_

#include <stdint.h>
#include "lal_typedefs.h"
#include "lal_translate_table.h"

#define LAL_NONEXISTENT (unsigned char)31

void seq2encodedseq(sequence_t inseq, sequence_t enseq, const char* encode);
void seq2encodedseq_trans(sequence_t inseq, sequence_t enseq, const char* encode, translate_table_t *tt);
void lal_reverse(const char * source, int len, char *dest, const char *reverse_tab);

#endif /* _LAL_ENCODING_H_ */