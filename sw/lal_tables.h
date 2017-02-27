/*
Copyright (C) 2016

Contact: Dmitry Sigaev <dima.sigaev@gmail.com>
*/

#ifndef _LAL_TABLES_H_
#define _LAL_TABLES_H_



/*
* Contains the numerical representation of the complementary nucleotide symbol.
*/

extern const char cns[31];
extern const unsigned char lal_encode31[];
extern const char lal_revers31[256];
extern const char lal_decode31[31];
extern const unsigned char lal_na2indx[256];

extern const char identity_nuc[];
extern const char blosum100[];
extern const char gaptest1[];
extern const char blosum62[];
extern const char human40[];

#endif /* _LAL_TABLES_H_ */
