#ifndef _CHECK_ALGORITHMS_H_
#define _CHECK_ALGORITHMS_H_

typedef	int pid_t;

#include <check.h>

void addEncodingTC(Suite *s);
void addGSSTC(Suite *s);
void addSWTC(Suite *s);
void addFPTC(Suite *s);
void addScoringMatrixTC(Suite *s);
#endif /* _CHECK_ALGORITHMS_H_ */
