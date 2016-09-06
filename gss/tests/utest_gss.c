/*
Copyright (C) 2016

Contact: Dmitry Sigaev <dima.sigaev@gmail.com>
*/

#include "../../utests/tests.h"

#include <stdint.h>
#include <stdio.h>
#include "../gss.h"


START_TEST(test_gss)
{
	int a[] = { -1 , -2 , 3 , 5 , 6 , -2 , -1 , 4 , -4 , 2 , -1 };
	int alength = sizeof(a) / sizeof(a[0]);

	range_t r = maxSubseq(a, alength);
	ck_assert_int_eq(maxSubseq(a, alength).sum, 15); /* Max sum */
	ck_assert_int_eq(r.start, 2); /* start position*/
	ck_assert_int_eq(r.end, 8); /* end position*/

}END_TEST

void addGSSTC(Suite *s) {
	TCase *tc_core = tcase_create("GSS");
	tcase_add_test(tc_core, test_gss);

	suite_add_tcase(s, tc_core);
}

