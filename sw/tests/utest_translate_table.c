/*
Copyright (C) 2016

Contact: Dmitry Sigaev <dima.sigaev@gmail.com>
*/

#include "../../utests/tests.h"

#include <stdint.h>
#include <stdio.h>
#include "../lal_tables.h"
#include "../lal_translate_table.h"


char badttable[] = { "#Genetic code file > Standard genetic code   \n \
#base A frequency 0.313393                   \n \
#base C frequency 0.206870                   \n \
#base G frequency 0.217274                   \n \
#base T frequency 0.262463                   \n \
#G + C = 0.424144                              \n \
#codon aa freq freq / aa odds  10 * bits         \n" };


START_TEST(test_bad_table)
{
	translate_table_t tt;
	int status = read_translate_table(&tt, badttable, strlen(badttable));
	ck_assert_int_eq(status, 0); /* bad format*/
}END_TEST

START_TEST(test_human40)
{
	translate_table_t tt;
	int status = read_translate_table(&tt, human40, strlen(human40));
	ck_assert_int_eq(status, 0); /* bad format*/
}END_TEST

void addTranslateTableTC(Suite *s) {
	TCase *tc_core = tcase_create("Translate_table");
	tcase_add_test(tc_core, test_bad_table);
	tcase_add_test(tc_core, test_human40);

	suite_add_tcase(s, tc_core);
}

