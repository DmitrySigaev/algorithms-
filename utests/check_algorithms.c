/*
Copyright (C) 2016

Contact: Dmitry Sigaev <dima.sigaev@gmail.com>
*/


#include "tests.h"

#include <stdlib.h>
#include <stdio.h>

Suite* algorithms_suite(void)
{
	Suite *s = suite_create("algorithms");

	/* Core test case */
	addEncodingTC(s);
	addGSSTC(s);
	addSWTC(s);
	addScoringMatrixTC(s);

	return s;
}

int main(void)
{
	printf("Using Check unit testing framework version %d.%d.%d\n", CHECK_MAJOR_VERSION, CHECK_MINOR_VERSION,
		CHECK_MICRO_VERSION);

	int number_failed;
	Suite *s = algorithms_suite();
	SRunner *sr = srunner_create(s);
	srunner_run_all(sr, CK_NORMAL);
	number_failed = srunner_ntests_failed(sr);
	srunner_free(sr);
	return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
