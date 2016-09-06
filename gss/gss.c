#include <stdio.h>
#include "gss.h"

/*
 * Given a sequence of integers, find a subsequence which maximizes
 * the sum of its elements, that is, the elements of no other single
 * subsequence add up to a value larger than this one.
 */
range_t maxSubseq(const int sequence[], const int len) {
	range_t r ;
	int maxSum = 0, thisSum = 0, i = 0;
	int start = 0, end = -1, j;

	for (j = 0; j < len; j++) {
		thisSum += sequence[j];
		if (thisSum < 0) {
			i = j + 1;
			thisSum = 0;
		}
		else if (thisSum > maxSum) {
			maxSum = thisSum;
			start = i;
			end = j;
		}
	}

	if (start <= end && start >= 0 && end >= 0) {
		r.start = start;
		r.end = end + 1;
		r.sum = maxSum;
	}
	else {
		r.start = 0;
		r.end = 0;
		r.sum = 0;
	}
	return r;
}

