/*
Copyright (C) 2016

Contact: Dmitry Sigaev <dima.sigaev@gmail.com>
*/

#ifndef _LAL_TYPEDEFS_H_
#define _LAL_TYPEDEFS_H_

#include <stdint.h>
#include "lal_scoring_matrix.h"
/**
 *  Structure of a sequence
 *  @typedef {struct} sequence_t
 *  @field ID   identifier of sequence ( currently used in fasta)
 *  @field seq   the sequence
 *  @field len   length of the sequence
 */
typedef struct tag_sequence {
	size_t ID;
	char * seq;
	size_t len;
} sequence_t;

/**
*  Structure of a search profile for smith-waterman algorithm with constant gap
*  @typedef {struct} search_profile_t
*  @field gap   gap opening cost and gap extension cost
*  @field score_matrix   score matrix
*/
typedef struct tag_search_swcg_profile
{
	double gap;
	scoring_matrix_t *mtx; /*can be NULL then use SCORE define*/
} search_swcg_profile_t;

/**
*  Structure of a search profile for smith-waterman algorithm with constant gap
*  @typedef {struct} search_profile_t
*  @field gap   gap opening cost and gap extension cost
*  @field score_matrix   score matrix
*/
typedef struct tag_search_swcg_profile_int
{
	int64_t gap;
	scoring_matrix_t *mtx;
} search_swcg_profile_int_t;

/**
*  Structure of a search profile for smith-waterman algorithm with affine gap
*  @typedef {struct} search_profile_t
*  @field gap   gap opening cost and gap extension cost
*  @field score_matrix   score matrix
*/
typedef struct tag_search_swag_profile
{
	double gapOpen;
	double gapExt;
	scoring_matrix_t *mtx;
} search_swag_profile_t;

#endif /* _LAL_TYPEDEFS_H_ */