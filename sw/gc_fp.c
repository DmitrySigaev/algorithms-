/*
Copyright (C) 2017

Contact: Dmitry Sigaev <dima.sigaev@gmail.com>
*/

#include <stdio.h>
#include <malloc.h>
#include <limits.h>
#include <math.h>

#include "lal_matrix.h"
#include "lal_report.h"
#include "gc_sw.h"

/* type of the profile and calc score for MS_ONE */
typedef int scoretype;

typedef struct {
	scoretype *data; /* array as if of [prof_length][this.line_size]*/
	int line_size;   /* size of the profile line */
					 /*  float scaleback;*/ /* factor to scale results back to float */
	double scaleback; /* factor to scale results back to float */
} ms_profile_t;

/* Type of the best result, current in and returned by the
internal score_calc function.
Contains score and end coordinates of best alignment.
*/
typedef struct best {
	scoretype score; /* alignment score */
					 /*  int state;*/
	int x; /* alignment's end coordinates */
	int y;
} best_t;

/* length of the match table section of a profile line.
In HMM2 these are match and insert sections of this length.
There are 32 scores in match section for 32 letters and symbols possible
in a sequence */
#define MS_PROFILE_MATCHTABLE_SIZE 32

/* each profile line id 38 floats
( 32 emission scores + 14 gap penalties : matchmax5 , matchmax6 , xopen , xext , fx, fopen , fext , 
	  yopen1 , yopen2 , yopen3 , yext, xopen1 , xext1 , fx1 ;
  profile : PROF_SCORES , $matchmax5 , $matchmax6 , $xopen , $xext , $fx, $fopen , $fext ,
	  $yopen1 , $yopen2 , $yopen3 , $yext, $xopen1 , $xext1 , $fx1 ;

 */
#define FP_PROFILE_LINE_SIZE 46

/* accessor for element j in line i of profile pr */
#define PROFILE_VAL(pr,i,j) pr[(i)*FP_PROFILE_LINE_SIZE+(j)]

#define PROFILE_EMISSIONS_NO 32


enum MS_FRAMEPLUS_OFFSETS {
	MS_FRAMEPLUS_PROF_SCORES = 0,
	MS_FRAMEPLUS_OFFSET_matchmax5 = MS_FRAMEPLUS_PROF_SCORES + MS_PROFILE_MATCHTABLE_SIZE,
	MS_FRAMEPLUS_OFFSET_matchmax6,

	/* @@ next 3 fields are not used */
	MS_FRAMEPLUS_OFFSET_xopen, MS_FRAMEPLUS_OFFSET_xext, MS_FRAMEPLUS_OFFSET_fx,

	MS_FRAMEPLUS_OFFSET_fopen,
	MS_FRAMEPLUS_OFFSET_fext, MS_FRAMEPLUS_OFFSET_yopen1,
	MS_FRAMEPLUS_OFFSET_yopen2, MS_FRAMEPLUS_OFFSET_yopen3,
	MS_FRAMEPLUS_OFFSET_yext, MS_FRAMEPLUS_OFFSET_xopen1,
	MS_FRAMEPLUS_OFFSET_xext1, MS_FRAMEPLUS_OFFSET_fx1
};

enum { MS_FRAMEPLUS_OFFSET_DEPTH_P2N = MS_FRAMEPLUS_OFFSET_fx1 + 1 }; /* depth of frame+_p2n profile in data4vec */
#define min_2(a, b)  (((a)<(b)) ? (a) : (b))/* b will be evaluated twice! */

/* Used as minus infinity in creating profile when a score of below
MS_MIN_PROFILE_CUTOFF is found in GenProf.
Also used in initializing transitions to minus infinity.
Type: scoretype
Need to be able to add 8 times (gw9)
Should be able to beat any MS_ONE_MAX_TOTAL
*/
#define MS_SCORE_MININF    ((scoretype)(INT_MIN/16))

/* Type of the sequence element.
Code of letter: (capital_char - 'A')
Range: 0..31
Used as offset in matchtable part of the profile line to get match score */
typedef int seqtype;

#define scoretype_max2(a, b) {register scoretype c=(b); if ((a)<c) (a)=c; }

/* Will be inlined only with linux/gcc/-O3
#define scoretype_max(a, b) {a = scoretype_max_(a, (b)); }
*/

#define scoretype_max(a, b) scoretype_max2(a, b)

#define checkbest(quality, xx, yy, max_v) \
                           if ((quality)>max_v.score) {\
                             max_v.score=(quality);   \
                             max_v.x=(xx);            \
                             max_v.y=(yy);            \
                           }


typedef struct prev_line_s {
	scoretype yskipmatch; /* delete state */
	scoretype prequal; /* match state */
	seqtype   nseq; /* letter of the sequence */
} prev_line_t;

/* this is the original MS mode */
static  best_t MS_Score_FP(seqtype* nseq, ms_profile_t *prof, int profLen, int seqLen)
{
	int x = 0, y = 0; /* x - coordinate in seq, y - coord in profile */


	scoretype xskipmatch[4], quality[4], fskipmatch[4]; /* insert, match, frame-shift insert states */
	scoretype ytest[4]; /* delete state */
	scoretype  q_new; /* match state score to be computed */
	scoretype  y_new; /* delete state to be computed */
	scoretype score; /* match score from profile*/

	prev_line_t *prev_line; /* array that stores state scores for whole line of DP matrix */
	prev_line_t *pre; /* pointer in prev_line array */
	scoretype *prof_line_p;

	scoretype matchmax5, matchmax6, /*xopen , xext , fx, */ f_open, fext,
		yopen1, yopen2, yopen3, yext, xopen1, xext1, fx1; /* gap scores from profile */
	int k0, k1, k2, k3, k; /* indeces for previous positions in insert&Co state scores buffers */


#define _COUNT 0
#if _COUNT
	static int count1 = 1, count2 = 1, count3 = 1, count12 = 1, count13 = 1, count123 = 1, count0 = 2;
#endif

	best_t max_v; /* best seen alignment score in a match state */


	max_v.score = MS_SCORE_MININF;
	max_v.x = -1;
	max_v.y = -1;

	/* return if shorter than 2 bases */
	if (seqLen <= 3) {
		return max_v;
	}

	prev_line = (prev_line_t *)calloc((seqLen + 1), sizeof(prev_line_t));

	prev_line[seqLen].nseq = -1; /* end of sequence token */

								 /* init prev_line by the scores for first line in the profile */

	for (x = 0; x<seqLen; x++)
	{
		prev_line[x].nseq = nseq[x];
		prev_line[x].prequal = prof->data[prev_line[x].nseq];
		prev_line[x].yskipmatch = MS_SCORE_MININF;
		checkbest(prev_line[x].prequal, x, 0, max_v);
	}

	/* loop on profile lines 1..(profLen-1) */
	/* prof_line_p already set to prof->data */
	prof_line_p = prof->data;

	for (y = 1; y<profLen; y++)
	{
		/* init sequence/state array */
		pre = prev_line + 3; /* first 3 are empty anyway */

							 /* gap penalties for line n are stored in line n-1 */
		xopen1 = prof_line_p[MS_FRAMEPLUS_OFFSET_xopen1];
		xext1 = prof_line_p[MS_FRAMEPLUS_OFFSET_xext1];
		fx1 = prof_line_p[MS_FRAMEPLUS_OFFSET_fx1];


		/* go to the next profile line */
		prof_line_p = prof->data + y * prof->line_size; /* prof_line_p+=prof->line_size; */

		matchmax5 = prof_line_p[MS_FRAMEPLUS_OFFSET_matchmax5];
		matchmax6 = prof_line_p[MS_FRAMEPLUS_OFFSET_matchmax6];
		yopen1 = prof_line_p[MS_FRAMEPLUS_OFFSET_yopen1];
		yopen2 = prof_line_p[MS_FRAMEPLUS_OFFSET_yopen2];
		yopen3 = prof_line_p[MS_FRAMEPLUS_OFFSET_yopen3];
		yext = prof_line_p[MS_FRAMEPLUS_OFFSET_yext];

		/* constant */
		fext = prof_line_p[MS_FRAMEPLUS_OFFSET_fext];
		f_open = prof_line_p[MS_FRAMEPLUS_OFFSET_fopen]; /* fopen is a std. function name too ? */



		for (k = 0; k<3; k++)
		{
			xskipmatch[k] = MS_SCORE_MININF;
			fskipmatch[k] = MS_SCORE_MININF;
			ytest[k] = MS_SCORE_MININF;
			quality[k] = MS_SCORE_MININF;
		}

		k0 = 0; k1 = 1; k2 = 2; k3 = 3;

		/*x=3, put into k2 */
		quality[k2] = pre->prequal;
		q_new = prof_line_p[pre->nseq];
		checkbest(q_new, 3, y, max_v);
		pre->prequal = q_new;

		ytest[k2] = pre->yskipmatch;
		y_new = ytest[k2] + yext;
		scoretype_max(y_new, quality[k2] + yopen1);
		pre->yskipmatch = y_new;

		pre++;

		/* loop on sequence positions */
		/* x = 4..(seqLen-1) */
		while (pre->nseq >= 0)
		{
			score = prof_line_p[pre->nseq];
			quality[k3] = pre->prequal;
			ytest[k3] = pre->yskipmatch;

#if _COUNT
			count0++;
			if (xskipmatch <= 0) count1++;
			if (ytest <= 0) count2++;
			if (xskipmatch <= 0 && ytest <= 0) count12++;
			if (quality <= pentest) count3++;
			if (quality <= pentest && ytest <= 0) count13++;
			if (xskipmatch <= 0 && ytest <= 0 && quality <= pentest) count123++;
#endif


			/* xskipmatch[k3] */

			xskipmatch[k3] = xskipmatch[k0] + xext1;
			scoretype_max(xskipmatch[k3], fskipmatch[k0] + fx1);
			scoretype_max(xskipmatch[k3], quality[k0] + xopen1);

			/* fskipmatch */
			fskipmatch[k3] = fskipmatch[k2] + fext;
			scoretype_max(fskipmatch[k3], quality[k2] + f_open);

			/* match */
			q_new = 0; /* V: Local alignment? */
			scoretype_max(q_new, quality[k0]);
			scoretype_max(q_new, ytest[k0]);
			scoretype_max(q_new, xskipmatch[k0]);
			scoretype_max(q_new, fskipmatch[k0]);

			q_new += score;

			/* These are not exactly Match transitions, do not add score */
			scoretype_max(q_new, quality[k1] + matchmax5);
			scoretype_max(q_new, quality[k2] + matchmax6);

			/* y - now can rewrite y */
			y_new = ytest[k3] + yext;
			scoretype_max(y_new, quality[k3] + yopen1);
			scoretype_max(y_new, quality[k2] + yopen2);
			scoretype_max(y_new, quality[k1] + yopen3);

			pre->yskipmatch = y_new;
			pre->prequal = q_new;

			checkbest(q_new, pre - prev_line, y, max_v);
			{
				int i_tmp;
				i_tmp = k0; k0 = k1; k1 = k2; k2 = k3; k3 = i_tmp;
				pre++;
			}
		} /* sequence loop */

		  /*   checkbest(quality[k2], seqLen, y+1, max_v)*/

	} /* profile loop */

	  /*   for(x=1, pre = prev_line+1;x<seqLen;x++) */
	  /*   { */
	  /*     checkbest(pre->prequal, x, profLen, max_v) */
	  /*     pre++; */
	  /*   } */

	  /*  max_v.x--;max_v.y--;*/

	free(prev_line);

	{
#if _COUNT    
		float a1, a2, a3, a12, a13;
		a1 = 1.0*count1 / count0;
		a2 = 1.0*count2 / count0;
		a12 = 1.0*count12 / count0;
		a3 = 1.0*count13 / count0;
		a13 = 1.0*count13 / count0;
		printf("\n\nSWAT:\ncount 1x\t%.2f%% out of %d, %d\n", (100.0*a1), count0, count1);
		printf("count 2y\t%.2f%% out of %d, %d\n", (100.0*a2), count0, count2);
		printf("count 12xy\t%.2f%% out of %d, %d\t=c %.4f%%\n", (100.0*a12), count0, count12, 100.0*(a12 - a1*a2) / sqrt(a1*(1 - a1)*a2*(1 - a2)));
		printf("count 3q\t%.2f%% out of %d, %d\n", (100.0*a3), count0, count3);
		printf("count 13xq\t%.2f%% out of %d, %d\t=c %.4f%%\n", (100.0*a13), count0, count13, 100.0*(a13 - a1*a3) / sqrt(a1*(1 - a1)*a3*(1 - a3)));
		printf("count 123xyq\t%.2f%% out of %d, %d\n\n", (100.0*count123 / count0), count0, count123);
#endif
	}

	return max_v;
} /* MS_Score_FP() */

typedef int score_t;

static score_t max2(score_t a, score_t b)
{
	return (a < b) ? b : a;
}
static score_t max3(score_t a, score_t b, score_t c)
{
	return max2(a, max2(b, c));
}
static score_t max4(score_t a, score_t b, score_t c, score_t d)
{
	return max2(max2(a, b), max2(c, d));
}

static score_t max5(score_t a, score_t b, score_t c, score_t d, score_t e)
{
	return max2(e, max2(max2(a, b), max2(c, d)));
}

/*Matrix of Frame plus algorithm*/
best_t score_frameplus_p2n(seqtype* nseq, ms_profile_t *prof, int profLen, int seqLen)
{
	/* What I know about the sequence:
	1. before we do anything, nucleic seq. is prepened by a blank symbol.
	2. Next, each position [i] is subst. with AA corresponding to [i-2,i-1, i]. It follows that the first 3 positions are always 'X' (no corr. AA).
	3. Stop codon has code 27 (corr to '*' in blosum & oth matrices)
	*/
	enum state_t { s_start = 0, s_fgap, s_xgap, s_ygap, s_match, s_tmpmatch, s_end, s_state_size };
	/* Datastructure: */
	const int seq_length = seqLen;
	const int *translated_sequence = nseq; /* translated_sequence[i] = AA of 3 nucl in original query */
	const int profile_len = profLen;

	score_t  ***mattr;
	static int called = 0;

	int i, j;
	score_t max_score;
	best_t max_v; /* best seen alignment score in a match state */

	mattr = malloc(seq_length * sizeof(score_t *));

	for (i = 0; i < seq_length; i++)
	{
		mattr[i] = malloc(profile_len * sizeof(mattr[0][0]));
		for (j = 0; j < profile_len; j++)
			mattr[i][j] = malloc(s_state_size * sizeof(mattr[0][0][0]));
	}

#define profile(prof, y, x) ( (prof)->data[ (y)*(prof)->line_size+(x) ] )

	max_v.score = 0;
	max_v.x = -1;
	max_v.y = -1;

	for (i = 0; i < seq_length; i++)
		for (j = 0; j < profile_len; j++)
			mattr[i][j][s_fgap] = mattr[i][j][s_xgap] = mattr[i][j][s_ygap] = mattr[i][j][s_match] = MS_SCORE_MININF;

	for (i = 3; i < seq_length; i++)
	{
		mattr[i][0][s_fgap] = mattr[i][0][s_xgap] = mattr[i][0][s_ygap] = mattr[i][0][s_match] = MS_SCORE_MININF;
		for (j = 1; j < profile_len; j++)
		{
			/* gap penalties for line n are stored in line n-1 */

			/* Note: the code below was modified to read profile at [j][xopen1] instead of [j-1][xopen1]
			the reasons are to be found out, but this seems to be the right place to read it to be consistent with om/ms runmodes */

			const score_t xopen1 = profile(prof, j, MS_FRAMEPLUS_OFFSET_xopen1); /*prof_line_p[MS_FRAMEPLUS_OFFSET_xopen1]; */
			const score_t xext1 = profile(prof, j, MS_FRAMEPLUS_OFFSET_xext1); /* prof_line_p[MS_FRAMEPLUS_OFFSET_xext1]; */

			const score_t fx1 = profile(prof, j, MS_FRAMEPLUS_OFFSET_fx1);  /* prof_line_p[MS_FRAMEPLUS_OFFSET_fx1]; */

																			/* go to the next profile line */
																			/* const score_t prof_line_p = prof->data+ y * prof->line_size; */ /* prof_line_p+=prof->line_size; */

			const score_t matchmax5 = profile(prof, j, MS_FRAMEPLUS_OFFSET_matchmax5); /* prof_line_p[MS_FRAMEPLUS_OFFSET_matchmax5]; */
			const score_t matchmax6 = profile(prof, j, MS_FRAMEPLUS_OFFSET_matchmax6); /* prof_line_p[MS_FRAMEPLUS_OFFSET_matchmax6]; */
			const score_t yopen1 = profile(prof, j, MS_FRAMEPLUS_OFFSET_yopen1); /* prof_line_p[MS_FRAMEPLUS_OFFSET_yopen1]; */
			const score_t yopen2 = profile(prof, j, MS_FRAMEPLUS_OFFSET_yopen2); /* prof_line_p[MS_FRAMEPLUS_OFFSET_yopen2]; */
			const score_t yopen3 = profile(prof, j, MS_FRAMEPLUS_OFFSET_yopen3); /* prof_line_p[MS_FRAMEPLUS_OFFSET_yopen3]; */
			const score_t yext = profile(prof, j, MS_FRAMEPLUS_OFFSET_yext); /* prof_line_p[MS_FRAMEPLUS_OFFSET_yext]; */
																			 /* constant */
			const score_t fext = profile(prof, j, MS_FRAMEPLUS_OFFSET_fext); /* prof_line_p[MS_FRAMEPLUS_OFFSET_fext]; */
			const score_t f_open = profile(prof, j, MS_FRAMEPLUS_OFFSET_fopen); /* prof_line_p[MS_FRAMEPLUS_OFFSET_fopen]; */ /* fopen is a std. function name too ? */

			const score_t match = profile(prof, j, nseq[i]);

			/* depends on [i-1][j] */
			mattr[i][j][s_fgap] = max2(mattr[i - 1][j][s_match] + f_open, 
				mattr[i - 1][j][s_fgap] + fext);

			/* depends on [i-3][j] */
			mattr[i][j][s_xgap] = max3(mattr[i - 3][j][s_xgap] + xext1,
				mattr[i - 3][j][s_match] + xopen1,
				/*	    xopen1,  -- match above is non-negative*/ /* from s_start */
				mattr[i - 3][j][s_fgap] + fx1);
			/* depends on [i-2..i][j-1] */
			mattr[i][j][s_ygap] = max5(mattr[i - 1][j - 1][s_match] + yopen2,
				mattr[i][j - 1][s_match] + yopen1,
				mattr[i - 2][j - 1][s_match] + yopen3,
				mattr[i][j - 1][s_ygap] + yext,
				0 /*from s_start */);
			/* depends on [i-3][j-1] */
			mattr[i][j][s_match] = match + max5(mattr[i - 3][j - 1][s_fgap],
				mattr[i - 3][j - 1][s_xgap],
				mattr[i - 3][j - 1][s_ygap],
				mattr[i - 3][j - 1][s_match],
				0);
			mattr[i][j][s_match] = max4(0,
				mattr[i][j][s_match],
				mattr[i - 2][j - 1][s_match] + matchmax5,
				mattr[i - 1][j - 1][s_match] + matchmax6);

			if (max_v.score < mattr[i][j][s_match])
			{
				max_v.score = mattr[i][j][s_match];
				max_v.x = i;
				max_v.y = j;
			}
		}
	}
	for (i = 0; i < seq_length; i++)
	{
		for (j = 0; j < profile_len; j++)
			free(mattr[i][j]);
		free(mattr[i]);
	}
	return max_v;
}



#define MIN_DSP_NEGATIVE -1000
  /*** Gencore profile ***/
  /* min_inf in GenProf profiles.
  This value is replaced with MS_SCORE_MININF in MS_ONE Profiles
  Set in Gencore to (-1000) */
#define MS_MIN_PROFILE_CUTOFF MIN_DSP_NEGATIVE


  /* convert x into scoretype value sTmp1, and compute relative error:
  if (x != 0)
  if (x-sTmp1) != 0)
  if (fabs((x-sTmp1)/x) > max_rel_error)
  max_rel_error = fabs((x-sTmp1)/x)

  */
#define to_scoretype_(ss)  ((scoretype)rint(ss))

static scoretype to_scoretype(float data, float scale, float *max_rel_error) {
	float ss = data * scale;
	scoretype sTmp1;

	sTmp1 = to_scoretype_(ss);
	if (ss != 0) {
		float fTmp2;
		/*     if (fabs(data) < 1e-6) {
		fTmp2 = 2048*fabs(data);
		} else */

		/* continuous function: but 0..1/2 map to 0..1, not straight 1 */
		if (sTmp1 == 0) {
			fTmp2 = 2 * fabs(ss);
		}
		else fTmp2 = fabs((ss - sTmp1) / ss);

		/* for data < 1e-6 reduce perceived error as rounding error is very big
		by now anyway. Still a continuous function. */
		if (fabs(data) < 1 / 0x10000) {
			fTmp2 *= fabs(data) * 0x10000;
		}

		if ((fTmp2 > *max_rel_error)
			/* 	 ((fabs(data) > 0.0009765625 ) /-* 1/1024 *-/ */
			/* 	 || (fabs(ss) > 0.0999 ))  /-* 10% *-/ */
			) {
			/*printf("DEBUG: to_scoretype() data=%.10f ss=%f fTmp2=%.10f max_rel_error=%.10f\n", data, ss, fTmp2, *max_rel_error); /-* debug */
			*max_rel_error = fTmp2;
		}
	}
	return sTmp1;
} /* to_scoretype()  */

  /* Return scoretype values into the new profile */
static scoretype get_scaled(float data, float scale, float *max_rel_error) {
	if (data <= MS_MIN_PROFILE_CUTOFF) {
		return MS_SCORE_MININF;
	}
	else {
		return to_scoretype(data, scale, max_rel_error);
	}
} /* get_scaled() */

#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define SCORE(a, b, c, d) (((a) == (b)) ? (c) : (d))
#define VDTABLE(a, b) (sp->mtx->sc_double_matrix.ddata[(a)][(b)])
#define VITABLE(a, b) (sp->mtx->sc_int_matrix.idata[(a)][(b)])

  /*
  *structure for beginning of calculation
  * four parameters should be passed to function
  */
typedef struct inpf {
	float score;
	int state;
	int x;
	int y;
} inpf_t;

double fp_gencore(const search_fp_profile_t * sp, const sequence_t * dseq, const sequence_t * qseq)
{
	/* model->profile_line_size = 66;  len of raw prof*/
	//MS_gcprofile_parts_bounds(model, &read_base, &line_size); gets own frame plus profile from model file
	int read_base = 20; /*shift in raw profile*/
	int line_size = 46; /*32 matrix + 14 gaps*/
	float *data_from;
	scoretype *data_to;
	float scale = sp->mtx->scale;
	float max_rel_error = 0;
	ms_profile_t * ms_profile = malloc(sizeof(ms_profile_t));
	ms_profile->scaleback = 1.0 / sp->mtx->scale;
	ms_profile->line_size = line_size;
	int prof_length = qseq->len + 1; // INT_MIN + seq(i-1)
	int prf_len = (prof_length)* FP_PROFILE_LINE_SIZE * sizeof(float);
	float *prof_data = (float *)malloc(prf_len);

	for (size_t j = 0; j < PROFILE_EMISSIONS_NO; j++)
		PROFILE_VAL(prof_data, 0/* *32+6=38*/, j) = INT_MIN; /*sigaev: global aligment*/

															 /* gap penalties copied from first line */
	prof_data[MS_FRAMEPLUS_OFFSET_matchmax5] = sp->matchMax5;	//32
	prof_data[MS_FRAMEPLUS_OFFSET_matchmax6] = sp->matchMax6;	//33
	prof_data[MS_FRAMEPLUS_OFFSET_xopen1] = sp->gapOpen;		//43
	prof_data[MS_FRAMEPLUS_OFFSET_xext1] = sp->gapExt;			//44
	prof_data[MS_FRAMEPLUS_OFFSET_fopen] = sp->gapFrameOpen;	//37
	prof_data[MS_FRAMEPLUS_OFFSET_fext] = sp->gapFrameExt;		//38
	prof_data[MS_FRAMEPLUS_OFFSET_yopen1] = sp->gapOpen;		//39
	prof_data[MS_FRAMEPLUS_OFFSET_yopen2] = sp->gapOpen2;		//40
	prof_data[MS_FRAMEPLUS_OFFSET_yopen3] = sp->gapOpen3;		//41
	prof_data[MS_FRAMEPLUS_OFFSET_yext] = sp->gapExt;			//42
	prof_data[MS_FRAMEPLUS_OFFSET_fx1] = sp->gapFrame;			//45 latest
	/* @@ next 3 fields are not used */
	//	MS_FRAMEPLUS_OFFSET_xopen, MS_FRAMEPLUS_OFFSET_xext, MS_FRAMEPLUS_OFFSET_fx,
	prof_data[MS_FRAMEPLUS_OFFSET_xopen] = 0.0;					//34
	prof_data[MS_FRAMEPLUS_OFFSET_xext] = 0.0;					//35
	prof_data[MS_FRAMEPLUS_OFFSET_fx] = 0.0;					//35


	for (size_t i = 1; i < prof_length; i++) {
		int k = qseq->seq[i - 1];
		for (size_t j = 0; j < PROFILE_EMISSIONS_NO; j++) { //#define PROFILE_EMISSIONS_NO 32
			double v;
			if (!sp->mtx)
				v = SCORE(k, j, 1.0, -1.0);
			else
				v = VDTABLE(k, j);
//			printf(" %f, ", v);
			PROFILE_VAL(prof_data, i/* *32+6=38*/, j) = v; //create profile line from  32 CM("A") to 38 (36-39: 1,1,1,1,0,0)
		}
																				  // prof_data[(i)*PROFILE_LINE_SIZE + (j)] == prof_data[(i)*38 + (j)]

		PROFILE_VAL(prof_data, i, MS_FRAMEPLUS_OFFSET_matchmax5) = sp->matchMax5;	//32
		PROFILE_VAL(prof_data, i, MS_FRAMEPLUS_OFFSET_matchmax6) = sp->matchMax6;	//33
		PROFILE_VAL(prof_data, i, MS_FRAMEPLUS_OFFSET_xopen1) = sp->gapOpen;		//43
		PROFILE_VAL(prof_data, i, MS_FRAMEPLUS_OFFSET_xext1) = sp->gapExt;			//44
		PROFILE_VAL(prof_data, i, MS_FRAMEPLUS_OFFSET_fopen) = sp->gapFrameOpen;	//37
		PROFILE_VAL(prof_data, i, MS_FRAMEPLUS_OFFSET_fext) = sp->gapFrameExt;		//38
		PROFILE_VAL(prof_data, i, MS_FRAMEPLUS_OFFSET_yopen1) = sp->gapOpen;		//39
		PROFILE_VAL(prof_data, i, MS_FRAMEPLUS_OFFSET_yopen2) = sp->gapOpen2;		//40
		PROFILE_VAL(prof_data, i, MS_FRAMEPLUS_OFFSET_yopen3) = sp->gapOpen3;		//41
		PROFILE_VAL(prof_data, i, MS_FRAMEPLUS_OFFSET_yext) = sp->gapExt;			//42
		PROFILE_VAL(prof_data, i, MS_FRAMEPLUS_OFFSET_fx1) = sp->gapFrame;			//45 latest
		/* @@ next 3 fields are not used */
		//	MS_FRAMEPLUS_OFFSET_xopen, MS_FRAMEPLUS_OFFSET_xext, MS_FRAMEPLUS_OFFSET_fx,
		PROFILE_VAL(prof_data, i, MS_FRAMEPLUS_OFFSET_xopen) = 0.0;					//34
		PROFILE_VAL(prof_data, i, MS_FRAMEPLUS_OFFSET_xext) = 0.0;					//35
		PROFILE_VAL(prof_data, i, MS_FRAMEPLUS_OFFSET_fx) = 0.0;					//35
	}

	ms_profile->data = (scoretype *)malloc(sizeof(scoretype) * ms_profile->line_size * prof_length);

	{
		FILE *output;
		char strbuf[200];
		sprintf(strbuf, "%s_%0.5d", "soft_api.dbg.a.", 0);
		output = fopen(strbuf, "w");

		for (size_t i = 0; i < prof_length; i++) { /* for each profile line */

											/* First model->max_dprf profile lines are not actual lines but are added
											so that algorithm could go beyond the boundary and first lines would
											not be special cases. MS computes first lines separately. */
			data_from = &(prof_data[(i + 0/*0== SoftInfo_g.model->max_dprf*/) * FP_PROFILE_LINE_SIZE]);
			data_to = &(ms_profile->data[i*ms_profile->line_size]);

			for (size_t j = 0; j < line_size; j++) { /* for each position in the profile */
											  /*float mre = max_rel_error;/-* debug */
				data_to[j] = get_scaled(data_from[j], scale, &max_rel_error);
				fprintf(output, "i= %zd j= %zd %f * %f -> %d\n", i, j, (double)data_from[j], scale, data_to[j]);
				//						printf(": %d %d data=%.10f mre=%.10f\n",i,j,data_from[j],max_rel_error);
				//					if (mre != max_rel_error)
				//							printf("DEBUG: MS_gcprofile2msone_profile(): %d %d data=%.10f mre=%.10f\n",i,j,data_from[j],max_rel_error);
			} /* for(j=0;j<line_size;j++) */

		} /* for(i=0;i<prof_length;i++) */

		fclose(output);
	}
	free(prof_data);

	int *ms_sequence = (int *)malloc(sizeof(int) * dseq->len+5);
	ms_sequence[0] = 13;
	for (size_t r = 0; r < dseq->len; r++) {
		ms_sequence[r+1] = dseq->seq[r];
	}
	inpf_t max_v;
	best_t max_s = MS_Score_FP(ms_sequence, ms_profile, prof_length, dseq->len+1);
	max_v.score = ((float)(max_s.score)) * ms_profile->scaleback;
	max_v.x = max_s.x;
	max_v.y = max_s.y;
	max_v.state = -1; /* unknown */
	free(ms_sequence);
	free(ms_profile->data);
	free(ms_profile);
	return max_v.score;
}

double fp_gencore_matrix(const search_fp_profile_t * sp, const sequence_t * dseq, const sequence_t * qseq)
{
	/* model->profile_line_size = 66;  len of raw prof*/
	//MS_gcprofile_parts_bounds(model, &read_base, &line_size); gets own frame plus profile from model file
	int read_base = 20; /*shift in raw profile*/
	int line_size = 46; /*32 matrix + 14 gaps*/
	float *data_from;
	scoretype *data_to;
	float scale = sp->mtx->scale;
	float max_rel_error = 0;
	ms_profile_t * ms_profile = malloc(sizeof(ms_profile_t));
	ms_profile->scaleback = 1.0 / sp->mtx->scale;
	ms_profile->line_size = line_size;
	int prof_length = qseq->len + 1; // INT_MIN + seq(i-1)
	int prf_len = (prof_length)* FP_PROFILE_LINE_SIZE * sizeof(float);
	float *prof_data = (float *)malloc(prf_len);

	for (size_t j = 0; j < PROFILE_EMISSIONS_NO; j++)
		PROFILE_VAL(prof_data, 0/* *32+6=38*/, j) = INT_MIN; /*sigaev: global aligment*/

															 /* gap penalties copied from first line */
	prof_data[MS_FRAMEPLUS_OFFSET_matchmax5] = sp->matchMax5;	//32
	prof_data[MS_FRAMEPLUS_OFFSET_matchmax6] = sp->matchMax6;	//33
	prof_data[MS_FRAMEPLUS_OFFSET_xopen1] = sp->gapOpen;		//43
	prof_data[MS_FRAMEPLUS_OFFSET_xext1] = sp->gapExt;			//44
	prof_data[MS_FRAMEPLUS_OFFSET_fopen] = sp->gapFrameOpen;	//37
	prof_data[MS_FRAMEPLUS_OFFSET_fext] = sp->gapFrameExt;		//38
	prof_data[MS_FRAMEPLUS_OFFSET_yopen1] = sp->gapOpen;		//39
	prof_data[MS_FRAMEPLUS_OFFSET_yopen2] = sp->gapOpen2;		//40
	prof_data[MS_FRAMEPLUS_OFFSET_yopen3] = sp->gapOpen3;		//41
	prof_data[MS_FRAMEPLUS_OFFSET_yext] = sp->gapExt;			//42
	prof_data[MS_FRAMEPLUS_OFFSET_fx1] = sp->gapFrame;			//45 latest
																/* @@ next 3 fields are not used */
																//	MS_FRAMEPLUS_OFFSET_xopen, MS_FRAMEPLUS_OFFSET_xext, MS_FRAMEPLUS_OFFSET_fx,
	prof_data[MS_FRAMEPLUS_OFFSET_xopen] = 0.0;					//34
	prof_data[MS_FRAMEPLUS_OFFSET_xext] = 0.0;					//35
	prof_data[MS_FRAMEPLUS_OFFSET_fx] = 0.0;					//35


	for (size_t i = 1; i < prof_length; i++) {
		int k = qseq->seq[i - 1];
		for (size_t j = 0; j < PROFILE_EMISSIONS_NO; j++) { //#define PROFILE_EMISSIONS_NO 32
			double v;
			if (!sp->mtx)
				v = SCORE(k, j, 1.0, -1.0);
			else
				v = VDTABLE(k, j);
			//			printf(" %f, ", v);
			PROFILE_VAL(prof_data, i/* *32+6=38*/, j) = v; //create profile line from  32 CM("A") to 38 (36-39: 1,1,1,1,0,0)
		}
		// prof_data[(i)*PROFILE_LINE_SIZE + (j)] == prof_data[(i)*38 + (j)]

		PROFILE_VAL(prof_data, i, MS_FRAMEPLUS_OFFSET_matchmax5) = sp->matchMax5;	//32
		PROFILE_VAL(prof_data, i, MS_FRAMEPLUS_OFFSET_matchmax6) = sp->matchMax6;	//33
		PROFILE_VAL(prof_data, i, MS_FRAMEPLUS_OFFSET_xopen1) = sp->gapOpen;		//43
		PROFILE_VAL(prof_data, i, MS_FRAMEPLUS_OFFSET_xext1) = sp->gapExt;			//44
		PROFILE_VAL(prof_data, i, MS_FRAMEPLUS_OFFSET_fopen) = sp->gapFrameOpen;	//37
		PROFILE_VAL(prof_data, i, MS_FRAMEPLUS_OFFSET_fext) = sp->gapFrameExt;		//38
		PROFILE_VAL(prof_data, i, MS_FRAMEPLUS_OFFSET_yopen1) = sp->gapOpen;		//39
		PROFILE_VAL(prof_data, i, MS_FRAMEPLUS_OFFSET_yopen2) = sp->gapOpen2;		//40
		PROFILE_VAL(prof_data, i, MS_FRAMEPLUS_OFFSET_yopen3) = sp->gapOpen3;		//41
		PROFILE_VAL(prof_data, i, MS_FRAMEPLUS_OFFSET_yext) = sp->gapExt;			//42
		PROFILE_VAL(prof_data, i, MS_FRAMEPLUS_OFFSET_fx1) = sp->gapFrame;			//45 latest
																					/* @@ next 3 fields are not used */
																					//	MS_FRAMEPLUS_OFFSET_xopen, MS_FRAMEPLUS_OFFSET_xext, MS_FRAMEPLUS_OFFSET_fx,
		PROFILE_VAL(prof_data, i, MS_FRAMEPLUS_OFFSET_xopen) = 0.0;					//34
		PROFILE_VAL(prof_data, i, MS_FRAMEPLUS_OFFSET_xext) = 0.0;					//35
		PROFILE_VAL(prof_data, i, MS_FRAMEPLUS_OFFSET_fx) = 0.0;					//35
	}

	ms_profile->data = (scoretype *)malloc(sizeof(scoretype) * ms_profile->line_size * prof_length);

	{
		FILE *output;
		char strbuf[200];
		sprintf(strbuf, "%s_%0.5d", "soft_api.dbg.a.", 0);
		output = fopen(strbuf, "w");

		for (size_t i = 0; i < prof_length; i++) { /* for each profile line */

												   /* First model->max_dprf profile lines are not actual lines but are added
												   so that algorithm could go beyond the boundary and first lines would
												   not be special cases. MS computes first lines separately. */
			data_from = &(prof_data[(i + 0/*0== SoftInfo_g.model->max_dprf*/) * FP_PROFILE_LINE_SIZE]);
			data_to = &(ms_profile->data[i*ms_profile->line_size]);

			for (size_t j = 0; j < line_size; j++) { /* for each position in the profile */
													 /*float mre = max_rel_error;/-* debug */
				data_to[j] = get_scaled(data_from[j], scale, &max_rel_error);
				fprintf(output, "i= %zd j= %zd %f * %f -> %d\n", i, j, (double)data_from[j], scale, data_to[j]);
				//						printf(": %d %d data=%.10f mre=%.10f\n",i,j,data_from[j],max_rel_error);
				//					if (mre != max_rel_error)
				//							printf("DEBUG: MS_gcprofile2msone_profile(): %d %d data=%.10f mre=%.10f\n",i,j,data_from[j],max_rel_error);
			} /* for(j=0;j<line_size;j++) */

		} /* for(i=0;i<prof_length;i++) */

		fclose(output);
	}
	free(prof_data);

	int *ms_sequence = (int *)malloc(sizeof(int) * dseq->len + 5);
	ms_sequence[0] = 13;
	for (size_t r = 0; r < dseq->len; r++) {
		ms_sequence[r + 1] = dseq->seq[r];
	}
	inpf_t max_v;
	best_t max_s = score_frameplus_p2n(ms_sequence, ms_profile, prof_length, dseq->len + 1);
	max_v.score = ((float)(max_s.score)) * ms_profile->scaleback;
	max_v.x = max_s.x;
	max_v.y = max_s.y;
	max_v.state = -1; /* unknown */
	free(ms_sequence);
	free(ms_profile->data);
	free(ms_profile);
	return max_v.score;
}