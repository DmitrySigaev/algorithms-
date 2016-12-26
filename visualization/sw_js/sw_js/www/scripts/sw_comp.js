/* sw.js
 * Copyright 2016 Dmitry Sigaev
 *
 * Released under the MIT license -- see MIT-LICENSE for details
 */

/* Create empty matrix */
var Matrix = function (rows, cols) {
	var arr = [], row = [];
	while (cols--) row.push(0);
	while (rows--) arr.push(row.slice());
	return arr;
};


/*
  Match/Mismatch Scoring model
  This model assigns a score of +1 and -1 respectively if a match or a mismatch occurs,
  whereas a value equal to -d in case of gaps (insertions or deletions).
*/
var score = function (substitution, x, y) {
	if (x.toUpperCase() == y.toUpperCase()) {
		return substitution.match || 1;
	} else {
		return substitution.mismatch || -1;
	}
};

/*
  Substitution Scoring model.
*/
var ScoreT = {
	identityNuc: {
		A: { A: 2.0, B: -1.0, C: -1.0, D: -1.0, G: -1.0, H: -1.0, K: -1.0, M: -1.0, R: -1.0, S: -1.0, T: -1.0, U: -1.0, V: -1.0, W: -1.0, Y: -1.0 },
		B: { A: -1.0, B: 2.0, C: -1.0, D: -1.0, G: -1.0, H: -1.0, K: -1.0, M: -1.0, R: -1.0, S: -1.0, T: -1.0, U: -1.0, V: -1.0, W: -1.0, Y: -1.0 },
		C: { A: -1.0, B: -1.0, C: 2.0, D: -1.0, G: -1.0, H: -1.0, K: -1.0, M: -1.0, R: -1.0, S: -1.0, T: -1.0, U: -1.0, V: -1.0, W: -1.0, Y: -1.0 },
		D: { A: -1.0, B: -1.0, C: -1.0, D: 2.0, G: -1.0, H: -1.0, K: -1.0, M: -1.0, R: -1.0, S: -1.0, T: -1.0, U: -1.0, V: -1.0, W: -1.0, Y: -1.0 },
		G: { A: -1.0, B: -1.0, C: -1.0, D: -1.0, G: 2.0, H: -1.0, K: -1.0, M: -1.0, R: -1.0, S: -1.0, T: -1.0, U: -1.0, V: -1.0, W: -1.0, Y: -1.0 },
		H: { A: -1.0, B: -1.0, C: -1.0, D: -1.0, G: -1.0, H: 2.0, K: -1.0, M: -1.0, R: -1.0, S: -1.0, T: -1.0, U: -1.0, V: -1.0, W: -1.0, Y: -1.0 },
		K: { A: -1.0, B: -1.0, C: -1.0, D: -1.0, G: -1.0, H: -1.0, K: 2.0, M: -1.0, R: -1.0, S: -1.0, T: -1.0, U: -1.0, V: -1.0, W: -1.0, Y: -1.0 },
		M: { A: -1.0, B: -1.0, C: -1.0, D: -1.0, G: -1.0, H: -1.0, K: -1.0, M: 2.0, R: -1.0, S: -1.0, T: -1.0, U: -1.0, V: -1.0, W: -1.0, Y: -1.0 },
		R: { A: -1.0, B: -1.0, C: -1.0, D: -1.0, G: -1.0, H: -1.0, K: -1.0, M: -1.0, R: 2.0, S: -1.0, T: -1.0, U: -1.0, V: -1.0, W: -1.0, Y: -1.0 },
		S: { A: -1.0, B: -1.0, C: -1.0, D: -1.0, G: -1.0, H: -1.0, K: -1.0, M: -1.0, R: -1.0, S: 2.0, T: -1.0, U: -1.0, V: -1.0, W: -1.0, Y: -1.0 },
		T: { A: -1.0, B: -1.0, C: -1.0, D: -1.0, G: -1.0, H: -1.0, K: -1.0, M: -1.0, R: -1.0, S: -1.0, T: 2.0, U: 1.0, V: -1.0, W: -1.0, Y: -1.0 },
		U: { A: -1.0, B: -1.0, C: -1.0, D: -1.0, G: -1.0, H: -1.0, K: -1.0, M: -1.0, R: -1.0, S: -1.0, T: 1.0, U: 2.0, V: -1.0, W: -1.0, Y: -1.0 },
		V: { A: -1.0, B: -1.0, C: -1.0, D: -1.0, G: -1.0, H: -1.0, K: -1.0, M: -1.0, R: -1.0, S: -1.0, T: -1.0, U: -1.0, V: 2.0, W: -1.0, Y: -1.0 },
		W: { A: -1.0, B: -1.0, C: -1.0, D: -1.0, G: -1.0, H: -1.0, K: -1.0, M: -1.0, R: -1.0, S: -1.0, T: -1.0, U: -1.0, V: -1.0, W: 2.0, Y: -1.0 },
		Y: { A: -1.0, B: -1.0, C: -1.0, D: -1.0, G: -1.0, H: -1.0, K: -1.0, M: -1.0, R: -1.0, S: -1.0, T: -1.0, U: -1.0, V: -1.0, W: -1.0, Y: 2.0 }
	}
};

var S = function (substitution, x, y) {
	var substitutional_matrix = substitution.submatrix || ScoreT.identityNuc;
	return substitutional_matrix[x.toUpperCase()][y.toUpperCase()];
};

var is_match = function (x, y) {
	if (x.toUpperCase() == y.toUpperCase()) {
		return true;
	} else {
		return false;
	}
};


/*
Gap Penalty models:
1. y = -d(insertion/deletion) 
Constant Gap model, that assigns a penalty equal to a d
value to each gap found during the alignment and so capable to evaluate
only the presence of a gap event but not its extension; 

Linear Gap model
that considers instead the gap length (g) to score the alignment , giving
the possibility to evaluate with different scores gaps of different sizes;
*/

var sw_linear_gap = function (search_profile, s1, s2) {
	var gap = search_profile.gap || -1;
	var substitution = search_profile.S || { method: score, match: 1.0, mismatch: -1.0 };
	var l1 = s1.length;
	var l2 = s2.length;
	var score_mat = Matrix(l1, l2);
	var trace_mat = Matrix(l1, l2);

	for (i = 0; i < l1; ++i) {

		for (j = 0; j < l2; ++j) {

			/* This is the first row / column which is all zeros */
			if (i == 0 || j == 0) {
				score_mat[i][j] = 0;
				trace_mat[i][j] = 3;
				continue;

			} else {
				var d_last = score_mat[i - 1][j - 1];
				var u_last = score_mat[i - 1][j];
				var l_last = score_mat[i][j - 1];
			}
			var s = substitution.method(substitution, s1[i], s2[j])
			var d_new = d_last + s;
			var u_new = u_last + gap;
			var l_new = l_last + gap;
			score_mat[i][j] = Math.max(d_new, u_new, l_new, 0);
			var arr = [score_mat[i][j], u_new, l_new, 0];
			var trace = arr.indexOf(Math.max.apply(Math, arr));
			trace_mat[i][j] = trace;
		}
	}
	/* console.log(score_mat); */
	var mscore = score_mat[0][0];
	for (i = 0; i < l1; ++i) {
		for (j = 0; j < l2; ++j) {
			if (mscore < score_mat[i][j])
				mscore = score_mat[i][j];
		}
	}
	/* console.log(mscore); */

	return [mscore, score_mat, trace_mat];
}

/*
 * Sequence alignments with constant and linear gap penalties can be computed in time O(n*m) for two
 * sequences of lingth m and n. With affine gap penalties the time increases to O(n*m*(n+m)),
 * since for each cell the algorithm has to checkif a gap is extended or a new one is opened.
 * In 1982 Gotoh described a method to compute optimal sequence alignments, with affine gap
 * penalties, in time O(n*m). His version uses two additional matrices(E and F) to keep track of open gap.
 * E keeps track of gaps in the query sequence and F if gaps in the database sequence.
 * The values in E are calculated as the maximum of the previous value in H plus the costs Q for opening
 * a gap an the previous balue in E plus the costs T for extending a gap.  The values in F are
 * computed the same way, except for gaps in the other sequence. The values in H are computed like
 * in the Smith-Waterman with linear gap costs, except the value in E and F are used ad the gap costs.
 */

var sw_affine_gotoh_gap = function (search_profile /* {match/mismatch or submatrix, gapOpen, gapExt } */, dseq, qseq) {
	var gapOpen = search_profile.gapOpen || -1;
	var gapExt = search_profile.gapExt || 0;
	var substitution = search_profile.S || { method: score, match: 1.0, mismatch: -1.0 };
	var l1 = dseq.length;
	var l2 = qseq.length;
	var score_mat = Matrix(l1, l2);
	var trace_mat = Matrix(l1, l2);
	var ee = Matrix(l1, l2);
	//  var ff = Matrix(l1, l2);
	var ff = Matrix(l2, l1);

	for (i = 0; i < l1; ++i) {

		for (j = 0; j < l2; ++j) {

			// This is the first row / column which is all zeros
			if (i == 0 || j == 0) {
				if (i == 0)
					ee[i][j] = 0;
				if (j == 0)
					ff[j][i] = 0; // j->i
				score_mat[i][j] = 0;
				trace_mat[i][j] = 0;
				continue;

			} else {
				var e_last = ee[i][j - 1];
				var f_last = ff[j][i - 1];
				var d_last = score_mat[i - 1][j - 1];
				var u_last = score_mat[i - 1][j];
				var l_last = score_mat[i][j - 1];
			}
			var s = substitution.method(substitution, dseq[i], qseq[j]);
			var d_new = d_last + s
			var u_new = u_last + gapOpen;
			var l_new = l_last + gapOpen;
			var e_new = e_last + gapExt;
			var f_new = f_last + gapExt;

			ee[i][j] = Math.max(e_new, l_new);
			ff[j][i] = Math.max(f_new, u_new);

			score_mat[i][j] = Math.max(d_new, ee[i][j], ff[j][i], 0);
			if (score_mat[i][j] == d_new && !is_match(dseq[i], qseq[j])) {
				trace_mat[i][j] |= (1 << 0);// #define LAL_MASK_MISMATCH      (1<<0)
				//				ee[i][j] = ee[i][j] + s;
				//				ff[j][i] = ff[j][i] + s;

			}
			if (score_mat[i][j] == ff[j][i])
				trace_mat[i][j] |= (1 << 2); // #define LAL_MASK_GAP_OPEN_LEFT (1<<3)
			if (score_mat[i][j] == ee[i][j])
				trace_mat[i][j] |= (1 << 3);// #define LAL_MASK_GAP_OPEN_UP   (1<<2)
			if (score_mat[i][j] == d_new && is_match(dseq[i], qseq[j]))
				trace_mat[i][j] |= (1 << 1);// #define LAL_MASK_MATCH         (1<<1)
			if (score_mat[i][j] == 0)
				trace_mat[i][j] = (1 << 6); // #define LAL_MASK_ZERO          (1<<6)
			/*
			var arr = [d_new, ff[i][j], ee[i][j], 0];
			var trace = arr.indexOf(Math.max.apply(Math, arr));
			trace_mat[i][j] = trace;
			*/
		}
	}
	/* console.log(score_mat); */
	var mscore = score_mat[0][0];
	for (i = 0; i < l1; ++i) {
		for (j = 0; j < l2; ++j) {
			if (mscore < score_mat[i][j])
				mscore = score_mat[i][j];
		}
	}
	/* console.log(mscore); */

	return [mscore, score_mat, trace_mat];
}

var sw_affine_gap_gotoh_v1 = function (search_profile /* {match/mismatch or submatrix, gapOpen, gapExt } */, dseq, qseq) {
	var gapOpen = search_profile.gapOpen || -1;
	var gapExt = search_profile.gapExt || 0;
	var substitution = search_profile.S || { method: score, match: 1.0, mismatch: -1.0 };
	var l1 = dseq.length;
	var l2 = qseq.length;
	var score_mat = Matrix(l1, l2);
	var trace_mat = Matrix(l1, l2);
	var ee = Matrix(l1, l2);
	//  var ff = Matrix(l1, l2);
	var ff = Matrix(l2, l1);

	for (i = 0; i < l1; ++i) {

		for (j = 0; j < l2; ++j) {

			// This is the first row / column which is all zeros
			if (i == 0 || j == 0) {
				if (i == 0)
					ee[i][j] = 0;
				if (j == 0)
					ff[j][i] = 0; // j->i
				score_mat[i][j] = 0;
				trace_mat[i][j] = 0;
				continue;

			} else {
				var e_last = ee[i][j - 1];
				var f_last = ff[j][i - 1];
				var d_last = score_mat[i - 1][j - 1];
				var u_last = score_mat[i - 1][j];
				var l_last = score_mat[i][j - 1];
			}
			var s = substitution.method(substitution, dseq[i], qseq[j]);
			var d_new = d_last + s
			var u_new = u_last + gapOpen + gapExt;
			var l_new = l_last + gapOpen + gapExt;
			var e_new = e_last + gapExt;
			var f_new = f_last + gapExt;

			ee[i][j] = Math.max(e_new, l_new);
			ff[j][i] = Math.max(f_new, u_new);

			score_mat[i][j] = Math.max(d_new, ee[i][j], ff[j][i], 0);
			if (score_mat[i][j] == d_new && !is_match(dseq[i], qseq[j])) {
				trace_mat[i][j] |= (1 << 0);// #define LAL_MASK_MISMATCH      (1<<0)
				ee[i][j] = ee[i][j] + s; /* compensation */
				ff[j][i] = ff[j][i] + s;

			}
			if (score_mat[i][j] == ff[j][i])
				trace_mat[i][j] |= (1 << 2); // #define LAL_MASK_GAP_OPEN_LEFT (1<<3)
			if (score_mat[i][j] == ee[i][j])
				trace_mat[i][j] |= (1 << 3);// #define LAL_MASK_GAP_OPEN_UP   (1<<2)
			if (score_mat[i][j] == d_new && is_match(dseq[i], qseq[j]))
				trace_mat[i][j] |= (1 << 1);// #define LAL_MASK_MATCH         (1<<1)
			if (score_mat[i][j] == 0)
				trace_mat[i][j] = (1 << 6); // #define LAL_MASK_ZERO          (1<<6)
			/*
			var arr = [d_new, ff[i][j], ee[i][j], 0];
			var trace = arr.indexOf(Math.max.apply(Math, arr));
			trace_mat[i][j] = trace;
			*/
		}
	}
	/* console.log(score_mat); */
	var mscore = score_mat[0][0];
	for (i = 0; i < l1; ++i) {
		for (j = 0; j < l2; ++j) {
			if (mscore < score_mat[i][j])
				mscore = score_mat[i][j];
		}
	}
	/* console.log(mscore); */

	return [mscore, score_mat, trace_mat];
}

var sw_affine_gap_gotoh_v2 = function (search_profile /* {match/mismatch or submatrix, gapOpen, gapExt } */, dseq, qseq) {
	var gapOpen = search_profile.gapOpen || -1;
	var gapExt = search_profile.gapExt || 0;
	var substitution = search_profile.S || { method: score, match: 1.0, mismatch: -1.0 };
	var l1 = dseq.length;
	var l2 = qseq.length;
	var score_mat = Matrix(l1, l2);
	var trace_mat = Matrix(l1, l2);
	var ee = Matrix(l1, l2);
	//  var ff = Matrix(l1, l2);
	var ff = Matrix(l2, l1);

	for (i = 0; i < l1; ++i) {

		for (j = 0; j < l2; ++j) {

			// This is the first row / column which is all zeros
			if (i == 0 || j == 0) {
				if (i == 0)
					ee[i][j] = 0;
				if (j == 0)
					ff[j][i] = 0; // j->i
				score_mat[i][j] = 0;
				trace_mat[i][j] = 0;
				continue;

			} else {
				var e_last = ee[i][j - 1];
				var f_last = ff[j][i - 1];
				var d_last = score_mat[i - 1][j - 1];
				var u_last = score_mat[i - 1][j];
				var l_last = score_mat[i][j - 1];
			}
			var s = substitution.method(substitution, dseq[i], qseq[j]);
			var d_new = d_last + s
			var u_new = u_last + gapOpen + gapExt;
			var l_new = l_last + gapOpen + gapExt;
			var e_new = e_last + gapExt;
			var f_new = f_last + gapExt;

			ee[i][j] = Math.max(e_new, l_new);
			ff[j][i] = Math.max(f_new, u_new);

			score_mat[i][j] = Math.max(d_new, ee[i][j], ff[j][i], 0);
			if (score_mat[i][j] == d_new && !is_match(dseq[i], qseq[j])) {
				trace_mat[i][j] |= (1 << 0);// #define LAL_MASK_MISMATCH      (1<<0)
				//				ee[i][j] = ee[i][j] + s;
				//				ff[j][i] = ff[j][i] + s;

			}
			if (score_mat[i][j] == ff[j][i])
				trace_mat[i][j] |= (1 << 2); // #define LAL_MASK_GAP_OPEN_LEFT (1<<3)
			if (score_mat[i][j] == ee[i][j])
				trace_mat[i][j] |= (1 << 3);// #define LAL_MASK_GAP_OPEN_UP   (1<<2)
			if (score_mat[i][j] == d_new && is_match(dseq[i], qseq[j]))
				trace_mat[i][j] |= (1 << 1);// #define LAL_MASK_MATCH         (1<<1)
			if (score_mat[i][j] == 0)
				trace_mat[i][j] = (1 << 6); // #define LAL_MASK_ZERO          (1<<6)
			/*
			var arr = [d_new, ff[i][j], ee[i][j], 0];
			var trace = arr.indexOf(Math.max.apply(Math, arr));
			trace_mat[i][j] = trace;
			*/
		}
	}
	/* console.log(score_mat); */
	var mscore = score_mat[0][0];
	for (i = 0; i < l1; ++i) {
		for (j = 0; j < l2; ++j) {
			if (mscore < score_mat[i][j])
				mscore = score_mat[i][j];
		}
	}
	/* console.log(mscore); */

	return [mscore, score_mat, trace_mat];
}

/*
   Other the Smith-Waterman algorithm implementation, which is described in:
   http://pages.cs.wisc.edu/~bsettles/ibs08/lectures/02-alignment.pdf
*/
var sw_affine_gap_v1 = function (search_profile, dseq, qseq) {
	var gapOpen = search_profile.gapOpen || -1;
	var gapExt = search_profile.gapExt || 0;
	var substitution = search_profile.S || { method: score, match: 1.0, mismatch: -1.0 };
	var l1 = dseq.length;
	var l2 = qseq.length;
	var score_mat = Matrix(l1, l2);
	var trace_mat = Matrix(l1, l2);
	var ee = Matrix(l1, l2);
	var ff = Matrix(l1, l2);

	for (i = 0; i < l1; ++i) {

		for (j = 0; j < l2; ++j) {

			/* initialization: http://pages.cs.wisc.edu/~bsettles/ibs08/lectures/02-alignment.pdf  ee(Iy) and ff(Ix) = -infinity 
               but in accordance with http://iwbbio.ugr.es/2014/papers/IWBBIO_2014_paper_143.pdf */
			if (i == 0 || j == 0) { /**/
				if (i == 0)
					ee[i][j] = 0;//-1000;
				if (j == 0)
					ff[i][j] = 0;//-1000;
				score_mat[i][j] = 0;
				trace_mat[i][j] = 0;
				continue;

			} else {
				var ey_last = ee[i][j - 1];
				var fx_last = ff[i - 1][j];
				var m_last = score_mat[i - 1][j - 1];
				var mx_last = score_mat[i - 1][j];
				var my_last = score_mat[i][j - 1];
			}
			var s = substitution.method(substitution, dseq[i], qseq[j]);
			var m_new = m_last + s;
			var mx_new = mx_last + gapOpen;
			var my_new = my_last + gapOpen;
			ee[i][j] = Math.max(ey_last + gapExt, my_new);
			ff[i][j] = Math.max(fx_last + gapExt, mx_new);

			var y = ee[i - 1][j - 1] + s;
			var x = ff[i - 1][j - 1] + s;
			score_mat[i][j] = Math.max(m_new, y/*y*//*ee[i - 1][j - 1] + s*/, x  /*x*//*ff[i - 1][j - 1] + s*/, 0);

			if (score_mat[i][j] == m_new && is_match(dseq[i], qseq[j]))
				trace_mat[i][j] |= (1 << 1);// #define LAL_MASK_MATCH         (1<<1)
			if (score_mat[i][j] == m_new && !is_match(dseq[i], qseq[j]))
				trace_mat[i][j] |= (1 << 0);// #define LAL_MASK_MISMATCH      (1<<0)
			if (score_mat[i][j] == x)
				trace_mat[i][j] |= (1 << 2); // #define LAL_MASK_GAP_OPEN_LEFT (1<<3)
			if (score_mat[i][j] == y)
				trace_mat[i][j] |= (1 << 3);// #define LAL_MASK_GAP_OPEN_UP   (1<<2)
			if (score_mat[i][j] == 0)
				trace_mat[i][j] |= (1 << 6); // #define LAL_MASK_ZERO          (1<<6)
			/*
            var arr = [m_new, ff[i - 1][j - 1] + s, ee[i - 1][j - 1] + s, 0];

            var trace = arr.indexOf(Math.max.apply(Math, arr));
            trace_mat[i][j] = trace;
            */
		}
	}
	/* console.log(score_mat); */
	var mscore = score_mat[0][0];
	for (i = 0; i < l1; ++i) {
		for (j = 0; j < l2; ++j) {
			if (mscore < score_mat[i][j])
				mscore = score_mat[i][j];
		}
	}
	/* console.log(mscore); */

	return [mscore, score_mat, trace_mat];
}


var sw_affine_gap_v1_comp = function (search_profile, dseq, qseq) {
	var gapOpen = search_profile.gapOpen || -1;
	var gapExt = search_profile.gapExt || 0;
	var substitution = search_profile.S || { method: score, match: 1.0, mismatch: -1.0 };
	var l1 = dseq.length;
	var l2 = qseq.length;
	var score_mat = Matrix(l1, l2);
	var trace_mat = Matrix(l1, l2);
	var ee = Matrix(l1, l2);
	var ff = Matrix(l1, l2);

	for (i = 0; i < l1; ++i) {

		for (j = 0; j < l2; ++j) {

			/* initialization: http://pages.cs.wisc.edu/~bsettles/ibs08/lectures/02-alignment.pdf  ee(Iy) and ff(Ix) = -infinity 
               but in accordance with http://iwbbio.ugr.es/2014/papers/IWBBIO_2014_paper_143.pdf */
			if (i == 0 || j == 0) { /**/
				if (i == 0)
					ee[i][j] = 0;//-1000;
				if (j == 0)
					ff[i][j] = 0;//-1000;
				score_mat[i][j] = 0;
				trace_mat[i][j] = 0;
				continue;

			} else {
				var ey_last = ee[i][j - 1];
				var fx_last = ff[i - 1][j];
				var m_last = score_mat[i - 1][j - 1];
				var mx_last = score_mat[i - 1][j];
				var my_last = score_mat[i][j - 1];
			}
			var s = substitution.method(substitution, dseq[i], qseq[j]);
			var m_new = m_last + s;
			var mx_new = mx_last + gapOpen;
			var my_new = my_last + gapOpen;
			ee[i][j] = Math.max(ey_last + gapExt, my_new);
			ff[i][j] = Math.max(fx_last + gapExt, mx_new);

			var y = ee[i - 1][j - 1] + s;
			var x = ff[i - 1][j - 1] + s;
			score_mat[i][j] = Math.max(m_new, y/*y*//*ee[i - 1][j - 1] + s*/, x  /*x*//*ff[i - 1][j - 1] + s*/, 0);

			if (score_mat[i][j] == m_new && is_match(dseq[i], qseq[j]))
				trace_mat[i][j] |= (1 << 1);// #define LAL_MASK_MATCH         (1<<1)
			if (score_mat[i][j] == m_new && !is_match(dseq[i], qseq[j]))
				trace_mat[i][j] |= (1 << 0);// #define LAL_MASK_MISMATCH      (1<<0)
			if (score_mat[i][j] == x)
				trace_mat[i][j] |= (1 << 2); // #define LAL_MASK_GAP_OPEN_LEFT (1<<3)
			if (score_mat[i][j] == y)
				trace_mat[i][j] |= (1 << 3);// #define LAL_MASK_GAP_OPEN_UP   (1<<2)
			if (score_mat[i][j] == 0)
				trace_mat[i][j] |= (1 << 6); // #define LAL_MASK_ZERO          (1<<6)
			/*
            var arr = [m_new, ff[i - 1][j - 1] + s, ee[i - 1][j - 1] + s, 0];

            var trace = arr.indexOf(Math.max.apply(Math, arr));
            trace_mat[i][j] = trace;
            */
		}
	}
	/* console.log(score_mat); */
	var mscore = score_mat[0][0];
	for (i = 0; i < l1; ++i) {
		for (j = 0; j < l2; ++j) {
			if (mscore < score_mat[i][j])
				mscore = score_mat[i][j];
		}
	}
	/* console.log(mscore); */

	return [mscore, score_mat, trace_mat, ee, ff];
}

/*
   Yet another the Smith-Waterman algorithm implementation, which is described in:
   http://iwbbio.ugr.es/2014/papers/IWBBIO_2014_paper_143.pdf
*/
var sw_affine_gap_v2 = function (search_profile /* {match/mismatch or submatrix, gapOpen, gapExt } */, dseq, qseq) {
	var gapOpen = search_profile.gapOpen || -1;
	var gapExt = search_profile.gapExt || 0;
	var substitution = search_profile.S || { method: score, match: 1.0, mismatch: -1.0 };
	var l1 = dseq.length;
	var l2 = qseq.length;
	var score_mat = Matrix(l1, l2);
	var trace_mat = Matrix(l1, l2);
	var ee = Matrix(l1, l2);
	var ff = Matrix(l1, l2);

	for (i = 0; i < l1; ++i) {

		for (j = 0; j < l2; ++j) {

			/* initialization in accordance with http://iwbbio.ugr.es/2014/papers/IWBBIO_2014_paper_143.pdf */
			if (i == 0 || j == 0) { /**/
				if (i == 0)
					ee[i][j] = 0;
				if (j == 0)
					ff[i][j] = 0;
				score_mat[i][j] = 0;
				trace_mat[i][j] = 3;
				continue;

			} else {
				var ey_last = ee[i][j - 1];
				var fx_last = ff[i - 1][j];
				var m_last = score_mat[i - 1][j - 1];
				var mx_last = score_mat[i - 1][j];
				var my_last = score_mat[i][j - 1];
			}
			var s = substitution.method(substitution, dseq[i], qseq[j]);  /* Substitutional Matrix */
			var m_new = m_last + s;
			var mx_new = mx_last + gapOpen; /* */
			var my_new = my_last + gapOpen;
			ee[i][j] = Math.max(ey_last + gapExt, my_new, 0);
			ff[i][j] = Math.max(fx_last + gapExt, mx_new, 0);

			score_mat[i][j] = Math.max(m_new, /*y*/ee[i - 1][j - 1] + s, /*x*/ff[i - 1][j - 1] + s, 0);
			var arr = [m_new, ff[i - 1][j - 1] + s, ee[i - 1][j - 1] + s, 0];
			var trace = arr.indexOf(Math.max.apply(Math, arr));
			trace_mat[i][j] = trace;
		}
	}
	/* console.log(score_mat); */
	var mscore = score_mat[0][0];
	for (i = 0; i < l1; ++i) {
		for (j = 0; j < l2; ++j) {
			if (mscore < score_mat[i][j])
				mscore = score_mat[i][j];
		}
	}
	/* console.log(mscore); */

	return [mscore, score_mat, trace_mat];
}


var sw_affine_gap_genc = function (search_profile /* {match/mismatch or submatrix, gapOpen, gapExt } */, dseq, qseq) {
	var gapOpen = search_profile.gapOpen || -1;
	var gapExt = search_profile.gapExt || 0;
	var substitution = search_profile.S || { method: score, match: 1.0, mismatch: -1.0 };
	var lx = dseq.length;
	var ly = qseq.length;
	var score_mat = Matrix(ly, lx);
	var trace_mat = Matrix(ly, lx);
	var score_mat2 = Matrix(lx, ly);
	var score_mat3 = Matrix(lx, ly);
	var trace_mat2 = Matrix(lx, ly);
	var yy = Matrix(ly, lx);
	var xx = Matrix(ly, lx);

	var prev_score = [];
	var prev_yskip = [];
	var max_score_perv = 0.0;
	var max_score_perv2 = 0.0;
	var d_xskipmatch = -100.0;
	var d_yskipmatch = -100.0;

	for (x = 0; x < lx; ++x) {
		prev_score[x] = 0;
		prev_yskip[x] = -100.0;
	}
	var d_quality;

	for (y = 0; y < ly; ++y) {
		for (x = 0; x < lx; ++x) {
			if (x == 0 || y == 0) { /**/
				if (x == 0) {
					xx[y][x] = -100.0;  /*xskipmatch = MS_SCORE_MININF;*/
				}
				if (y == 0)
					yy[y][x] = -100.0; // prev_yskip = [];
				score_mat[y][x] = 0; // prev_score[x];
				trace_mat[y][x] = 0;
				trace_mat2[y][x] = 0;
				continue;

			} else {
				var m_last = score_mat[y - 1][x - 1]; // prev_score[x];???

				var ex_last = xx[y][x - 1];
				var fy_last = yy[y - 1][x];
				var my_last = score_mat[y - 1][x];
				var mx_last = score_mat[y][x - 1];
			}
			var s = substitution.method(substitution, dseq[x], qseq[y]);  /* Substitutional Matrix */

			if (x == 1) {
				var d_quality = s;
				var d_lastquality = d_quality;
				max_score_perv = (s > max_score_perv2) ? s : max_score_perv; /*checkbest_m(v, 1, y + 1, max_v);*/
				d_xskipmatch = -10000.0; /*xskipmatch = MS_SCORE_MININF;*/
			}
			else {

				d_score = s;               // score = prof_line_p[pre->nseq];
				d_quality = prev_score[x]; // quality = pre->prequal;
				d_yskipmatch = prev_yskip[x]; // yskipmatch = pre->yskipmatch;

				if (d_yskipmatch <= 0) {
					if (d_xskipmatch <= 0) {
						/* yskipmatch <= 0 and xskipmatch <= 0 */
						if (d_quality <= 0) {
							d_quality = 0;
						} else {
							if (d_quality > gapOpen/*pentest*/) {
								d_xskipmatch = d_quality + gapOpen;
								prev_yskip[x] = d_quality + gapOpen /* gapOpenY == penopen1 */;
							}
							max_score_perv = (d_quality > max_score_perv) ? d_quality : max_score_perv; /* checkbest_m(quality, pre - prev_line, y, max_v) */
						}
					} else {
						/* yskipmatch <= 0 && xskipmatch > 0 */
						// 1 - x
						if (d_quality <= 0) {
							d_quality = d_xskipmatch;
							d_xskipmatch += gapExt;
						} else {
							var d_bestjumpY = d_quality + gapOpen; /*penopen1*/
							prev_yskip[x] = d_bestjumpY;
							if (d_quality < d_xskipmatch) { /* q<0 handled here */
								d_quality = d_xskipmatch;
								d_xskipmatch += gapExt;  /* assuming penopen<penext @@ test input */
							} else {
								var bestjumpX = d_quality + gapOpen;
								d_xskipmatch += gapExt;
								if (bestjumpX > d_xskipmatch) d_xskipmatch = bestjumpX;
								max_score_perv = (d_quality > max_score_perv) ? d_quality : max_score_perv; /* checkbest_m(quality, pre - prev_line, y, max_v) */
							}
						}

					}
				} else {
					/* yskip > 0 */
					if (d_xskipmatch <= 0) {
						/* yest > 0 and xskipmatch <= 0 */
						// 1 - y 
						if (d_quality <= 0) {
							d_quality = d_yskipmatch;
							d_yskipmatch += gapExt; // sp->gapExt1;
						}
						else {
							var bestjumpX = d_quality + gapOpen;
							/*if (bestjump > 0)*/ d_xskipmatch = bestjumpX; /*becames > 0 or steal < 0 in case d_quality < sp->gapOpen */

							if (d_quality < d_yskipmatch) { /* q<0 handled here */
								d_quality = d_yskipmatch;
								d_yskipmatch += gapExt; // penext1; /* assuming penopen1<penext1 */
							}
							else {
								var bestjumpY = d_quality + gapOpen; // penopen1;
								d_yskipmatch += gapExt;
								if (bestjumpY > d_yskipmatch) d_yskipmatch = bestjumpY;
								max_score_perv = (d_quality > max_score_perv) ? d_quality : max_score_perv; /* checkbest_m(quality, pre - prev_line, y, max_v) */
							}
						}
					} else {
						/* xskipmatch > 0 && yskipmatch > 0 */
						if (d_quality < d_xskipmatch) {
							if (d_quality < d_yskipmatch) { /* q<0 handled here */
								if (d_yskipmatch > d_xskipmatch)
									d_quality = d_yskipmatch;
								else
									d_quality = d_xskipmatch;
								d_yskipmatch += gapExt; //penext1; /* assuming penopen1<penext1 */
							}
							else {
								var bestjumpY = d_quality + gapOpen;//penopen1;
								d_yskipmatch += gapExt;// penext1;
								if (bestjumpY > d_yskipmatch) d_yskipmatch = bestjumpY;
								d_quality = d_xskipmatch;
							}
							d_xskipmatch += gapExt;// penext;  /* assuming penopen<penext */
						} else {
							/* quality>=xskipmatch */
							var bestjumpX = d_quality + gapOpen;
							d_xskipmatch += gapExt;
							if (bestjumpX > d_xskipmatch) d_xskipmatch = bestjumpX;
							if (d_quality < d_yskipmatch) {
								d_quality = d_yskipmatch;
								d_yskipmatch += gapExt; // penext1; /* assuming penopen1<penext1 */
							} else {
								var bestjumpY = d_quality + gapOpen;// penopen1;
								d_yskipmatch += gapExt;// penext1;
								if (bestjumpY > d_yskipmatch) d_yskipmatch = bestjumpY;
								max_score_perv = (d_quality > max_score_perv) ? d_quality : max_score_perv; /* checkbest_m(quality, pre - prev_line, y, max_v) */
							}
						} /* quality>=xskipmatch */
					} /* yskipmatch test */
					prev_yskip[x] = d_yskipmatch;
				} /* xskipmatch test */

				d_quality += d_score;
				prev_score[x] = d_lastquality;
				d_lastquality = d_quality;
			}
			var m_new = m_last + s;
			var mx_new = mx_last + gapOpen; /* */
			var my_new = my_last + gapOpen;
			xx[y][x] = Math.max(ex_last + gapExt, mx_new, 0);
			yy[y][x] = Math.max(fy_last + gapExt, my_new, 0);
			//		ee[i][j] = Math.max(e_new, l_new);
			//		ff[j][i] = Math.max(f_new, u_new);

			score_mat[y][x] = Math.max(m_new, xx[y][x], yy[y][x], 0);
			if (score_mat[y][x] == m_new && !is_match(dseq[x], qseq[y])) {
				trace_mat[y][x] |= (1 << 0);// #define LAL_MASK_MISMATCH      (1<<0)
			}
			if (score_mat[y][x] == xx[y][x])
				trace_mat[y][x] |= (1 << 2); // #define LAL_MASK_GAP_OPEN_LEFT (1<<3)
			if (score_mat[y][x] == yy[y][x])
				trace_mat[y][x] |= (1 << 3);// #define LAL_MASK_GAP_OPEN_UP   (1<<2)
			if (score_mat[y][x] == m_new && is_match(dseq[x], qseq[y]))
				trace_mat[y][x] |= (1 << 1);// #define LAL_MASK_MATCH         (1<<1)
			if (score_mat[y][x] == 0)
				trace_mat[y][x] = (1 << 6); // #define LAL_MASK_ZERO          (1<<6)
			/*
			var arr = [d_new, ff[i][j], ee[i][j], 0];
			var trace = arr.indexOf(Math.max.apply(Math, arr));
			trace_mat[i][j] = trace;
			*/

			var test = Math.max(d_quality, d_yskipmatch, d_xskipmatch, 0);
			score_mat3[x][y] = test;
			if (test == d_quality && !is_match(dseq[x], qseq[y])) {
				trace_mat2[x][y] |= (1 << 0);// #define LAL_MASK_MISMATCH      (1<<0)
			}
			if (test == d_xskipmatch)
				trace_mat2[x][y] |= (1 << 2); // #define LAL_MASK_GAP_OPEN_LEFT (1<<3)
			if (test == d_yskipmatch)
				trace_mat2[x][y] |= (1 << 3);// #define LAL_MASK_GAP_OPEN_UP   (1<<2)
			if (test == d_quality && is_match(dseq[x], qseq[y]))
				trace_mat2[x][y] |= (1 << 1);// #define LAL_MASK_MATCH         (1<<1)
			if (test == 0)
				trace_mat2[x][y] = (1 << 6); // #define LAL_MASK_ZERO          (1<<6)

		}
		if (y) {
			max_score_perv = (d_quality > max_score_perv) ? d_quality : max_score_perv; /* checkbest_m(quality, seqLen, y+1, max_v) */
			for (x = 0; x < lx - 1; ++x)
				score_mat2[x][y] = prev_score[x + 1];
			score_mat2[lx - 1][y] = d_quality;

		}
	}

	for (x = 0; x < lx; ++x) {
		max_score_perv = (prev_score[x] > max_score_perv) ? prev_score[x] : max_score_perv; /* checkbest_m(quality, x, yprofLen, max_v) */
		if (x != (lx - 1))
			score_mat2[x][y] = prev_score[x + 1];
		else
			score_mat2[lx - 1][y] = max_score_perv;

	}

	/* console.log(score_mat); */
	var mscore = score_mat[0][0];
	for (i = 0; i < ly; ++i) {
		for (j = 0; j < lx; ++j) {
			if (mscore < score_mat[i][j])
				mscore = score_mat[i][j];
		}
	}
	/* console.log(mscore); */

	//	return [mscore, score_mat, trace_mat];
	return [max_score_perv, score_mat2, trace_mat2];
}

/*
   Other the Smith-Waterman algorithm implementation, which is described in:
   http://pages.cs.wisc.edu/~bsettles/ibs08/lectures/02-alignment.pdf
*/
var sw_affine_gap_sg_v1 = function (search_profile, dseq, qseq) {
	var gapOpen = search_profile.gapOpen || -1;
	var gapExt = search_profile.gapExt || 0;
	var substitution = search_profile.S || { method: score, match: 1.0, mismatch: -1.0 };
	var l1 = dseq.length;
	var l2 = qseq.length;
	var h = Matrix(l1, l2);
	var trace_mat = Matrix(l1, l2);
	var ee = Matrix(l1, l2);
	var ff = Matrix(l1, l2);
	var er;
	var er0 = [];
	var er1 = [];
	var fr0 = [];
	var fr1 = [];
	var hr0 = [];
	var hr1 = [];

	for (var j = 0; j < l2; ++j) {
		er1[j] = 0;
		fr1[j] = 0;
		hr1[j] = 0;
	}

	for (var i = 0; i < l1; ++i) {
		fr1[0] = 0;
		er = 0;
		for (var j = 0; j < l2; ++j) {
			/* initialization: http://pages.cs.wisc.edu/~bsettles/ibs08/lectures/02-alignment.pdf  ee(Iy) and ff(Ix) = -infinity 
               but in accordance with http://iwbbio.ugr.es/2014/papers/IWBBIO_2014_paper_143.pdf */
			if (i == 0 || j == 0) { /**/
				if (i == 0) {
					ee[i][j] = 0;//-1000;
				}
				if (j == 0)
					ff[i][j] = 0;//-1000;
				h[i][j] = 0;
				trace_mat[i][j] = 0;
				hr = 0;
				continue;
			}
			var s = substitution.method(substitution, dseq[i], qseq[j]);
			var mx_new = hr0[j] + gapOpen;
			var my_new = hr + gapOpen;

			hr1[j] = hr = h[i][j] = Math.max(hr0[j - 1] + s, er0[j - 1] + s, fr0[j - 1] + s, 0);

			fr1[j] = ff[i][j] = Math.max(fr0[j] + gapExt, mx_new);
			er1[j - 1] = er;
			er = ee[i][j] = Math.max(er + gapExt, my_new);

			if (h == (hr0[j - 1] + s) && is_match(dseq[i], qseq[j]))
				trace_mat[i][j] |= (1 << 1);// #define LAL_MASK_MATCH         (1<<1)
			if (h == (hr0[j - 1] + s) && !is_match(dseq[i], qseq[j]))
				trace_mat[i][j] |= (1 << 0);// #define LAL_MASK_MISMATCH      (1<<0)
			if (h == fr0[j - 1] + s)
				trace_mat[i][j] |= (1 << 2); // #define LAL_MASK_GAP_OPEN_LEFT (1<<3)
			if (h == er0[j - 1] + s)
				trace_mat[i][j] |= (1 << 3);// #define LAL_MASK_GAP_OPEN_UP   (1<<2)
			if (h == 0)
				trace_mat[i][j] |= (1 << 6); // #define LAL_MASK_ZERO          (1<<6)
			/*
            var arr = [m_new, ff[i - 1][j - 1] + s, ee[i - 1][j - 1] + s, 0];

            var trace = arr.indexOf(Math.max.apply(Math, arr));
            trace_mat[i][j] = trace;
            */
		}
		er0 = er1;
		for (var j = 0; j < l2; ++j) {
			fr0[j] = fr1[j];
		}
		for (var j = 0; j < l2; ++j) {
			hr0[j] = hr1[j];
		}
		//		fr0 = fr1;
	}
	/* console.log(score_mat); */
	var mscore = h[0][0];
	for (i = 0; i < l1; ++i) {
		for (j = 0; j < l2; ++j) {
			if (mscore < h[i][j])
				mscore = h[i][j];
		}
	}
	/* console.log(mscore); */

	return [mscore, h, trace_mat];
}

/*
   Other the Smith-Waterman algorithm implementation, which is described in:
   http://pages.cs.wisc.edu/~bsettles/ibs08/lectures/02-alignment.pdf
*/
var sw_affine_gap_sg = function (search_profile, dseq, qseq) {
	var gapOpen = search_profile.gapOpen || -1;
	var gapExt = search_profile.gapExt || 0;
	var substitution = search_profile.S || { method: score, match: 1.0, mismatch: -1.0 };
	var l1 = dseq.length;
	var l2 = qseq.length;
	var h = Matrix(l1, l2);
	var trace_mat = Matrix(l1, l2);
	var ee = Matrix(l1, l2);
	var ff = Matrix(l1, l2);
	var er;
	var er0 = [];
	var er1 = [];
	var fr0 = [];
	var fr1 = [];
	var hr0 = [];
	var hr1 = [];

	for (var j = 0; j < l2; ++j) {
		er1[j] = 0;
		fr1[j] = 0;
		hr1[j] = 0;
	}

	for (var i = 0; i < l1; ++i) {
		fr1[0] = 0;
		er = 0;
		for (var j = 0; j < l2; ++j) {
			/* initialization: http://pages.cs.wisc.edu/~bsettles/ibs08/lectures/02-alignment.pdf  ee(Iy) and ff(Ix) = -infinity 
               but in accordance with http://iwbbio.ugr.es/2014/papers/IWBBIO_2014_paper_143.pdf */
			if (i == 0 || j == 0) { /**/
				if (i == 0) {
					ee[i][j] = 0;//-1000;
				}
				if (j == 0)
					ff[i][j] = 0;//-1000;
				h[i][j] = 0;
				trace_mat[i][j] = 0;
				hr = 0;
				continue;
			}
			var s = substitution.method(substitution, dseq[i], qseq[j]);
			var mx_new = hr0[j] + gapOpen;
			var my_new = hr + gapOpen;

			hr1[j] = hr = h[i][j] = Math.max(hr0[j - 1] + s, er0[j - 1] + s, fr0[j - 1] + s, 0);

			fr1[j] = ff[i][j] = Math.max(fr0[j] + gapExt, mx_new);
			er1[j - 1] = er;
			er = ee[i][j] = Math.max(er + gapExt, my_new);

			if (h == (hr0[j - 1] + s) && is_match(dseq[i], qseq[j]))
				trace_mat[i][j] |= (1 << 1);// #define LAL_MASK_MATCH         (1<<1)
			if (h == (hr0[j - 1] + s) && !is_match(dseq[i], qseq[j]))
				trace_mat[i][j] |= (1 << 0);// #define LAL_MASK_MISMATCH      (1<<0)
			if (h == fr0[j - 1] + s)
				trace_mat[i][j] |= (1 << 2); // #define LAL_MASK_GAP_OPEN_LEFT (1<<3)
			if (h == er0[j - 1] + s)
				trace_mat[i][j] |= (1 << 3);// #define LAL_MASK_GAP_OPEN_UP   (1<<2)
			if (h == 0)
				trace_mat[i][j] |= (1 << 6); // #define LAL_MASK_ZERO          (1<<6)
			/*
            var arr = [m_new, ff[i - 1][j - 1] + s, ee[i - 1][j - 1] + s, 0];

            var trace = arr.indexOf(Math.max.apply(Math, arr));
            trace_mat[i][j] = trace;
            */
		}
		er0 = er1;
		for (var j = 0; j < l2; ++j) {
			fr0[j] = fr1[j];
		}
		for (var j = 0; j < l2; ++j) {
			hr0[j] = hr1[j];
		}
		//		fr0 = fr1;
	}
	/* console.log(score_mat); */
	var mscore = h[0][0];
	for (i = 0; i < l1; ++i) {
		for (j = 0; j < l2; ++j) {
			if (mscore < h[i][j])
				mscore = h[i][j];
		}
	}
	/* console.log(mscore); */

	return [mscore, h, trace_mat];
}

function CalculateSWandDrawComp(seq_1, seq_2, matrix, gapOpen, gapExt) {
	var sequence_1 = seq_1.split(" ");
	var sequence_2 = seq_2.split(" ");
	var subtitution = { method: score, match: 1.0, mismatch: -1.0 }; /* define Match/Mismatch Scoring model */
	if (matrix == 'identityNuc')
		subtitution = { method: S, submatrix: ScoreT[matrix] }; /* define Substitution Scoring model */
	/* Append dash in beginning for first row / column */
	sequence_1 = ["-"].concat(sequence_1);
	sequence_2 = ["-"].concat(sequence_2);

	//    var search_profile = { S: subtitution, gap: gapOpen }; /*define search profile*/
	//    ret = sw_linear_gap(search_profile, sequence_1, sequence_2);
	//console.log('max score: ' + ret[0]);

	//    var search_profile = { S: subtitution, gapOpen: gapOpen, gapExt: gapExt };  /*define search profile*/
	//    ret = sw_affine_gotoh_gap(search_profile, sequence_1, sequence_2);
	//    console.log('max score of affine: ' + ret[0]);

	// var search_profile = { S: subtitution, gapOpen: gapOpen, gapExt: gapExt };  /*define search profile*/
	//	ret = sw_affine_gap_gotoh_v2(search_profile, sequence_1, sequence_2);
	//	console.log('max score of affine: ' + ret[0]);

	var search_profile = { S: subtitution, gapOpen: gapOpen, gapExt: gapExt };  /*define search profile*/
	ret = sw_affine_gap_genc(search_profile, sequence_1, sequence_2);
	console.log('max score of affine: ' + ret[0]);

	var search_profile = { S: subtitution, gapOpen: gapOpen, gapExt: gapExt };  /*define search profile*/
	ret2 = sw_affine_gap_v1(search_profile, sequence_1, sequence_2);
	console.log('max score of v1: ' + ret2[0]);
	ret = ret2;
	var search_profile = { S: subtitution, gapOpen: gapOpen, gapExt: gapExt };  /*define search profile*/
	ret3 = sw_affine_gap_sg(search_profile, sequence_1, sequence_2);
	console.log('max score of opt: ' + ret3[0]);
	ret2 = ret3;
	var search_profile = { S: subtitution, gapOpen: -1.0, gapExt: 0.0 };  /*define search profile*/
	ret4 = sw_affine_gap_sg(search_profile, sequence_1, sequence_2);
	console.log('max score of opt -1.0. 0.0: ' + ret4[0]);


	//if (ret2[0] == ret[0])
	//	return;

	if (1) { // compare score matrix
		var lx = sequence_1.length;
		var ly = sequence_2.length;

		for (y = 0; y < ly; ++y) {
			for (x = 0; x < lx; ++x) {
				if (ret[1][x][y] < 0.0)
					ret[1][x][y] = 0.0;
				else
					ret[1][x][y] = ret[1][x][y].toFixed(5)
				if (ret2[1][x][y] < 0.0)
					ret2[1][x][y] = 0.0;
				else
					ret2[1][x][y] = ret2[1][x][y].toFixed(5)
			}
		}

		var max_score_check = 0.0;
		for (y = 0; y < ly; ++y) {
			for (x = 0; x < lx; ++x) {
				ret[1][x][y] = ret[1][x][y] - ret2[1][x][y];
				if (Math.abs(ret[1][x][y]) < 0.00000001)
					ret[1][x][y] = 0.0;
				if (max_score_check < Math.abs(ret[1][x][y]))
					max_score_check = Math.abs(ret[1][x][y]);
			}
		}
		console.log("check: " + max_score_check);

		if (max_score_check <= 0)
			return;
	}


	return;

	//    var search_profile = { S: subtitution, gapOpen: gapOpen, gapExt: gapExt };  /*define search profile*/
	//    ret = sw_affine_gap_v2(search_profile, sequence_1, sequence_2);
	//    console.log('max score of v2: ' + ret[0]);



	var ready_to_sketch = {
		seq1: sequence_1,
		seq2: sequence_2,
		score: ret[0],
		scorematrix: ret[1],
		directions: ret[2]
	};
	var ojal = { "score": 34, "scorematrix": [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 2, 2, 1, 1, 1, 1, 1, 2, 1, 2, 2, 1, 2, 2, 2, 2, 1], [0, 2, 1, 1, 3, 2, 2, 2, 2, 3, 2, 3, 3, 2, 2, 2, 3, 2, 3, 2, 2, 2, 2, 2, 4, 3, 3, 3, 3, 3, 3, 3, 3, 4, 3, 3, 3, 3, 3], [0, 2, 1, 1, 3, 2, 2, 2, 2, 4, 3, 4, 5, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 4, 4, 4, 4, 4], [0, 1, 4, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 3, 3, 3, 3, 3, 3, 6, 5, 6, 6, 5, 5, 5, 5, 5, 6, 5, 6, 6, 5, 7, 6, 6, 6, 6], [0, 1, 3, 3, 2, 5, 4, 4, 5, 4, 4, 4, 4, 4, 6, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 7, 7, 6, 6, 8, 7, 7, 7, 7, 7, 7, 7, 8], [0, 1, 3, 2, 2, 4, 4, 3, 6, 5, 5, 5, 5, 5, 6, 5, 5, 7, 6, 6, 6, 6, 6, 6, 6, 6, 7, 9, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 9], [0, 1, 3, 2, 2, 4, 3, 3, 5, 5, 4, 4, 4, 4, 5, 5, 4, 6, 6, 5, 8, 7, 8, 8, 7, 7, 7, 8, 8, 10, 9, 10, 10, 9, 10, 10, 10, 10, 9], [0, 1, 3, 2, 2, 4, 3, 3, 5, 4, 4, 4, 4, 4, 5, 4, 4, 6, 5, 5, 7, 7, 9, 10, 9, 9, 9, 9, 9, 10, 9, 11, 12, 11, 11, 12, 12, 12, 11], [0, 1, 3, 2, 2, 4, 3, 3, 5, 4, 4, 4, 4, 4, 6, 5, 5, 6, 5, 5, 7, 6, 8, 9, 9, 8, 11, 11, 10, 10, 12, 11, 11, 11, 11, 11, 11, 11, 14], [0, 1, 3, 5, 4, 4, 6, 5, 5, 5, 6, 5, 5, 6, 5, 8, 7, 7, 7, 7, 7, 9, 8, 9, 8, 11, 10, 10, 13, 12, 12, 12, 12, 12, 12, 12, 12, 12, 13], [0, 1, 3, 5, 4, 4, 6, 8, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 9, 8, 9, 8, 9, 8, 10, 10, 10, 12, 12, 11, 11, 11, 11, 11, 11, 11, 11, 13], [0, 1, 3, 4, 4, 4, 5, 7, 7, 6, 6, 6, 6, 6, 6, 7, 6, 6, 6, 8, 11, 10, 11, 10, 10, 10, 10, 10, 12, 14, 13, 13, 13, 13, 13, 13, 13, 13, 13], [0, 1, 3, 5, 4, 4, 6, 7, 6, 6, 8, 7, 7, 8, 7, 8, 7, 7, 7, 8, 10, 13, 12, 12, 12, 12, 12, 12, 12, 13, 13, 12, 12, 12, 12, 12, 12, 12, 13], [0, 1, 3, 4, 4, 4, 5, 7, 6, 6, 7, 7, 6, 7, 7, 7, 7, 6, 6, 8, 10, 12, 15, 14, 14, 14, 14, 14, 14, 14, 14, 15, 14, 14, 14, 14, 14, 14, 14], [0, 1, 3, 4, 3, 6, 5, 7, 9, 8, 8, 8, 8, 8, 9, 8, 8, 9, 8, 8, 10, 12, 14, 14, 13, 13, 16, 16, 15, 15, 16, 15, 15, 15, 15, 15, 15, 15, 16], [0, 2, 3, 4, 6, 5, 5, 7, 8, 11, 10, 10, 10, 10, 10, 10, 10, 10, 11, 10, 10, 12, 14, 13, 16, 15, 15, 15, 15, 15, 15, 15, 15, 17, 16, 16, 16, 16, 16], [0, 1, 3, 4, 5, 8, 7, 7, 9, 10, 10, 9, 9, 9, 12, 11, 11, 12, 11, 11, 11, 12, 14, 13, 15, 15, 17, 17, 16, 16, 17, 16, 16, 16, 16, 16, 16, 16, 18], [0, 1, 3, 4, 5, 7, 7, 7, 8, 10, 9, 9, 9, 9, 11, 11, 10, 11, 11, 10, 13, 12, 14, 16, 15, 15, 16, 16, 16, 18, 17, 19, 18, 18, 18, 18, 18, 18, 18], [0, 2, 3, 4, 6, 7, 6, 7, 8, 10, 9, 11, 11, 10, 11, 10, 13, 12, 13, 12, 12, 12, 14, 15, 18, 17, 17, 17, 17, 17, 17, 18, 18, 20, 19, 19, 19, 19, 19], [0, 1, 4, 4, 5, 7, 6, 7, 8, 10, 9, 10, 10, 10, 11, 10, 12, 12, 12, 12, 14, 13, 14, 16, 17, 17, 16, 16, 16, 19, 18, 19, 20, 19, 22, 21, 21, 21, 21], [0, 1, 3, 4, 5, 7, 6, 7, 9, 10, 9, 10, 10, 9, 12, 11, 12, 14, 13, 13, 13, 13, 14, 15, 17, 16, 19, 18, 18, 18, 21, 20, 20, 20, 21, 21, 20, 20, 23], [0, 1, 3, 4, 5, 7, 6, 7, 9, 10, 9, 10, 10, 9, 11, 11, 12, 14, 13, 13, 13, 13, 14, 15, 17, 16, 18, 21, 20, 20, 20, 20, 20, 20, 21, 20, 20, 20, 22], [0, 1, 3, 4, 5, 7, 6, 7, 8, 10, 9, 10, 10, 9, 11, 10, 12, 13, 13, 12, 15, 14, 15, 16, 17, 16, 18, 20, 20, 22, 21, 22, 22, 21, 22, 23, 22, 22, 22], [0, 1, 3, 4, 5, 7, 6, 7, 9, 10, 9, 10, 10, 9, 11, 10, 12, 14, 13, 13, 14, 14, 14, 15, 17, 16, 18, 20, 19, 21, 24, 23, 23, 23, 23, 23, 23, 23, 24], [0, 2, 3, 4, 6, 7, 6, 7, 8, 11, 10, 11, 12, 11, 11, 11, 12, 13, 16, 15, 15, 15, 15, 15, 17, 16, 18, 20, 19, 21, 23, 23, 22, 25, 24, 24, 24, 24, 24], [0, 1, 4, 4, 5, 7, 6, 7, 8, 10, 10, 10, 11, 11, 11, 10, 12, 13, 15, 15, 17, 16, 17, 17, 17, 16, 18, 20, 19, 21, 23, 25, 25, 24, 27, 26, 26, 26, 26], [0, 1, 3, 4, 5, 7, 6, 7, 8, 10, 9, 10, 11, 10, 11, 10, 12, 13, 15, 14, 17, 16, 18, 19, 18, 18, 18, 20, 19, 21, 23, 25, 27, 26, 26, 29, 28, 28, 28], [0, 2, 3, 4, 6, 7, 6, 7, 8, 10, 9, 11, 12, 11, 11, 11, 12, 13, 15, 14, 16, 16, 17, 18, 21, 20, 20, 20, 20, 21, 23, 24, 26, 29, 28, 28, 28, 28, 28], [0, 1, 3, 4, 5, 8, 7, 7, 9, 10, 9, 10, 11, 11, 13, 12, 12, 14, 15, 14, 16, 15, 17, 18, 20, 20, 22, 22, 21, 21, 23, 24, 26, 28, 28, 28, 27, 27, 30], [0, 1, 3, 4, 5, 7, 7, 7, 8, 10, 9, 10, 11, 10, 12, 12, 12, 13, 15, 14, 16, 15, 17, 19, 20, 19, 21, 21, 21, 23, 23, 25, 26, 28, 30, 30, 30, 29, 29], [0, 1, 3, 5, 5, 7, 9, 9, 8, 10, 12, 11, 11, 13, 12, 14, 13, 13, 15, 17, 16, 18, 17, 18, 20, 22, 21, 21, 23, 22, 23, 24, 26, 28, 29, 29, 29, 29, 29], [0, 1, 3, 4, 5, 7, 8, 8, 11, 10, 11, 11, 11, 12, 15, 14, 14, 15, 15, 16, 16, 17, 17, 18, 20, 21, 24, 23, 23, 23, 24, 24, 26, 28, 29, 29, 29, 28, 31], [0, 1, 3, 4, 5, 7, 8, 8, 10, 10, 11, 10, 11, 12, 14, 14, 13, 14, 15, 16, 18, 17, 19, 19, 20, 21, 23, 23, 22, 25, 24, 26, 26, 28, 30, 31, 31, 31, 30], [0, 1, 3, 4, 5, 7, 8, 8, 10, 10, 11, 10, 11, 12, 14, 13, 13, 15, 15, 16, 17, 17, 18, 18, 20, 21, 23, 25, 24, 24, 27, 26, 26, 28, 29, 30, 30, 30, 33], [0, 1, 3, 5, 5, 7, 9, 10, 10, 10, 12, 11, 11, 13, 14, 16, 15, 15, 15, 17, 17, 19, 18, 18, 20, 22, 23, 24, 27, 26, 26, 26, 26, 28, 29, 30, 30, 30, 32], [0, 1, 3, 4, 5, 7, 8, 9, 10, 10, 11, 11, 11, 12, 14, 15, 15, 14, 15, 16, 19, 18, 21, 20, 20, 21, 23, 24, 26, 29, 28, 28, 28, 28, 30, 31, 32, 32, 32], [0, 2, 3, 4, 6, 7, 8, 9, 10, 12, 11, 13, 13, 12, 14, 15, 17, 16, 16, 16, 18, 18, 20, 20, 22, 21, 23, 24, 26, 28, 28, 27, 27, 30, 29, 30, 31, 31, 32], [0, 1, 3, 5, 5, 7, 9, 10, 10, 11, 14, 13, 13, 15, 14, 16, 16, 16, 15, 18, 18, 20, 20, 19, 21, 24, 23, 24, 26, 28, 27, 27, 27, 29, 29, 30, 31, 31, 32], [0, 1, 3, 4, 5, 7, 8, 9, 10, 11, 13, 13, 12, 14, 14, 15, 16, 15, 15, 17, 20, 19, 22, 22, 21, 23, 23, 24, 26, 28, 27, 29, 29, 29, 31, 31, 32, 33, 32], [0, 1, 3, 5, 5, 7, 9, 10, 10, 11, 13, 12, 12, 14, 14, 16, 16, 15, 15, 17, 19, 22, 21, 21, 21, 23, 23, 24, 26, 28, 27, 28, 28, 29, 30, 30, 31, 32, 32], [0, 1, 3, 4, 5, 7, 8, 9, 12, 11, 13, 12, 12, 14, 16, 15, 16, 18, 17, 17, 19, 21, 21, 21, 21, 23, 25, 25, 26, 28, 30, 29, 29, 29, 30, 30, 31, 32, 34]], "directions": [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 77, 2, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 2, 9, 2, 2, 9, 9, 9, 9, 9, 2, 9, 2, 2, 9, 2, 2, 2, 2, 9], [0, 2, 13, 9, 2, 9, 9, 9, 9, 2, 9, 2, 2, 9, 9, 9, 2, 9, 2, 9, 9, 9, 9, 9, 2, 9, 9, 9, 9, 9, 9, 9, 9, 2, 9, 9, 9, 9, 9], [0, 2, 13, 9, 2, 9, 9, 9, 9, 2, 9, 2, 2, 9, 9, 9, 10, 9, 10, 9, 9, 9, 9, 9, 10, 9, 9, 9, 9, 9, 9, 9, 9, 2, 9, 9, 9, 9, 9], [0, 5, 2, 9, 9, 9, 9, 9, 9, 13, 9, 13, 5, 1, 13, 13, 13, 13, 13, 13, 2, 9, 2, 2, 9, 9, 9, 9, 9, 2, 9, 2, 2, 9, 2, 10, 10, 10, 9], [0, 5, 5, 1, 13, 2, 9, 9, 2, 9, 9, 9, 13, 9, 2, 9, 9, 10, 9, 9, 13, 9, 13, 13, 9, 9, 2, 2, 9, 9, 2, 9, 9, 9, 9, 9, 9, 9, 2], [0, 5, 5, 13, 13, 6, 1, 13, 2, 9, 9, 9, 9, 9, 2, 9, 9, 2, 9, 9, 9, 9, 9, 9, 9, 9, 2, 2, 9, 9, 10, 9, 9, 9, 9, 9, 9, 9, 2], [0, 5, 6, 13, 13, 5, 13, 13, 5, 1, 13, 13, 13, 13, 5, 1, 13, 5, 1, 13, 2, 9, 2, 2, 9, 9, 9, 5, 1, 2, 9, 2, 2, 9, 2, 2, 2, 2, 9], [0, 5, 6, 13, 13, 5, 13, 13, 5, 13, 13, 13, 13, 13, 5, 13, 13, 5, 13, 13, 6, 1, 2, 2, 9, 9, 9, 9, 9, 2, 9, 2, 2, 9, 10, 2, 2, 2, 9], [0, 5, 5, 13, 13, 6, 13, 13, 6, 13, 13, 13, 13, 13, 2, 9, 9, 6, 13, 13, 5, 13, 5, 5, 1, 13, 2, 2, 9, 9, 2, 9, 13, 9, 9, 13, 13, 13, 2], [0, 5, 5, 2, 9, 13, 2, 10, 13, 9, 2, 9, 9, 2, 13, 2, 9, 9, 9, 10, 13, 2, 13, 5, 13, 2, 13, 13, 2, 9, 9, 9, 9, 9, 9, 9, 9, 9, 5], [0, 5, 5, 2, 9, 13, 2, 2, 9, 9, 10, 9, 9, 10, 9, 14, 9, 9, 9, 2, 9, 2, 13, 5, 13, 6, 5, 5, 6, 1, 13, 13, 13, 13, 13, 13, 13, 13, 5], [0, 5, 6, 5, 1, 5, 5, 5, 1, 13, 13, 13, 13, 13, 13, 5, 13, 13, 13, 5, 2, 9, 2, 10, 9, 13, 13, 13, 5, 2, 9, 10, 10, 9, 10, 10, 10, 10, 13], [0, 5, 5, 2, 9, 13, 2, 6, 13, 13, 2, 9, 9, 2, 9, 2, 9, 9, 9, 6, 5, 2, 9, 9, 9, 10, 9, 9, 14, 5, 1, 13, 13, 13, 13, 13, 13, 13, 5], [0, 5, 6, 5, 1, 5, 5, 5, 13, 13, 5, 1, 13, 5, 1, 5, 1, 13, 13, 5, 6, 5, 2, 10, 9, 9, 9, 9, 9, 10, 9, 2, 10, 9, 10, 10, 10, 10, 9], [0, 5, 5, 5, 13, 2, 13, 5, 2, 9, 9, 9, 9, 9, 2, 9, 9, 2, 9, 13, 5, 5, 5, 1, 13, 13, 2, 2, 9, 9, 2, 9, 9, 9, 9, 9, 9, 9, 2], [0, 2, 5, 5, 2, 13, 13, 5, 5, 2, 9, 10, 10, 9, 9, 9, 10, 9, 2, 9, 13, 5, 5, 13, 2, 9, 13, 13, 9, 9, 13, 9, 9, 2, 9, 9, 9, 9, 9], [0, 5, 5, 5, 5, 2, 9, 13, 2, 5, 1, 13, 13, 13, 2, 9, 9, 2, 9, 9, 9, 5, 5, 13, 5, 1, 2, 2, 9, 9, 2, 9, 9, 13, 9, 9, 9, 9, 2], [0, 5, 6, 5, 5, 5, 1, 5, 5, 5, 13, 13, 13, 13, 5, 1, 13, 5, 1, 13, 2, 13, 6, 2, 13, 9, 5, 5, 1, 2, 9, 2, 10, 9, 10, 10, 10, 10, 9], [0, 2, 5, 5, 2, 5, 13, 5, 5, 6, 13, 2, 2, 9, 5, 13, 2, 9, 2, 9, 13, 13, 5, 5, 2, 9, 9, 9, 9, 13, 9, 5, 1, 2, 9, 9, 9, 9, 9], [0, 5, 2, 5, 5, 5, 13, 5, 5, 5, 13, 5, 5, 1, 5, 13, 5, 1, 5, 1, 2, 9, 6, 2, 5, 1, 13, 13, 13, 2, 9, 2, 2, 13, 2, 10, 10, 10, 9], [0, 5, 5, 5, 5, 6, 13, 5, 2, 5, 13, 5, 5, 13, 2, 9, 5, 2, 9, 9, 13, 9, 5, 5, 5, 13, 2, 10, 9, 13, 2, 9, 9, 9, 5, 1, 13, 13, 2], [0, 5, 5, 5, 5, 6, 13, 5, 2, 5, 13, 5, 5, 13, 6, 1, 5, 2, 9, 9, 13, 9, 5, 5, 5, 13, 6, 2, 9, 9, 14, 9, 9, 9, 5, 13, 13, 13, 6], [0, 5, 6, 5, 5, 5, 13, 5, 5, 5, 13, 5, 5, 13, 5, 13, 5, 5, 1, 13, 2, 9, 2, 2, 5, 13, 5, 5, 1, 2, 9, 2, 2, 9, 2, 2, 10, 10, 13], [0, 5, 5, 5, 5, 6, 13, 5, 2, 5, 13, 5, 5, 13, 6, 13, 5, 2, 9, 9, 5, 1, 5, 5, 5, 13, 6, 6, 13, 5, 2, 9, 9, 9, 9, 9, 9, 9, 2], [0, 2, 5, 5, 2, 5, 13, 5, 5, 2, 9, 2, 2, 9, 13, 9, 6, 5, 2, 9, 9, 9, 9, 13, 6, 13, 5, 5, 13, 5, 5, 1, 13, 2, 9, 9, 9, 9, 9], [0, 5, 2, 5, 5, 5, 13, 5, 5, 5, 1, 5, 5, 1, 5, 13, 5, 5, 5, 1, 2, 9, 2, 2, 5, 13, 5, 5, 13, 6, 5, 2, 2, 13, 2, 10, 10, 10, 9], [0, 5, 6, 5, 5, 5, 13, 5, 5, 5, 13, 5, 5, 13, 5, 13, 5, 5, 5, 13, 2, 9, 2, 2, 9, 9, 13, 5, 13, 6, 5, 2, 2, 9, 14, 2, 10, 10, 9], [0, 2, 5, 5, 2, 5, 13, 5, 5, 6, 13, 2, 2, 9, 13, 9, 6, 5, 6, 13, 5, 1, 5, 5, 2, 9, 9, 13, 9, 5, 5, 5, 5, 2, 9, 13, 9, 9, 9], [0, 5, 5, 5, 5, 2, 9, 13, 2, 5, 13, 5, 5, 1, 2, 9, 13, 2, 5, 13, 5, 13, 5, 5, 5, 1, 2, 2, 9, 13, 6, 5, 5, 5, 1, 5, 13, 13, 2], [0, 5, 6, 5, 5, 5, 1, 5, 5, 5, 13, 5, 5, 13, 5, 1, 5, 5, 5, 13, 6, 13, 6, 2, 5, 13, 5, 5, 1, 2, 5, 2, 6, 5, 2, 2, 2, 10, 13], [0, 5, 5, 2, 5, 5, 2, 2, 13, 5, 2, 9, 13, 2, 13, 2, 9, 13, 5, 2, 13, 2, 13, 5, 5, 2, 13, 13, 2, 13, 5, 5, 5, 5, 5, 5, 5, 1, 5], [0, 5, 5, 5, 5, 6, 5, 5, 2, 13, 5, 1, 5, 5, 2, 9, 9, 2, 5, 5, 5, 5, 5, 5, 5, 5, 2, 10, 9, 9, 2, 5, 5, 5, 5, 5, 5, 13, 2], [0, 5, 6, 5, 5, 5, 5, 5, 5, 5, 5, 13, 5, 5, 5, 1, 13, 5, 5, 5, 2, 13, 2, 2, 5, 5, 5, 1, 13, 2, 9, 2, 6, 5, 2, 2, 2, 2, 13], [0, 5, 5, 5, 5, 6, 5, 5, 6, 5, 5, 13, 5, 5, 6, 13, 13, 2, 5, 5, 5, 5, 5, 5, 5, 5, 6, 2, 9, 13, 2, 9, 13, 5, 5, 5, 5, 5, 2], [0, 5, 5, 2, 5, 5, 2, 2, 5, 5, 2, 9, 13, 2, 5, 2, 9, 9, 13, 2, 5, 2, 13, 13, 5, 2, 5, 5, 2, 9, 13, 9, 13, 5, 5, 5, 5, 5, 5], [0, 5, 6, 5, 5, 5, 5, 5, 5, 5, 5, 1, 5, 5, 5, 5, 1, 13, 5, 5, 2, 13, 2, 10, 13, 5, 5, 5, 5, 2, 9, 10, 10, 13, 2, 2, 2, 2, 5], [0, 2, 5, 5, 2, 5, 5, 5, 5, 2, 13, 2, 2, 13, 5, 5, 2, 9, 10, 13, 5, 5, 5, 1, 2, 13, 5, 5, 5, 5, 1, 13, 13, 2, 13, 5, 5, 5, 5], [0, 5, 5, 2, 5, 5, 2, 2, 5, 5, 2, 9, 9, 2, 13, 2, 5, 1, 13, 2, 5, 2, 5, 13, 5, 2, 13, 5, 6, 5, 13, 13, 13, 5, 5, 5, 5, 5, 5], [0, 5, 6, 5, 5, 5, 5, 5, 5, 5, 5, 1, 13, 5, 5, 5, 5, 13, 13, 5, 2, 13, 2, 2, 13, 5, 5, 5, 5, 6, 13, 2, 2, 5, 2, 2, 2, 2, 13], [0, 5, 5, 2, 5, 5, 2, 2, 5, 5, 6, 13, 13, 6, 5, 2, 5, 13, 13, 6, 5, 2, 13, 13, 13, 6, 5, 5, 6, 5, 13, 5, 5, 5, 5, 5, 5, 5, 5], [0, 5, 5, 5, 5, 6, 5, 5, 2, 13, 5, 13, 13, 5, 2, 13, 5, 2, 9, 13, 5, 5, 5, 5, 5, 5, 2, 2, 5, 5, 2, 9, 9, 13, 5, 5, 5, 5, 2]], "seq1": "CAACTTCCTGGCGCTATCACTTCTACCATCGTCTGCAGCGT", "seq2": "acgatggtagaagtgatagcgccagttgctccacccct" }

	var seq_10 = [];
	for (var i = 0 ; i < ojal.seq1.length; i++) {
		seq_10.push(ojal.seq1.charAt(i));
	}
	var seq_20 = [];
	for (var i = 0, seq_20; i < ojal.seq2.length; i++) {
		seq_20.push(ojal.seq2.charAt(i));
	}
	seq_10 = ["-"].concat(seq_10);
	seq_20 = ["-"].concat(seq_20);

	ojal.seq1 = seq_10;
	ojal.seq2 = seq_20;

	// ready_to_sketch = ojal;

	sketch_matrixes(ready_to_sketch);
}