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
            var d_new = d_last + substitution.method(substitution, s1[i], s2[j]);
            var u_new = u_last + gap;
            var l_new = l_last + gap;
            score_mat[i][j] = Math.max(d_new, u_new, l_new, 0);
            var arr = [d_new, u_new, l_new, 0];
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
    var ff = Matrix(l1, l2);

    for (i = 0; i < l1; ++i) {

        for (j = 0; j < l2; ++j) {

            // This is the first row / column which is all zeros
            if (i == 0 || j == 0) {
                if (i == 0)
                    ee[i][j] = 0;
                if (j == 0)
                    ff[i][j] = 0;
                score_mat[i][j] = 0;
                trace_mat[i][j] = 3;
                continue;

            } else {
                var e_last = ee[i][j - 1];
                var f_last = ff[i - 1][j];
                var d_last = score_mat[i - 1][j - 1];
                var u_last = score_mat[i - 1][j];
                var l_last = score_mat[i][j - 1];
            }
            var d_new = d_last + substitution.method(substitution, dseq[i], qseq[j]);
            var u_new = u_last + gapOpen;
            var l_new = l_last + gapOpen;
            var e_new = e_last + gapExt;
            var f_new = f_last + gapExt;
            ee[i][j] = Math.max(e_new, l_new);
            ff[i][j] = Math.max(f_new, u_new);

            score_mat[i][j] = Math.max(d_new, ee[i][j], ff[i][j], 0);
            var arr = [d_new, ff[i][j], ee[i][j], 0];
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

function CalculateSWandDraw(seq_1, seq_2, matrix, gapOpen, gapExt) {
    var sequence_1 = seq_1.split(" ");
    var sequence_2 = seq_2.split(" ");
    var subtitution = { method: score, match: 1.0, mismatch: -1.0 }; /* define Match/Mismatch Scoring model */
    if (matrix == 'identityNuc')
        subtitution = { method: S, submatrix: ScoreT[matrix] }; /* define Substitution Scoring model */
    /* Append dash in beginning for first row / column */
    sequence_1 = ["-"].concat(sequence_1);
    sequence_2 = ["-"].concat(sequence_2);

    var search_profile = { S: subtitution, gap: gapOpen }; /*define search profile*/
    ret = sw_linear_gap(search_profile, sequence_1, sequence_2);
    console.log('max score: ' + ret[0]);

    var search_profile = { S: subtitution, gapOpen: gapOpen, gapExt: gapExt };
    ret = sw_affine_gotoh_gap(search_profile, sequence_1, sequence_2);
    console.log('max score of affine: ' + ret[0]);
}