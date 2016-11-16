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
var score = function (search_profile, x, y) {
    if (x.toUpperCase() == y.toUpperCase()) {
        return search_profile.match || 1;
    } else {
        return search_profile.mismatch || -1;
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
    var substitution_function = search_profile.S || score;
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
            var d_new = d_last + substitution_function(search_profile, s1[i], s2[j]);
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

function CalculateSWandDraw(seq_1, seq_2, matrix, gapOpen, gapExt) {
    var sequence_1 = seq_1.split(" ");
    var sequence_2 = seq_2.split(" ");
    var subtitution = { subfun: score, match: 1.0, mismatch: -1.0 }
    /* Append dash in beginning for first row / column */
    sequence_1 = ["-"].concat(sequence_1);
    sequence_2 = ["-"].concat(sequence_2);

    var search_profile = { S: score, match: 1.0, mismatch: -1.0, gap: -1.0 };
    ret = sw_linear_gap(search_profile, sequence_1, sequence_2);
    console.log('max score: ' + ret[0]);
}