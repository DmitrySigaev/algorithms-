// Create empty matrix
var Matrix = function (rows, cols) {
    var arr = [], row = [];
    while (cols--) row.push(0);
    while (rows--) arr.push(row.slice());
    return arr;
};

var score = function (x, y, match, mism) {
    if (x == y) {
        return match;
    } else {
        return mism;
    }
};

/*
 * Sequence alignments are used in different area of computer science.
 * Main feature of alignment is a gap function because complexity of
 * computation depends on this function. The most common and simple
 * case is using a constant gap penalties.
 */
var sw_constant_gap = function (s1, s2, match, mismatch, gap) {
    var l1 = s1.length;
    var l2 = s2.length;
    var score_mat = Matrix(l1, l2);
    var trace_mat = Matrix(l1, l2);

    for (i = 0; i < l1; ++i) {

        for (j = 0; j < l2; ++j) {

            // This is the first row / column which is all zeros
            if (i == 0 || j == 0) {
                score_mat[i][j] = 0;
                trace_mat[i][j] = 3;
                continue;

            } else {
                var d_last = score_mat[i - 1][j - 1];
                var u_last = score_mat[i - 1][j];
                var l_last = score_mat[i][j - 1];
            }
            var d_new = d_last + score(s1[i], s2[j], match, mismatch);
            var u_new = u_last + gap;
            var l_new = l_last + gap;
            score_mat[i][j] = Math.max(d_new, u_new, l_new, 0);
            var arr = [d_new, u_new, l_new, 0];
            var trace = arr.indexOf(Math.max.apply(Math, arr));
            trace_mat[i][j] = trace;
        }
    }
    console.log(score_mat);
    var mscore = 0;
    for (i = 0; i < l1; ++i) {
        for (j = 0; j < l2; ++j) {
            if (mscore < score_mat[i][j])
                mscore = score_mat[i][j];
        }
    }
    console.log(mscore);

    return [mscore, score_mat, trace_mat];
}

function CalculateSWandDraw(seq_1, seq_2, match, misMatch, gapOpen, gapExt) {
    var sequence_1 = seq_1.split(" ");
    var sequence_2 = seq_2.split(" ");

    // Append dash in beginning for first row / column
    sequence_1 = ["-"].concat(sequence_1);
    sequence_2 = ["-"].concat(sequence_2);

    //Number of rows and columns
    var nrow = sequence_1.length;
    var ncol = sequence_2.length;

    ret = sw_constant_gap(sequence_1, sequence_1, 1.0, -1.0, -1.0);
    console.log(ret[0] == 4);
}