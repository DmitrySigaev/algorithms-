/* sketch_matrix.js
 * Copyright 2016 Dmitry Sigaev
 *
 * Released under the MIT license -- see MIT-LICENSE for details
 */
/* visualization of scoring matrix */

var upp = (1 << 2);// #define LAL_MASK_GAP_OPEN_UP   (1<<2)
var lef = (1 << 3); // #define LAL_MASK_GAP_OPEN_LEFT (1<<3)
var dig = (1 << 1);// #define LAL_MASK_MATCH         (1<<1)
var misdig = (1 << 0);// #define LAL_MASK_MISMATCH         (1<<0)
/*
#define LAL_MASK_MISMATCH      (1<<0)
#define LAL_MASK_MATCH         (1<<1)
#define LAL_MASK_GAP_OPEN_UP   (1<<2)
#define LAL_MASK_GAP_OPEN_LEFT (1<<3)
#define LAL_MASK_GAP_EXT_UP    (1<<4)
#define LAL_MASK_GAP_EXT_LEFT  (1<<5)
#define LAL_MASK_ZERO          (1<<6)
*/

var getMaxSymbolsSize = function (sequence, font, fontSize) {
	/* Measure symbol size in sequence. */
	var testSvg = d3.select("body").append("svg").attr("id", "getMaxSymbolsSize");
	/*
	 The <g> SVG element is a container used to group other SVG elements. 
	 Transformations applied to the <g> element are performed on all of its child elements, 
	 and any of its attributes are inherited by its child elements.
	 It can also group multiple elements to be referenced later with the <use> element.
	*/
	var text = testSvg.selectAll("g")
				 .data(sequence)
				 .enter()
				 .append("g")
				 .attr("id", "test")
				 .append("text")
				 .text(
				 function (d) {
				 	return d;
				 }
				 )
				 .attr("font-size", fontSize)
				 .style("font-family", font);

	var maxWidth = 0;
	var maxHeight = 0;
	var counter = 0;
	text.each(function () {
		var width = this.getBBox().width;
		var height = this.getBBox().height;
		if (width > maxWidth) {
			maxWidth = width;
		}
		if (height > maxHeight) {
			maxHeight = height;
		}
		counter++;
	});

	if (counter !== sequence.length)
		console.error("Please check the input sequence: " + sequence);
	d3.select("#getMaxSymbolsSize").remove();
	return { w: maxWidth, h: maxHeight, font: font, font_size: fontSize };

}

var setSvgArea = function (yseq, xseq, symSize) {
	/* Number of rows and columns */
	var nrows = yseq.length;
	var ncols = xseq.length;
	var containerWidth = parseFloat(d3.select('#main_frame').style("width"));

	var visualWidth = 0.8 * containerWidth;
	var size = 0.8 * visualWidth;

	if (nrows == ncols) {
		var w = size;
		var h = size;
	} else if (nrows > ncols) {
		var h = size;
		var w = size * ncols / nrows;
	} else if (nrows < ncols) {
		var w = size;
		var h = size * nrows / ncols;
	}

	var left_margin = 1.1 * Math.max(symSize.w, symSize.h); /*can contains size of string box instead of symbol box */
	var top_margin = left_margin;
	/* calculte sixe of cell */
	var cell_width = (w - left_margin) / ncols;
	var cell_height = (h - top_margin) / nrows;

	/* Padding between matix cells */
	var padding = 0.05;
	var row_padding = padding * cell_height
	var col_padding = padding * cell_width;
	var color = "orchid"; /*Iana's favorite color */ /* http://www.w3schools.com/colors/colors_names.asp */
	var color = "orange"; /*George's favorite color */
	return { svg: { nrow: nrows, ncol: ncols, w: w, h: h, margin: { left: left_margin, top: top_margin } }, cell: { w: cell_width, h: cell_height, padding: { row: row_padding, col: col_padding }, color: color }, symbol: symSize };
}

/* 
 Get the coordinates of the upper left corner of a cell from the cell index. 
 Sequentially indexed row by row
 */
var get_cell_coord = function (sa, i, coord) {
	var col = i % sa.svg.ncol;
	var row = Math.floor((i / sa.svg.ncol));
	if (coord == 0 || coord == 2 /*"x" */) {
		var returnValue = col * sa.cell.w + sa.cell.padding.col + sa.svg.margin.left;
	} else if (coord == 1 || coord == 3/* "y" */) {
		var returnValue = row * sa.cell.h + sa.cell.padding.row + sa.svg.margin.top;
	}
	return returnValue;
}

/* Get the endpoint of a trace line depending on the trace */
function get_line_coord(sa, d, i, coord) {
	var halfCellPadW = sa.cell.padding.col + 0.5 * (sa.cell.w - sa.cell.padding.col);
	var halfCellPadH = sa.cell.padding.row + 0.5 * (sa.cell.h - sa.cell.padding.row);

	var additoin = [/*zero      (0 << 0)*/[0.0, 0.0, 0.0, 0.0],
					/*misdig    (1 << 0)*/[0.0, 0.0, 0.6 * halfCellPadW, 0.6 * halfCellPadH],
					/*diag_dir  (1 << 1)*/[0.0, 0.0, 0.6 * halfCellPadW, 0.6 * halfCellPadH],
					/*up_dir    (1 << 2)*/[halfCellPadW, 0.0, halfCellPadW, 0.6 * halfCellPadH],
					/*left_dir  (1 << 3)*/[0.0, halfCellPadH, 0.6 * halfCellPadW, halfCellPadH],
					/*up_dir_e  (1 << 4)*/[halfCellPadW, 0.0, halfCellPadW, 0.6 * halfCellPadH],
					/*left_dir_e(1 << 5)*/[0.0, halfCellPadH, 0.6 * halfCellPadW, halfCellPadH],
					/*zero_dir  (1 << 6)*/[0.0, 0.0, 0.0, 0.0]];

	for (var sh = 0; sh <= 6; sh++) {
		if (!!(d['trace'] & (1 << sh)) && d['status'][sh] == 1) {
			if (coord == 3)
				d['status'][sh] = 0;
			return get_cell_coord(sa, i, coord) + additoin[sh + 1][coord];
		}
		else
			d['status'][sh] = 0;
	}
	return get_cell_coord(sa, i, coord);
}

/* redraw trace line */
var redraw = function (sa, matrix) {
	matrix.selectAll("g")
				 .append("line")
				 .attr("x1", function (d, i) { return get_line_coord(sa, d, i, 0); })
				 .attr("y1", function (d, i) { return get_line_coord(sa, d, i, 1); })
				 .attr("x2", function (d, i) { return get_line_coord(sa, d, i, 2); })
				 .attr("y2", function (d, i) { return get_line_coord(sa, d, i, 3); })
				 .attr("stroke-width", 1)
				 .attr("stroke", function (d, i) {
				 	if (d['trace'] > 0)
				 		console.log('');
				 	//Don't draw a trace line if the score is zero or the choice was zero
				 	if (!!(d['trace'] & 15)) {
				 		if (!!(d['trace'] & 3) && !!(d['trace'] & upp))
				 			return "blue";
				 		if (!!(d['trace'] & upp) && !!(d['trace'] & lef))
				 			return "yellow";
				 		else
				 			return "none";
				 	} else {
				 		return "none";
				 	}
				 });
}

/* Draw trace line */
var DrawTraceLine = function (sa, matrix, data) {
	matrix.selectAll("g")
	.data(data)
	.enter()
	.append("g")
	.append("line")
	.attr("x1", function (d, i) {
		return get_line_coord(sa, d, i, 0);
	})
	.attr("y1", function (d, i) {
		return get_line_coord(sa, d, i, 1);
	})
	.attr("x2", function (d, i) {
		return get_line_coord(sa, d, i, 2);
	})
	.attr("y2", function (d, i) {
		return get_line_coord(sa, d, i, 3);
	})
	.attr("stroke-width", 1)
	.attr("stroke", function (d, i) {
		if (d['trace'] > 0)
			console.log('');
		/* Don't draw a trace line if the score is zero or the choice was zero */
		if (!!(d['trace'] & 15)) {
			if (!!(d['trace'] & 1))
				return "red";
			if (!!(d['trace'] & 2))
				return "green";
			if (!!(d['trace'] & upp))
				return "grey";
			else
				return "black";
		} else {
			return "none";
		}
	});
	redraw(sa, matrix);
	redraw(sa, matrix);
}

/* Draw the scores */
var DrawScores = function (sa, matrix) {
	var halfCellPadW = sa.cell.padding.col + 0.5 * (sa.cell.w - sa.cell.padding.col);
	var halfCellPadH = sa.cell.padding.row + 0.5 * (sa.cell.h - sa.cell.padding.row);

	matrix.selectAll("g")
		.append("text")
		.text(function (d, i) { return d["score"]; })
		.attr("x", function (d, i) {
			var out = get_cell_coord(sa, i, 0) + halfCellPadW - 0.2 * sa.symbol.font_size;
			return out;
		})
		.attr("y", function (d, i) {
			var out = get_cell_coord(sa, i, 1) + halfCellPadH + 0.2 * sa.symbol.font_size;
			return out;
		})
		.style("text-anchor", "middle")
		.style("font-size", 0.3 * sa.cell.w)
		.attr("font-family", "Arial")
		.style("fill", "grey");
}

/* Draw the row names */
var DrawRowNames = function (sa, matrix, yseqence) {
	var names = matrix.selectAll("text.row")
		.data(yseqence)
		.enter()
		.append("text");

	names.text(function (d) { return d; })
		.attr("y", function (d, i) {
			var out = sa.svg.margin.top + 0.5 * sa.cell.h + i * sa.cell.h + sa.cell.padding.row;
			return out;
		})
		.attr("x", 0)
		.attr("id", function (d, i) { return "row_" + i; })
		.style("text-anchor", "right")
		.style("fill", "black")
		.attr("font-family", sa.symbol.font)
		.style("font-size", Math.min(0.3 * sa.cell.h, sa.symbol.font_size));
}

/* Draw the column names */
var DrawColNames = function (sa, matrix, xseqence) {
	var names = matrix.selectAll("text.col")
		.data(xseqence)
		.enter()
		.append("text");

	names.text(function (d) { return d; })
		.attr("x", function (d, i) {
			var out = sa.svg.margin.left + 0.5 * sa.cell.w + i * sa.cell.w + sa.cell.padding.col;
			return out;
		})
		.attr("y", 0)
		.attr("id", function (d, i) { return "col_" + i; })
		.style("text-anchor", "left")
		.style("writing-mode", "vertical-rl")
		.style("fill", "black")
		.attr("font-family", sa.symbol.font)
		.style("font-size", Math.min(0.3 * sa.cell.h, sa.symbol.font_size));
}

/* Draw cells */
var DrawCells = function (sa, matrix, max_score) {
	var tip = d3.select(".g-tip");
	var tipMetric = tip.selectAll(".g-tip-metric")
		.datum(function () { return this.getAttribute("data-name"); });


	matrix.selectAll("g")
		.append("rect")
		.attr("x", function (d, i) { return get_cell_coord(sa, i, 0); })
		.attr("y", function (d, i) { return get_cell_coord(sa, i, 1); })
		.attr("width", (sa.cell.w - sa.cell.padding.col))
		.attr("height", (sa.cell.h - sa.cell.padding.row))
		.style("fill", sa.cell.color)
		.style("stroke", "white")
		.attr("active", false)
		.attr("id", function (d, i) { return "cell_" + i })
		.style("opacity", function (d) {
			var s = d['score'];
			return (((s / max_score) + 0.1) * 0.6);
		})
		.on('mouseover', function (d, i) {
			d3.select(this).style("stroke", "grey")
				.style("stroke-width", 2);

			tip.style("display", null)
				.style("top", get_cell_coord(sa, i, 1) - sa.cell.h + "px")
				.style("left", get_cell_coord(sa, i, 0) + sa.cell.w + "px");

			tipMetric.select(".g-tip-metric-value").text(function (name) {
				switch (name) {
					case "j_1_i_1": return d['score'];
					case "j_i_1": return d['score'] - 1;
					case "j_1_i": return d['score'] - 2;
				}
			});
		})
		.on('mouseout', function (d) {
			d3.select(this).style("stroke", "white")
				.style("stroke-width", 1)
			tip.style("display", "none");
		});
}

function sketch_matrixes(object_alignment) {
	var symSize = getMaxSymbolsSize(object_alignment.seq1, "Arial", 12);
	var sa = setSvgArea(object_alignment.seq1, object_alignment.seq2, symSize);
	console.log(sa);
	// Calculate the matrix. Imported from sw.js

	var score_mat = object_alignment.scorematrix;
	var trace_mat = object_alignment.directions;


	var data = new Array();
	var c = 0;
	var max_score = 0;
	for (i = 0; i < sa.svg.nrow; i++) { /* var nrows = yseq.length = object_alignment.seq1.length */
		for (j = 0; j < sa.svg.ncol; j++) { /* var ncols = xseq.length = object_alignment.seq2 */
			/* get max score for color gradient scaling */
			if (score_mat[i][j] > max_score) {
				max_score = score_mat[i][j];
			}
			/* Store all relevant data in json object array in order to append to DOMS in d3 */
			var dom = {
				"score": score_mat[i][j],
				"trace": trace_mat[i][j],
				"status": [1, 1, 1, 1, 1, 1, 1],
				"active": false,
				"counter": c
			}
			data.push(dom)
			c = c + 1;
		}
	}
	console.log(data);
	if (object_alignment.score !== max_score)
		console.error("Please check the object_alignment!");
	/* Create SVG matrix element */
	var scoreMatrix = d3.select("#main_frame")
						.append("svg")
						.attr("width", sa.svg.w)
						.attr("height", sa.svg.h);
	/* Draw trace line */
	DrawTraceLine(sa, scoreMatrix, data);
	/* Draw the scores */
	DrawScores(sa, scoreMatrix);
	DrawRowNames(sa, scoreMatrix, object_alignment.seq1);
	DrawColNames(sa, scoreMatrix, object_alignment.seq2);
	DrawCells(sa, scoreMatrix, object_alignment.score); /*or max_score */
}

