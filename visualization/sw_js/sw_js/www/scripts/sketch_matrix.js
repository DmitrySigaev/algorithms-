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
	var color = "blue";
	return { svg: { w: w, h: h, margin: { left: left_margin, top: top_margin } }, cell: { w: cell_width, h: cell_height, padding: { row: row_padding, col: col_padding }, color: color } };
}

function sketch_matrixes(object_alignment) {
	var symSize = getMaxSymbolsSize(object_alignment.seq1, "Arial", 12);
	var svgArea = setSvgArea(object_alignment.seq1, object_alignment.seq2, symSize);
	console.log(svgArea);
	// Calculate the matrix. Imported from sw.js

	var score_mat = object_alignment.scorematrix;
	var trace_mat = object_alignment.directions;


	var data = new Array();
	var c = 0;
	var max_score = 0;
	for (i = 0; i < object_alignment.seq1.length; i++) { /* var nrows = yseq.length = object_alignment.seq1.length */
		for (j = 0; j < object_alignment.seq2.length; j++) { /* var ncols = xseq.length = object_alignment.seq2 */
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
}

