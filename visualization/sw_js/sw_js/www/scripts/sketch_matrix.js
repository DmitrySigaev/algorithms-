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

function sketch_matrixes(object_alignment) {
	var symSize = getMaxSymbolsSize(object_alignment.seq1, "Arial", 12);
	console.log(symSize);
}

