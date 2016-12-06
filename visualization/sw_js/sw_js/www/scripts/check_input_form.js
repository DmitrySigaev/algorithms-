/* check_input_form.js
 * Copyright 2016 Dmitry Sigaev
 *
 * Released under the MIT license -- see MIT-LICENSE for details
 */
/* Function to handle the form submission */
function handleForm(event) {
	var seq_1 = document.getElementById("seq_1").value;
	var seq_10 = [];
	for (var i = 0 ; i < seq_1.length; i++) {
		seq_10.push(seq_1.charAt(i));
	}
	var seq_20 = [];
	var seq_2 = document.getElementById("seq_2").value;
	for (var i = 0, seq_20; i < seq_2.length; i++) {
		seq_20.push(seq_2.charAt(i));
	}
	var matrix = document.getElementById("matrix_ident").value;
	var gapOpen = parseFloat(document.getElementById("gapOpen").value);
	var gapExt = parseFloat(document.getElementById("gapExt").value);
	CalculateSWandDraw(seq_10.join(' '), seq_20.join(' '), matrix, gapOpen, gapExt);
	return false;
}

