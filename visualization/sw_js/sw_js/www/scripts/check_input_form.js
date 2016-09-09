/* Function to handle the form submission */
function handleForm(event){
    //console.log(document.getElementById("seq_q").value)
    //draw(document.getElementById("myVal").value)

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
    var gapOpen =
        parseFloat(document.getElementById("gapOpen").value);
    var gapExt =
        parseFloat(document.getElementById("gapExt").value);

    Calculate(seq_10.join(' '), seq_20.join(' '), matrix, gapOpen, gapExt);
    return false;
}