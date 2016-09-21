/*
Copyright (C) 2016

Contact: Dmitry Sigaev <dima.sigaev@gmail.com>
*/

const unsigned char lal_encode31[] = {
  31,31,31,31, 31,31,31,31, 31,31,31,31, 31,31,31,31,
  31,31,31,31, 31,31,31,31, 31,31,31,31, 31,31,31,31,
  31,31,31,31, 31,31,29,31, 31,31,27,28, 31,26,26,31,
  31,31,31,31, 31,31,31,31, 31,31,31,31, 31,31,31,31,

  30, 0, 1, 2,  3, 4, 5, 6,  7, 8, 9,10, 11,12,13,14,
  15,16,17,18, 19,20,21,22, 23,24,25,31, 31,31,31,31,
  31, 0, 1, 2,  3, 4, 5, 6,  7, 8, 9,10, 11,12,13,14,
  15,16,17,18, 19,20,21,22, 23,24,25,31, 31,31,31,31,

  31,31,31,31, 31,31,31,31, 31,31,31,31, 31,31,31,31,
  31,31,31,31, 31,31,31,31, 31,31,31,31, 31,31,31,31,
  31,31,31,31, 31,31,31,31, 31,31,31,31, 31,31,31,31,
  31,31,31,31, 31,31,31,31, 31,31,31,31, 31,31,31,31,

  31,31,31,31, 31,31,31,31, 31,31,31,31, 31,31,31,31,
  31,31,31,31, 31,31,31,31, 31,31,31,31, 31,31,31,31,
  31,31,31,31, 31,31,31,31, 31,31,31,31, 31,31,31,31,
  31,31,31,31, 31,31,31,31, 31,31,31,31, 31,31,31,31 };


/*
* Contains the numerical representation of the complementary nucleotide symbol.
*/
const char cns[31] = { 19,21,6,7,26,26,2,3,26,26,12,26,10,13,26,26,26,24,18,0,0,1,22,13,17,26,26,26,26,26,26 };
const char lal_decode31_s[31] = { 'A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z','.','*','+','&','@' };


char identity_nuc[] = { " #  Matrix                                           \n \
    A    B    C    D    G    H    K    M    R    S    T    U    V    W    Y  \n \
A  1.0 -0.6 -0.6  0.2 -0.6  0.2 -0.6  0.6  0.6 -0.6 -0.6 -0.6  0.2  0.6 -0.6 \n \
B -0.6  0.2  0.2  0.2  0.2  0.2  0.2  0.2  0.2  0.2  0.2  0.2  0.2  0.2  0.2 \n \
C -0.6  0.2  1.0 -0.6 -0.6  0.2 -0.6  0.6 -0.6  0.6 -0.6 -0.6  0.2 -0.6  0.6 \n \
D  0.2  0.2 -0.6  0.2  0.2  0.2  0.2  0.2  0.2  0.2  0.2  0.2  0.2  0.2  0.2 \n \
G -0.6  0.2 -0.6  0.2  1.0 -0.6  0.6 -0.6  0.6  0.6 -0.6 -0.6  0.2 -0.6 -0.6 \n \
H  0.2  0.2  0.2  0.2 -0.6  0.2  0.2  0.2  0.2  0.2  0.2  0.2  0.2  0.2  0.2 \n \
K -0.6  0.2 -0.6  0.2  0.6  0.2  0.6 -0.6  0.6  0.6  0.6  0.6  0.2  0.6  0.6 \n \
M  0.6  0.2  0.6  0.2 -0.6  0.2 -0.6  0.6  0.6  0.6 -0.6 -0.6  0.2  0.6  0.6 \n \
R  0.6  0.2 -0.6  0.2  0.6  0.2  0.6  0.6  0.6  0.6 -0.6 -0.6  0.2  0.6 -0.6 \n \
S -0.6  0.2  0.6  0.2  0.6  0.2  0.6  0.6  0.6  0.6 -0.6 -0.6  0.2 -0.6  0.6 \n \
T -0.6  0.2 -0.6  0.2 -0.6  0.2  0.6 -0.6 -0.6 -0.6  1.0  1.0 -0.6  0.6  0.6 \n \
U -0.6  0.2 -0.6  0.2 -0.6  0.2  0.6 -0.6 -0.6 -0.6  1.0  1.0 -0.6  0.6  0.6 \n \
V  0.2  0.2  0.2  0.2  0.2  0.2  0.2  0.2  0.2  0.2 -0.6 -0.6  0.2  0.2  0.2 \n \
W  0.6  0.2 -0.6  0.2 -0.6  0.2  0.6  0.6  0.6 -0.6  0.6  0.6  0.2  0.6  0.6 \n \
Y -0.6  0.2  0.6  0.2 -0.6  0.2  0.6  0.6 -0.6  0.6  0.6  0.6  0.2  0.6  0.6 " };


