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

char blosum100[] = { "#  Matrix made by matblas from blosum100_3.iij         \n \
#  * column uses minimum score                                               \n \
#  BLOSUM Clustered Scoring Matrix in 1/3 Bit Units                          \n \
#  Blocks Database = /data/blocks_5.0/blocks.dat                             \n \
#  Cluster Percentage: >= 100                                                \n \
#  Entropy =   1.4516, Expected =  -1.0948                                   \n \
   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *    \n \
A  8 -3 -4 -5 -2 -2 -3 -1 -4 -4 -4 -2 -3 -5 -2  1 -1 -6 -5 -2 -4 -2 -2 -10   \n \
R -3 10 -2 -5 -8  0 -2 -6 -1 -7 -6  3 -4 -6 -5 -3 -3 -7 -5 -6 -4 -1 -3 -10   \n \
N -4 -2 11  1 -5 -1 -2 -2  0 -7 -7 -1 -5 -7 -5  0 -1 -8 -5 -7  5 -2 -3 -10   \n \
D -5 -5  1 10 -8 -2  2 -4 -3 -8 -8 -3 -8 -8 -5 -2 -4 -10 -7 -8  6  0 -4 -10  \n \
C -2 -8 -5 -8 14 -7 -9 -7 -8 -3 -5 -8 -4 -4 -8 -3 -3 -7 -6 -3 -7 -8 -5 -10   \n \
Q -2  0 -1 -2 -7 11  2 -5  1 -6 -5  2 -2 -6 -4 -2 -3 -5 -4 -5 -2  5 -2 -10   \n \
E -3 -2 -2  2 -9  2 10 -6 -2 -7 -7  0 -5 -8 -4 -2 -3 -8 -7 -5  0  7 -3 -10   \n \
G -1 -6 -2 -4 -7 -5 -6  9 -6 -9 -8 -5 -7 -8 -6 -2 -5 -7 -8 -8 -3 -5 -4 -10   \n \
H -4 -1  0 -3 -8  1 -2 -6 13 -7 -6 -3 -5 -4 -5 -3 -4 -5  1 -7 -2 -1 -4 -10   \n \
I -4 -7 -7 -8 -3 -6 -7 -9 -7  8  2 -6  1 -2 -7 -5 -3 -6 -4  4 -8 -7 -3 -10   \n \
L -4 -6 -7 -8 -5 -5 -7 -8 -6  2  8 -6  3  0 -7 -6 -4 -5 -4  0 -8 -6 -3 -10   \n \
K -2  3 -1 -3 -8  2  0 -5 -3 -6 -6 10 -4 -6 -3 -2 -3 -8 -5 -5 -2  0 -3 -10   \n \
M -3 -4 -5 -8 -4 -2 -5 -7 -5  1  3 -4 12 -1 -5 -4 -2 -4 -5  0 -7 -4 -3 -10   \n \
F -5 -6 -7 -8 -4 -6 -8 -8 -4 -2  0 -6 -1 11 -7 -5 -5  0  4 -3 -7 -7 -4 -10   \n \
P -2 -5 -5 -5 -8 -4 -4 -6 -5 -7 -7 -3 -5 -7 12 -3 -4 -8 -7 -6 -5 -4 -4 -10   \n \
S  1 -3  0 -2 -3 -2 -2 -2 -3 -5 -6 -2 -4 -5 -3  9  2 -7 -5 -4 -1 -2 -2 -10   \n \
T -1 -3 -1 -4 -3 -3 -3 -5 -4 -3 -4 -3 -2 -5 -4  2  9 -7 -5 -1 -2 -3 -2 -10   \n \
W -6 -7 -8 -10 -7 -5 -8 -7 -5 -6 -5 -8 -4  0 -8 -7 -7 17  2 -5 -9 -7 -6 -10  \n \
Y -5 -5 -5 -7 -6 -4 -7 -8  1 -4 -4 -5 -5  4 -7 -5 -5  2 12 -5 -6 -6 -4 -10   \n \
V -2 -6 -7 -8 -3 -5 -5 -8 -7  4  0 -5  0 -3 -6 -4 -1 -5 -5  8 -7 -5 -3 -10   \n \
B -4 -4  5  6 -7 -2  0 -3 -2 -8 -8 -2 -7 -7 -5 -1 -2 -9 -6 -7  6  0 -4 -10   \n \
Z -2 -1 -2  0 -8  5  7 -5 -1 -7 -6  0 -4 -7 -4 -2 -3 -7 -6 -5  0  6 -2 -10   \n \
X -2 -3 -3 -4 -5 -2 -3 -4 -4 -3 -3 -3 -3 -4 -4 -2 -2 -6 -4 -3 -4 -2 -3 -10   \n \
* -10 -10 -10 -10 -10 -10 -10 -10 -10 -10 -10 -10 -10 -10 -10 -10 -10 -10 -10 -10 -10 -10 -10  1 " };


char gaptest1[] = { " #gaptest1.table                                          \n \
    A    B    C    D    G    H    K    M    R    S    T    U    V    W    Y    \n \
A  2.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0   \n \
B -1.0  2.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0   \n \
C -1.0 -1.0  2.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0   \n \
D -1.0 -1.0 -1.0  2.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0   \n \
G -1.0 -1.0 -1.0 -1.0  2.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0   \n \
H -1.0 -1.0 -1.0 -1.0 -1.0  2.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0   \n \
K -1.0 -1.0 -1.0 -1.0 -1.0 -1.0  2.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0   \n \
M -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0  2.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0   \n \
R -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0  2.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0   \n \
S -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0  2.0 -1.0 -1.0 -1.0 -1.0 -1.0   \n \
T -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0  2.0  1.0 -1.0 -1.0 -1.0   \n \
U -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0  1.0  2.0 -1.0 -1.0 -1.0   \n \
V -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0  2.0 -1.0 -1.0   \n \
W -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0  2.0 -1.0   \n \
Y -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0  2.0" };
