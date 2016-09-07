#include <utility>   // for std::pair
#include <iterator>  // for std::iterator_traits
#include <iostream>  // for std::cout
#include <ostream>   // for output operator and std::endl
#include <algorithm> // for std::copy
#include <iterator>  // for std::output_iterator
#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <future>
 
// Function template max_subseq
//
// Given a sequence of integers, find a subsequence which maximizes
// the sum of its elements, that is, the elements of no other single
// subsequence add up to a value larger than this one.
//
// Requirements:
// * ForwardIterator is a forward iterator
// * ForwardIterator's value_type is less-than comparable and addable
// * default-construction of value_type gives the neutral element
//   (zero)
// * operator+ and operator< are compatible (i.e. if a>zero and
//   b>zero, then a+b>zero, and if a<zero and b<zero, then a+b<zero)
// * [begin,end) is a valid range
//
// Returns:
//   a pair of iterators describing the begin and end of the
//   subsequence
template<typename ForwardIterator>
 std::pair<ForwardIterator, ForwardIterator>
 max_subseq(ForwardIterator begin, ForwardIterator end)
{
  typedef typename std::iterator_traits<ForwardIterator>::value_type
    value_type;
 
  ForwardIterator seq_begin = begin, seq_end = seq_begin;
  value_type seq_sum = value_type();
  ForwardIterator current_begin = begin;
  value_type current_sum = value_type();
 
  value_type zero = value_type();
 
  for (ForwardIterator iter = begin; iter != end; ++iter)
  {
    value_type value = *iter;
    if (zero < value)
    {
      if (current_sum < zero)
      {
        current_sum = zero;
        current_begin = iter;
      }
    }
    else
    {
      if (seq_sum < current_sum)
      {
        seq_begin = current_begin;
        seq_end = iter;
        seq_sum = current_sum;
      }
    }
    current_sum += value;
  }
 
  if (seq_sum < current_sum)
  {
    seq_begin = current_begin;
    seq_end = end;
    seq_sum = current_sum;
  }
 
  return std::make_pair(seq_begin, seq_end);
}
 
// the test array
int array[] = { -1, -2, 3, 5, 6, -2, -1, 4, -4, 2, -1 };
 

using namespace std;

struct S1 {};
struct S2 {};

template<typename X = S1>
struct T;

template<>
struct T<S1> {
	static void f() { cout << 1 << endl; }
};

template<>
struct T<S2> {
	static void f() { cout << 2 << endl; }
};

template< template<typename X = S2> class C>
void g() {
	C<>::f();
}

#include <iostream>
uint8_t swap_bits(uint8_t byte_to_reverse)
{
	uint8_t b = byte_to_reverse;
	b = ((b * 0x0802LU & 0x22110LU) | (b * 0x8020LU & 0x88440LU)) * 0x10101LU >> 16;
	return b;
}
#define N 10
int32_t average(int32_t m[N])
{
	int32_t sum = 0;
	for (size_t i = 0; i < N; i++) {
		sum += m[i];
	}
	return sum / N;
}

int main()
{

	int32_t n[] = { 1,2,3,4,5,6,7,8,9,10 };
	int32_t av = average(n);
	float  a = 1 / 3., b = -1. / 7.;
	float c = 3.;
	float  i = 1 / 3., j = 1. / 7.;
	float k = 3.;

	if ((a + b)*c != (a*c + b*c))
	{
		printf("Distributivity is not always true");
	}

	if ((i + j)*k == (i*k + j*k))
	{
		printf("Distributivity is true");
	}




}


/*
int main() {
	g<T>();
	return 0;
}
*/
 /*
int main()
{
  // find the subsequence
  std::pair<int*, int*> seq = max_subseq(array, end(array));
 
  // output it
  std::copy(seq.first, seq.second, std::ostream_iterator<int>(std::cout, " "));
  std::cout << std::endl;
 
  return 0;
}
*/