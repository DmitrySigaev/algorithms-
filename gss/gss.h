/*
Copyright (C) 2016 Dmitry Sigaev

Contact: Dmitry Sigaev <dima.sigaev@gmail.com>
*/


#ifndef _GSS_H_
#define _GSS_H_


typedef struct tag_range{
	int start, end, sum;
}range_t;

range_t maxSubseq(const int sequence[], const int len);

#endif /* _GSS_H_ */