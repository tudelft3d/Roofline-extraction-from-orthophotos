#pragma once
#include "defs.h"


class LineItem
{
public:
	LineItem();

	LineItem(double & _x1, double & _y1, double & _x2, double & _y2, uchar _r, uchar _g, uchar _b);

	~LineItem();

public:
	double x1;
	double y1;
	double x2;
	double y2;

	uchar r;
	uchar g;
	uchar b;
};

