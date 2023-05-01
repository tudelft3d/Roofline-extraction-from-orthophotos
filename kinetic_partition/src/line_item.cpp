#include "line_item.h"



LineItem::LineItem()
{
	x1 = y1 = x2 = y2 = 0;
	r = g = b = 0;
}


LineItem::LineItem(double & _x1, double & _y1, double & _x2, double & _y2, uchar _r, uchar _g, uchar _b)
{
	x1 = _x1;
	y1 = _y1;
	x2 = _x2;
	y2 = _y2;

	r = _r;
	g = _g;
	b = _b;
}


LineItem::~LineItem()
{
}
