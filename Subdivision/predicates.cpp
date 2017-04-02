#include "predicates.h"
#include "Point.h"

double orient2d(Point const& a, Point const& b, Point const& c) {
	return orient2d((double*)&a, (double*)&b, (double*)&c);
}

double incircle(Point const& a, Point const& b, Point const& c, Point const& d) {
	return incircle((double*)&a, (double*)&b, (double*)&c, (double*)&d);
}