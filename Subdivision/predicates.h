#pragma once
#include "Point.h"

// from robust predicates lib
void exactinit();
double orient2d(double* pa, double* pb, double* pc);
double incircle(double* pa, double* pb, double* pc, double* pd);

double orient2d(Point const& a, Point const& b, Point const& c);
double incircle(Point const& a, Point const& b, Point const& c, Point const& d);