#pragma once
#include "Point.h"
#include <vector>
#include <array>

struct Rect {
	Point origin;
	Point dir;
};

struct Circle {
	Point center;
	double radius;
};

std::vector<Point> rectHull(Rect r, int Nx, int Ny);
std::vector<Point> circleHull(Point cen, double rad, int N = 20);
std::vector<Point> rectUniform(Rect, int);
std::array<Point, 3> triangleCover(std::vector<Point>&);

double quality_measure(Point const& p1, Point const& p2, Point const& p3);
Point circumCenter(Point const& a, Point const& b, Point const& c);
Circle circumCircle(Point const& a, Point const& b, Point const& c);
double circumRadius(Point const& p1, Point const& p2, Point const& p3);