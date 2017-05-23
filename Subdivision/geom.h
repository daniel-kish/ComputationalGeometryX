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


template <class Fun>
std::vector<Point> anyHull(Fun& f, int N)
{
	using boost::math::double_constants::pi;

	std::vector<Point> hull;
	hull.reserve(N);
	double step = 2.0*pi / N;

	for (int i = 0; i < N; ++i)
	{
		double phi{step*i};
		double rad = f(phi);
		Point p{rad*cos(phi),rad*sin(phi)};
		hull.push_back(p);
	}
	return hull;
}


double quality_measure(Point const& p1, Point const& p2, Point const& p3);
Point circumCenter(Point const& a, Point const& b, Point const& c);
Circle circumCircle(Point const& a, Point const& b, Point const& c);
double circumRadius(Point const& p1, Point const& p2, Point const& p3);

template <class T, std::size_t N>
T KahanSum(std::array<T,N> const& v)
{
	T sum{};
	T c{};
	const auto size = v.size();
	for (std::size_t i{}; i < size; ++i)
	{
		T y = v[i] - c;
		T t = sum + y;
		c = (t - sum) - y;
		sum = t;
	}
	return sum;
}

double triangleArea(Point const& p1, Point const& p2, Point const& p3);

Point off_center(Point a, Point b, Point c, double beta);