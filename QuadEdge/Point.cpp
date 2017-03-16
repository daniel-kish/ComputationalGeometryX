#include "Point.h"
#include <tuple>

bool operator==(Point const& p, Point const& q)
{
	return p.x == q.x && p.y == q.y; // TODO approx
}

bool operator!=(Point const& p, Point const& q)
{
	return !(p == q);
}

bool operator<(Point const& p, Point const& q)
{
	return std::make_tuple(p.x, p.y) < std::make_tuple(q.x, q.y);
}

Point operator+(Point const& p, Point const& q)
{
	return{p.x + q.x, p.y + q.y};
}

Point operator-(Point const& p, Point const& q)
{
	return{p.x - q.x, p.y - q.y};
}

Point operator*(Point const& p, double f)
{
	return{p.x*f, p.y*f};
}

double operator*(Point const& p, Point const& q)
{
	return p.x*q.x + p.y*q.y;
}

double crossProdZ(Point const & a, Point const & b)
{
	return a.x*b.y - b.x*a.y;
}

std::ostream& operator<< (std::ostream& os, Point const& p)
{
	return os << '(' << p.x << ' ' << p.y << ')';
}

double norm(Point const& p)
{
	return hypot(p.x, p.y);
}

double sqNorm(Point const& p)
{
	return p.x*p.x + p.y*p.y;
}

double dist(Point const& p, Point const& q)
{
	return norm(p - q);
}