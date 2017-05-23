#include "geom.h"
#include "Point.h"
#include "boost\math\constants\constants.hpp"
#include <algorithm>
#include <random>
#include "predicates.h"

std::vector<Point> rectHull(Rect r, int Nx, int Ny)
{
	std::vector<Point> hull;

	double W = r.dir.x;
	double H = r.dir.y;

	double x0 = r.origin.x;
	double y0 = r.origin.y;
	double xend = x0 + W;
	double yend = y0 + H;
	double dx = W / Nx;
	double dy = H / Ny;

	hull.push_back(Point{x0,y0}); // SW

	for (int i = 1; i <= Nx - 1; ++i)
		hull.push_back(Point{x0 + dx*i,y0});

	hull.push_back(Point{xend,y0}); // SE

	for (int i = 1; i <= Ny - 1; ++i)
		hull.push_back(Point{xend,y0 + dy*i});

	hull.push_back(Point{xend,yend}); // NE

	for (int i = 1; i <= Nx - 1; ++i)
		hull.push_back(Point{xend - dx*i,yend});

	hull.push_back(Point{x0,yend}); // NW

	for (int i = 1; i <= Ny - 1; ++i)
		hull.push_back(Point{x0, yend - dy*i});

	return hull;
}

std::vector<Point> circleHull(Point cen, double rad, int N)
{
	using boost::math::double_constants::pi;

	std::vector<Point> hull;
	hull.reserve(N);
	double step = 2.0*pi / N;

	for (int i = 0; i < N; ++i)
	{
		double phi{step*i};
		Point p{rad*cos(phi),rad*sin(phi)};
		p = p + cen;
		hull.push_back(p);
	}
	return hull;
}

std::vector<Point> rectUniform(Rect r, int N)
{
	std::random_device rd;
	std::mt19937 mt{0/*rd()*/};
	std::uniform_real_distribution<double> x_dist(r.origin.x, r.origin.x + r.dir.x);
	std::uniform_real_distribution<double> y_dist(r.origin.y, r.origin.y + r.dir.y);

	std::vector<Point> points(N);
	std::generate_n(points.begin(), points.size(), [&]() {
		return Point{x_dist(mt),y_dist(mt)};
	});
	return points;
}

std::array<Point, 3> triangleCover(std::vector<Point>& pts)
{
	
	auto lr = std::minmax_element(pts.begin(), pts.end(), [](Point const& p, Point const& q) {
		return p.x < q.x;
	});
	double x_min{lr.first->x}, x_max{lr.second->x};

	auto bt = std::minmax_element(pts.begin(), pts.end(), [](Point const& p, Point const& q) {
		return p.y < q.y;
	});
	double y_min{bt.first->y}, y_max{bt.second->y};

	Rect bounds{{x_min,y_min},{x_max - x_min, y_max - y_min}};

	// inscribed circle
	Point center = bounds.origin + bounds.dir*0.5;
	double rad = dist(bounds.origin, center);
	double a = 2.0*sqrt(3.0)*rad;
	Point D = center - Point{0,rad};
	
	Point A = D - Point{0.5*a,0.0};
	Point B = D + Point{0.5*a,0.0};
	
	double h = 3.0*rad;
	Point C = D + Point{0.0,h};

	Point barycenter = (A + B + C) * (1.0/3.0);

	std::array<Point, 3> trian{A,B,C};

	for (Point& p : trian)
		p = p + (p - barycenter)*0.2;

	return trian;
}

Point circumCenter(Point const& a, Point const& b, Point const& c)
{
	double x = sqNorm(a)*(b.y - c.y) + sqNorm(b)*(c.y - a.y) + sqNorm(c)*(a.y - b.y);
	double y = sqNorm(a)*(c.x - b.x) + sqNorm(b)*(a.x - c.x) + sqNorm(c)*(b.x - a.x);
	double D = 2.0f * (a.x*(b.y - c.y) + b.x*(c.y - a.y) + c.x*(a.y - b.y));

	return Point{x / D, y / D};
}

Circle circumCircle(Point const& a, Point const& b, Point const& c)
{
	Point p = circumCenter(a, b, c);
	return Circle{p,dist(p,a)};
}

double circumRadius(Point const& p1, Point const& p2, Point const& p3)
{
	double a = norm(p1 - p2);
	double b = norm(p2 - p3);
	double c = norm(p3 - p1);

	double d = a*b*c / sqrt((a + b + c)*(-a + b + c)*(a - b + c)*(a + b - c));
	return d;
}

double quality_measure(Point const& p1, Point const& p2, Point const& p3)
{
	double a = norm(p1 - p2);
	double b = norm(p2 - p3);
	double c = norm(p3 - p1);
	double min_edge_len = std::min({a,b,c});

	double circumradius = a*b*c / sqrt((a + b + c)*(-a + b + c)*(a - b + c)*(a + b - c));

	return circumradius / min_edge_len;
}

double triangleArea(Point const& p1, Point const& p2, Point const& p3)
{
	double a = norm(p1 - p2);
	double b = norm(p2 - p3);
	double c = norm(p3 - p1);

	std::array<double, 3> lens{a,b,c};
	double  p = KahanSum(lens) / 2.0;
	double area = sqrt(p*(p-a)*(p-b)*(p-c));
	
	return area;
}

Point off_center(Point A, Point B, Point C, double beta_req)
{
	double det = orient2d(A, B, C);
	// precondition: C lies on {A,B} bisector
	double l = dist(A, C);
	double L = dist(A, B);
	double h = dist(C, (A + B) * 0.5);
	double R = 0.5*l*l / h;
	double beta = R / L;

	if (beta <= beta_req) {
		//std::cout << "C\n";
		return C;
	}
	beta_req *= 0.9;
	double Rnew = L*beta_req;
	double a = 1;
	double b = -2 * Rnew;
	double c = L*L / 4;

	double disc = b*b - 4*a*c;
	double newh1 = (-b + sqrt(disc))*0.5;
	double newh2 = (-b - sqrt(disc))*0.5;
	double newh = std::max({newh1,newh2});

	double t = 1 - newh / h;
	if (t < 0 || t > 1) {
		std::cerr << "t out of range\n";
		std::exit(1);
	}
	Point CM = (A + B)*0.5 - C;
	Point newC = C + CM*t;

	// check:
	double newl = dist(newC, A);
	double newR = 0.5*newl*newl / newh;
	double newbeta = newR / L;
	//std::cout << "newC\n";
	return newC;
}