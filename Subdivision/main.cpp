#include <iostream>
#include <list>
#include "Point.h"
#include <tuple>
#include <vector>
#include <fstream>
#include <sstream>
#include <random>
#include <chrono>
#include "Subdivision.h"
#include "predicates.h"
#include "delaunay.h"
#include "boost/math/constants/constants.hpp"

struct Rect {
	Point origin;
	Point dir;
};
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

std::vector<Point> circleHull(Point cen, double rad, int N = 20)
{
	using boost::math::double_constants::pi;

	std::vector<Point> hull;
	hull.reserve(N);
	double step = 2.0*pi/N;

	for (int i = 0; i < N; ++i)
	{
		double phi{step*i};
		Point p{rad*cos(phi),rad*sin(phi)};
		hull.push_back(p);
	}
	return hull;
}

struct Graphics
{
	double scale{30};
	double wid, height;
	double X0, Y0;
	std::string main_code;
	std::string ending;
	double pointRad{2.0};
	double lineWidth{1.0};
	Graphics(double w = 1000.0, double h = 500.0, double X = 0.0, double Y = 0.0)
		: X0{X}, Y0{Y}, wid{w}, height{h}
	{
		main_code =
			R"(<?xml version="1.0" encoding="UTF-8" standalone="no"?>)"
			"\n"
			R"(<svg version = "1.1" )" "\n"
			R"(baseProfile="full" )" "\n"
			R"(xmlns = "http://www.w3.org/2000/svg" )" "\n"
			R"(xmlns:xlink = "http://www.w3.org/1999/xlink" )" "\n"
			R"(xmlns:ev = "http://www.w3.org/2001/xml-events" )" "\n"
			R"(height = "1000px"  width = "1200px">)" "\n";

		std::ostringstream tmp;
		tmp << R"(<rect x=")" << X0 << R"(" y=")" << Y0 << R"(" width=")" << wid
			<< R"(" height=")" << height << R"(" fill="none" stroke="black" stroke-width="5px"/>)";
		main_code.append(tmp.str() + "\n");
		tmp.str(""); tmp.clear();
		tmp << R"(<g transform="translate)" << '(' << wid / 2 << ',' << height / 2 << ')'
			<< R"q( scale)q" << '(' << scale << ',' << -scale << ')' << R"q("> )q";
		main_code.append(tmp.str() + "\n");
		ending = "\n</g>\n</svg>";
	}

	void addLine(Point const& p, Point const& q)
	{
		std::ostringstream str;
		str << R"( <line x1=")" << p.x << R"(" y1=")" << p.y
			<< R"(" x2=")" << q.x << R"(" y2=")" << q.y
			<< R"(" style="stroke:black;stroke-width:)" << lineWidth / scale << R"("/> )"
			<< '\n';
		main_code.append(str.str());
	}
	void addPoint(Point const& p)
	{
		std::ostringstream str;
		str << R"(<circle cx=")" << p.x
			<< R"(" cy=")" << p.y << R"(" r=")" << pointRad / scale << R"("  fill="black" />)"
			<< '\n';
		main_code.append(str.str());
	}
	void output(std::ofstream& svgfile)
	{
		svgfile << main_code << ending;
	}
};

int main()
{
	using namespace std::chrono;
	exactinit();

	std::vector<Point> points = rectHull(Rect{{-5,-5},{10,10}}, 10, 10);

	//std::random_device rd;
	//std::mt19937 mt{rd()};
	//std::uniform_real_distribution<double> dist(-4, 4);

	//std::generate_n(points.begin(), points.size(), [&]() {
	//	return Point{dist(mt),dist(mt)};
	//});
	//std::sort(points.begin(), points.end());

	
	Subdivision dt;
	EdgeRef l, r;
	std::tie(dt, l, r) = delaunay_dnc(points.begin(), points.end());

	VertexRef v = insertSite(dt, {0,0.1});
	if (v == dt.vertices.end())
		std::cerr << "404 not found\n";
	else {
		EdgeRef iter = v->leaves;
		EdgeRef first{iter};
		do {
			std::cout << Org(iter)->point << ' ' << Dest(iter)->point << '\n';
			iter = iter.Onext();
		} while (iter != first);
	}

	Graphics g;
	for (Point const& p : points)
		g.addPoint(p);

	for (auto qref = dt.edges.begin(); qref != dt.edges.end(); ++qref)
	{
		EdgeRef e(qref);
		g.addLine(Org(e)->point, Dest(e)->point);
	}

	std::ofstream xml{"trian.xml"};
	g.output(xml);
}