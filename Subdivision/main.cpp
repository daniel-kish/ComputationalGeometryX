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

struct Rect {
	Point origin;
	Point dir;
};
std::vector<Point> rectHull(Rect r, int x_pts, int y_pts)
{
	double wid = std::abs(r.dir.x);
	double height = std::abs(r.dir.y);
	double xStep = wid / x_pts;
	double yStep = height / y_pts;

	std::vector<Point> pts; pts.reserve(2 * x_pts + 2 * x_pts);
	for (double x = 0.0; x < wid; x += xStep)
		pts.push_back(Point{x,0.0});
	for (double y = 0.0; y < height; y += yStep)
		pts.push_back(Point{wid,y});
	for (double x = wid; x > 0.0; x -= xStep)
		pts.push_back(Point{x,height});
	for (double y = height; y > 0.0; y -= yStep)
		pts.push_back(Point{0.0,y});

	pts.shrink_to_fit();

	if (r.origin != Point{0.0,0.0})
		for (Point& p : pts)
			p = p + r.origin;

	return pts;
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

	Rect rect{{-5,-5},{3.00001,10.000023}};
	std::vector<Point> S = rectHull(rect, 50, 50);
	//std::vector<Point> S(200);
	//std::random_device rd;
	//std::mt19937 mt(rd());
	//std::uniform_real_distribution<double> d(-6,6);

	//std::generate(S.begin(), S.end(), [&]() {return Point{d(mt),d(mt)}; });


	//std::sort(S.begin(), S.end());

	Subdivision dt;
	Edge l,r;

	auto t1 = high_resolution_clock::now();

	std::tie(dt,l,r) = delaunay_dnc(S.begin(), S.end());
	
	auto t2 = high_resolution_clock::now();

	std::cout << duration_cast<milliseconds>(t2 - t1).count() << '\n';

	Graphics g;

	for (auto qiter = dt.edges.begin(); qiter != dt.edges.end(); ++qiter)
	{
		Edge e(qiter);
		g.addLine(Org(e)->point, Dest(e)->point);
	}
	for (Point const& p : S)
		g.addPoint(p);
	std::ofstream xml{"trian.xml"};
	g.output(xml);

	std::cout << dt.edges.size() << '\n';
}