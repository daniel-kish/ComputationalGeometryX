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
#include "geom.h"
#include "mesh.h"
#include <valarray>

double step = 1.0 / 5.0;

double red (double t)
{
	while (t > 1.0) t = -1.0;
	while (t < 0.0) t += 1.0;
	if (t < step)
		return 1.0;
	else if (t < 2 * step)
		return -2 * t + 2;
	else if (t < 4 * step)
		return 0.0;
	else if (t < 5 * step)
		return 5 * t - 4;
}

double green (double t)
{
	return 1 - red(t + step);
}

double blue (double t)
{
	return 1 - red(t - step);
}

struct Graphics
{
	std::random_device rdev;
	std::mt19937 mt;
	double scale{120};
	double wid, height;
	double X0, Y0;
	std::string main_code;
	std::string ending;
	double pointRad{1.0};
	double lineWidth{0.5};
	Graphics(double w = 1000.0, double h = 500.0, double X = 0.0, double Y = 0.0)
		: X0{X}, Y0{Y}, wid{w}, height{h}, mt{rdev()}
	{
		main_code =
			R"(<?xml version="1.0" encoding="UTF-8" standalone="no"?>)"
			"\n"
			R"(<svg version = "1.1" )" "\n"
			R"(baseProfile="full" )" "\n"
			R"(xmlns = "http://www.w3.org/2000/svg" )" "\n"
			R"(xmlns:xlink = "http://www.w3.org/1999/xlink" )" "\n"
			R"(xmlns:ev = "http://www.w3.org/2001/xml-events" )" "\n"
			R"(height = "5000px"  width = "5000px">)" "\n";

		std::ostringstream tmp;
		/*tmp << R"(<rect x=")" << X0 << R"(" y=")" << Y0 << R"(" width=")" << wid
			<< R"(" height=")" << height << R"(" fill="none" stroke="black" stroke-width="5px"/>)";*/
		main_code.append(tmp.str() + "\n");
		tmp.str(""); tmp.clear();
		tmp << R"(<g transform="translate)" << '(' << wid / 2 + X0 << ',' << height / 2 + Y0 << ')'
			<< R"q( scale)q" << '(' << scale << ',' << -scale << ')' << R"q("> )q";
		main_code.append(tmp.str() + "\n");
		ending = "\n</g>\n</svg>";
	}

	void addLine(Point const& p, Point const& q, std::string color = {}, double width=0.0)
	{
		if (width == 0.0) width = lineWidth;
		if (color.empty()) color = "black";
		std::ostringstream str;
		str << R"( <line x1=")" << p.x << R"(" y1=")" << p.y
			<< R"(" x2=")" << q.x << R"(" y2=")" << q.y
			<< R"(" style="stroke:)" << color << R"(;stroke-width:)" << width / scale << R"("/> )"
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
	void add_polygon(Point a, Point b, Point c, std::string color = std::string{}, double opacity = 1.0)
	{
		std::ostringstream str;
		int R = 255, G = 255, B = 255;
		if (color.empty()) {
			std::uniform_int_distribution<int> dist_int(0, 255);
			R = dist_int(mt);
			G = dist_int(mt);
			B = dist_int(mt);
			str << "<polygon points = \"" <<
				a.x << ',' << a.y << ' ' <<
				b.x << ',' << b.y << ' ' <<
				c.x << ',' << c.y << "\" style=\"fill: rgb("
				<< R << ',' << G << ',' << B << "); opacity:" << opacity << "; stroke-width:0\" />\n";
		}
		else {
			str << "<polygon points = \"" <<
				a.x << ',' << a.y << ' ' <<
				b.x << ',' << b.y << ' ' <<
				c.x << ',' << c.y << "\" style=\"fill:" << color
				<< "; opacity:" << opacity << "; stroke-width:0\" />\n";
		}
		main_code.append(str.str());
	}

	void add_polygon(Point a, Point b, Point c, int R, int G, int B, double opacity = 1.0)
	{
		std::ostringstream str;
		str << "<polygon points = \"" <<
			a.x << ',' << a.y << ' ' <<
			b.x << ',' << b.y << ' ' <<
			c.x << ',' << c.y << "\" style=\"fill: rgb("
			<< R << ',' << G << ',' << B << "); opacity:" << opacity << "; stroke-width:0\" />\n";
		main_code.append(str.str());
	}

	void output(std::ofstream& svgfile)
	{
		svgfile << main_code << ending;
	}
};

std::valarray<double> scale(double t)
{
	std::valarray<double> v{red(-t),green(-t),blue(-t)};
	return v * 255.0;
}

int main()
{
	using namespace std::chrono;
	using namespace std::literals;
	exactinit();

	/*std::vector<Point> model{
		{-4,4},{-4,1},{-3,1},{-3,-1},{-4,-1},{-4,-2},{3,-2},{4,-1},{4,2}
	};*/
	std::vector<Point> model = rectHull({{-5,-5},{10,10}}, 1, 1);
	
	auto trian = triangleCover(model);
	std::vector<Point> cover(trian.begin(), trian.end());
	std::sort(cover.begin(), cover.end());
	Subdivision dt;
	EdgeRef l, r;
	std::tie(dt, l, r) = delaunay_dnc(cover.begin(), cover.end());


	insertClosedLoop(dt, model);

	std::vector<Point> hole = anyHull([](double phi) {return 1.0+sin(phi)*sin(phi); }, 100);
	insertClosedLoop(dt, hole);

	hole = circleHull({4,4}, 0.5, 50);
	insertClosedLoop(dt, hole);

	hole = circleHull({-4,4}, 0.5, 50);
	insertClosedLoop(dt, hole);

	hole = circleHull({-4,-4}, 0.5, 50);
	insertClosedLoop(dt, hole);

	hole = circleHull({4,-4}, 0.5, 50);
	insertClosedLoop(dt, hole);

	init_faces(dt);
	mark_outer_faces(dt, r);


	FaceRef face; double max_ratio, max_area, min_area;
	std::tie(face, max_ratio) = find_worst(dt);
	std::tie(face, max_area) = find_biggest(dt);
	std::tie(face, min_area) = find_smallest(dt);
	std::cout << "before: " << max_ratio << ' ' << max_area << '\n';

	ruppert_refinement(dt, 1.0, 3000);


	double ratio_after, area_after;
	std::tie(face, ratio_after) = find_worst(dt);
	std::tie(face, area_after) = find_biggest(dt);
	std::cout << "after: " << ratio_after << ' ' << area_after << '\n';
	

	// connectivity checks
	for (auto v = dt.vertices.begin(); v != dt.vertices.end(); ++v)
	{
		if (Org(v->leaves) != v) {
			std::cerr << "vertex connectivity check failed:\n";
			std::exit(1);
		}
		auto eiter = v->leaves;
		auto end = eiter;
		do {
			if (Org(eiter) != v) {
				std::cerr << "Onext data fatal\n";
				std::exit(1);
			}
			eiter = eiter.Onext();
		} while (eiter != end);
	}
	std::cerr << "vertex connectivity check OK\n";
	for (auto f = dt.faces.begin(); f != dt.faces.end(); ++f)
	{
		auto edge = f->bounds;
		auto face = Left(edge);
		if (Left(f->bounds) != f) {
			std::cerr << "faces connectivity check failed\n";
			std::exit(1);
		}
		auto eiter = f->bounds;
		auto end = eiter;
		do {
			if (Left(eiter) != f) {
				std::cerr << "Lnext data fatal\n";
				std::exit(1);
			}
			eiter = eiter.Lnext();
		} while (eiter != end);
	}
	std::cerr << "faces connectivity check OK\n";

	Graphics g(1000,500,300,400);

	for (auto face = dt.faces.begin(); face != dt.faces.end(); ++face)
	{
		if (face == dt.outer_face) continue;
		auto e = face->bounds;
		if (!face->mark) {
			//g.add_polygon(Org(e)->point, Dest(e)->point, Dest(e.Onext())->point, "green", 0.5);
		}
		else {
			//double ratio = quality_measure(Org(e)->point, Dest(e)->point, Dest(e.Onext())->point);
			//auto col = scale(ratio*0.8 / max_ratio);
			//g.add_polygon(Org(e)->point, Dest(e)->point, Dest(e.Onext())->point,
			//	col[0], col[1], col[2]);
			g.add_polygon(Org(e)->point, Dest(e)->point, Dest(e.Onext())->point, "blue",0.25);
		}
	}
	for (auto qref = dt.edges.begin(); qref != dt.edges.end(); ++qref)
	{
		EdgeRef e(qref);
		if (e.data().var.which() == 1) {
			e = e.Rot();
		}
		if (e.data().fixed) {
			g.addLine(Org(e)->point, Dest(e)->point, "red"s, 2.5);
		}
		if (Left(e)->mark || Right(e)->mark)
			g.addLine(Org(e)->point, Dest(e)->point);
	}
	for (Subdivision::Vertex const& p : dt.vertices)
		g.addPoint(p.point);

	std::ofstream xml{"trian.xml"};
	g.output(xml);

	std::cout << "Euler invariant: " << dt.vertices.size() - dt.edges.size() + dt.faces.size() << '\n';
}

/*std::vector<Point> hole = {
{-1,-1}, {1,-1}, {0.5, 0}, {1, 1}, {-1,1}
};*/