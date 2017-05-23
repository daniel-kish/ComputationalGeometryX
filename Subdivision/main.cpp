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
	double scale{100};
	double wid, height;
	double X0, Y0;
	std::string main_code;
	std::string ending;
	double pointRad{2.0};
	double lineWidth{2.0};
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

	void addLine(Point const& p, Point const& q, std::string color = {}, double width = 0.0)
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

	void addPoint(Point const& p, double rad = 0.0)
	{
		if (rad == 0.0) rad = pointRad;
		std::ostringstream str;
		str << R"(<circle cx=")" << p.x
			<< R"(" cy=")" << p.y << R"(" r=")" << rad / scale << R"("  fill="black" />)"
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

	void add_polygon(FaceRef face, std::string color = std::string{}, double opacity = 1.0)
	{
		std::ostringstream str;
		str << "<polygon points = \"";
		auto ei = face->bounds;
		auto end = ei;
		do {
			Point const& p = Org(ei)->point;
			str << p.x << ',' << p.y << ' ';
			ei = ei.Lnext();
		} while (ei != end);
		str << "\" style=\"fill:" << color
			<< "; opacity:" << opacity << "; stroke-width:0\" />\n";
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

EdgeRef insertEdge2(Subdivision& s, VertexRef a, VertexRef b)
{
	EdgeRef e = a->leaves;

	EdgeRef p = e;
	do {
		if (Dest(e) == b) {
			e.data().fixed = true;
			e.Sym().data().fixed = true;
			return e;
		}
		e = e.Onext();
	} while (e != p);

	p = e;
	do {
		if (rightOf(Dest(e), {a,b}) && leftOf(Dest(e.Onext()), {a,b}))
			break;
		e = e.Onext();
	} while (e != p);

	assert(Dest(e.Onext()) == Dest(e.Lnext()));

	EdgeRef first = e.Onext().Sym();
	assert(Dest(first) == a);

	e = e.Lnext();
	p = e;
	do {
		if (Dest(e.Oprev()) == b) {
			s.deleteEdge(e);
			break;
		}
		VertexRef n = Dest(e.Oprev());
		double det = orient2d(a->point, b->point, n->point);
		if (det > 0.0) // leftOf
			p = e.Oprev();
		else if (det < 0.0) // rightOf
			p = e.Dnext();
		s.deleteEdge(e);
		e = p;
	} while (true);

	auto ei = first;
	do
		ei = ei.Lnext();
	while (Dest(ei) != b);

	EdgeRef last = ei;

	auto c = s.connect(first, last.Lnext());

	c.data().fixed = true;
	c.Sym().data().fixed = true;

	triangulatePseudoPolygon(s, c);
	triangulatePseudoPolygon(s, c.Sym());

	return c;
}

int nr_triangles(double q)
{
	std::vector<Point> model{
		{-4,2},{-4,1},{-3,1},{-3,-1},{-4,-1},{-4,-2},{3,-2},{4,-1},{4,2}
	};

	auto trian = triangleCover(model);
	std::vector<Point> cover(trian.begin(), trian.end());
	std::sort(cover.begin(), cover.end());
	Subdivision dt;
	EdgeRef l, r;
	std::tie(dt, l, r) = delaunay_dnc(cover.begin(), cover.end());

	insertClosedLoop(dt, model);

	std::vector<Point> hole = {
		{-1,-1},{1,-1},{0.5, 0},{1, 1},{-1,1}
	};
	insertClosedLoop(dt, hole);

	insertClosedLoop(dt, circleHull({2.5,0}, 0.25, 20));

	auto a = insertSite(dt, {-3.5,1.5});
	auto b = insertSite(dt, {3.5,1.5});
	insertEdge(dt, a, b);


	/*  ONLY USE *_WF FROM HERE ON */
	init_faces(dt);
	mark_outer_faces(dt, r);

	chew_2nd_refinement_alper(dt, q, 1.0, 0.5);

	int trs = std::count_if(dt.faces.begin(), dt.faces.end(), [](Subdivision::Face const& f) {
		return f.mark == 1;
	});

	return trs;
}

int main()
{
	using namespace std::chrono;
	using namespace std::literals;

	exactinit();

	//std::vector<Point> model{
	//{-4,2},{-4,1},{-3,1},{-3,-1},{-4,-1},{-4,-2},{3,-2},{4,-1},{4,2}
	//};
	auto model = rectHull({{-4,-4},{8,8}}, 1, 1);

	auto trian = triangleCover(model);
	std::vector<Point> cover(trian.begin(), trian.end());
	std::sort(cover.begin(), cover.end());
	Subdivision dt;
	EdgeRef l, r;
	std::tie(dt, l, r) = delaunay_dnc(cover.begin(), cover.end());

	insertClosedLoop(dt, model);

	std::vector<Point> hole = {
	{-1,-1}, {1,-1}, {0.5, 0}, {1, 1}, {-1,1}
	};

	//insertClosedLoop(dt, hole);
	insertClosedLoop(dt, circleHull({0,0}, 1.0, 20));

	insertClosedLoop(dt, circleHull({3.5,3.5}, 0.25, 10));
	insertClosedLoop(dt, circleHull({-3.5,3.5}, 0.25, 10));
	insertClosedLoop(dt, circleHull({-3.5,-3.5}, 0.25, 10));
	insertClosedLoop(dt, circleHull({3.5,-3.5}, 0.25, 10));

	/*  ONLY USE *_WF FROM HERE ON */
	init_faces(dt);
	mark_outer_faces(dt, r);

	FaceRef face; double max_ratio, max_area, min_area;
	std::tie(face, max_ratio) = find_worst(dt);
	std::tie(face, max_area) = find_biggest(dt);
	std::tie(face, min_area) = find_smallest(dt);
	std::cout << "before: " << max_ratio << ' ' << max_area << '\n';


	chew_2nd_refinement_alper(dt, 0.895);
	//chew_2nd_refinement(dt, 0.895);
	//ruppert_refinement(dt, 0.895);
	/*for (auto face = dt.faces.begin(); face != dt.faces.end(); ++face)
	{
		if (!face->mark || face == dt.outer_face) continue;
		auto e = face->bounds;
		double ratio = quality_measure(Org(e)->point, Dest(e)->point, Dest(e.Onext())->point);
		if (ratio > 1.0) {
			auto ei = e;
			double min_l=std::numeric_limits<double>::max();
			EdgeRef min_e = e;
			do {
				double l = dist(Org(ei)->point, Dest(ei)->point);
				if (l < min_l) {
					min_e = ei;
					min_l = l;
				}
				ei = ei.Lnext();
			} while (ei != e);
			Point C = circumCenter(Org(e)->point, Dest(e)->point, Dest(e.Onext())->point);
			Point off = off_center(Org(min_e)->point, Dest(min_e)->point, C, 1.0);
			std::cout << off << '\n';
			break;
		}
	}*/


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

	Graphics g(1000, 500, 100, 200);

	for (auto face = dt.faces.begin(); face != dt.faces.end(); ++face)
	{
		if (face == dt.outer_face) continue;
		auto e = face->bounds;
		if (!face->mark)
			; // g.add_polygon(face, "green"s, 0.25);
		else
			g.add_polygon(face, "white"s, 0.25);
	}
	for (auto qref = dt.edges.begin(); qref != dt.edges.end(); ++qref)
	{
		EdgeRef e(qref);
		if (e.data().var.which() == 1) {
			e = e.Rot();
		}
		if (e.data().fixed)
			g.addLine(Org(e)->point, Dest(e)->point, "black"s, 2.0);
		if (Left(e)->mark || Right(e)->mark)
			g.addLine(Org(e)->point, Dest(e)->point, "black"s, 1.0);
	}
	for (Subdivision::Vertex const& v : dt.vertices)
	{
		double rad = 1.5;
		if (!v.circumcenter)
			rad = 3.0;
		g.addPoint(v.point, rad);
	}


	int trs = std::count_if(dt.faces.begin(), dt.faces.end(), [](Subdivision::Face const& f) {
		return f.mark == 1;
	});
	std::cout << "triangles:\t" << trs << '\n';
	std::ofstream xml{"alper.xml"};
	g.output(xml);

	std::cout << "Euler invariant: " << dt.vertices.size() - dt.edges.size() + dt.faces.size() << '\n';
}