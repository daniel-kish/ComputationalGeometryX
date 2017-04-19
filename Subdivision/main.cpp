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
	double lineWidth{1.0};
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

	void addLine(Point const& p, Point const& q, double width = 0.0)
	{
		if (width == 0.0)
			width = lineWidth;
		std::ostringstream str;
		str << R"( <line x1=")" << p.x << R"(" y1=")" << p.y
			<< R"(" x2=")" << q.x << R"(" y2=")" << q.y
			<< R"(" style="stroke:black;stroke-width:)" << width / scale << R"("/> )"
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
				<< R << ',' << G << ',' << B << "); opacity:" << opacity << "; stroke-width:0\" />";
		}
		else {
			str << "<polygon points = \"" <<
				a.x << ',' << a.y << ' ' <<
				b.x << ',' << b.y << ' ' <<
				c.x << ',' << c.y << "\" style=\"fill:" << color
				<< "; opacity:" << opacity << "; stroke-width:0\" />";
		}
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
	using namespace std::literals;
	exactinit();

	std::vector<Point> model{
		{-4,2}, {-4,1}, {-3,1}, {-3,-1}, {-4,-1}, {-4,-2}, {3,-2}, {4,-1},{4,2}
	};
	auto trian = triangleCover(model);
	std::vector<Point> cover(trian.begin(), trian.end());


	std::sort(cover.begin(), cover.end());

	Subdivision dt;
	EdgeRef l, r;

	std::tie(dt, l, r) = delaunay_dnc(cover.begin(), cover.end());


	std::vector<VertexRef> inserted;
	inserted.reserve(model.size());
	for (Point const& p : model)
		inserted.push_back(insertSite(dt, p));

	for (unsigned i = 0; i < inserted.size(); ++i) {
		auto ins = insertEdge(dt, inserted[i], inserted[(i + 1) % inserted.size()]);
		ins.data().boundary = true;
		ins.Sym().data().boundary = true;
	}

	std::vector<Point> hole{
		{-1,-1},{1,-1},{-0.8,0},{1,1},{-1,1}
	};
	inserted.clear();
	inserted.reserve(hole.size());
	for (Point const& p : hole)
		inserted.push_back(insertSite(dt, p));
	for (unsigned i = 0; i < inserted.size(); ++i) {
		auto ins = insertEdge(dt, inserted[i], inserted[(i + 1) % inserted.size()]);
		ins.data().boundary = true;
		ins.Sym().data().boundary = true;
	}

	VertexRef a = insertSite(dt, {-3.5,1.5});
	VertexRef b = insertSite(dt, {3.5,1.5});
	insertEdge(dt, a, b);

	init_faces(dt);
	mark_outer_faces(dt, r);

	std::cout << dt.edges.size() << '\n';
	int n = 60;
	while (n--)
	{
		for (auto qref = dt.edges.begin(); qref != dt.edges.end(); ++qref)
		{
			EdgeRef e(qref);
			if (e.data().var.which() != 0)
				e = e.Rot();

			if (!e.data().boundary) continue;
			if (encroaches(e, Dest(e.Onext())) || encroaches(e, Dest(e.Oprev()))) {
				splitBoundaryEdge(dt, e);
				break;
			}
		}
	}
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
		/*std::cout << "start\n";
		auto eiter = f->bounds;
		auto end = eiter;
		do {
			printEdge(eiter);
			eiter = eiter.Lnext();
		} while (eiter != end);
		std::cout << "end\n";*/

		if (Left(f->bounds) != f) {
			std::cerr << "faces connectivity check failed, boundary of 'f' :\n";
			auto eiter = f->bounds;
			auto end = eiter;
			do {
				printEdge(eiter);
				eiter = eiter.Lnext();
			} while (eiter != end);

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

	Graphics g;

	for (auto face = dt.faces.begin(); face != dt.faces.end(); ++face)
	{
		if (face == dt.outer_face) continue;
		auto e = face->bounds;
		if (!face->mark) {
			g.add_polygon(Org(e)->point, Dest(e)->point, Dest(e.Onext())->point, "cornsilk");
		}
		else
			g.add_polygon(Org(e)->point, Dest(e)->point, Dest(e.Onext())->point, {}, 0.5);
	}
	for (auto qref = dt.edges.begin(); qref != dt.edges.end(); ++qref)
	{
		EdgeRef e(qref);
		if (e.data().var.which() == 1) {
			e = e.Rot();
		}
		if (e.data().fixed) {
			g.addLine(Org(e)->point, Dest(e)->point, 2.0);
		}
		else  // if (Left(e)->mark || Right(e)->mark)
			g.addLine(Org(e)->point, Dest(e)->point);
	}
	for (Subdivision::Vertex const& p : dt.vertices)
		g.addPoint(p.point);

	std::ofstream xml{"trian.xml"};
	g.output(xml);

	std::cout << "Euler invariant: " << dt.vertices.size() - dt.edges.size() + dt.faces.size() << '\n';
}