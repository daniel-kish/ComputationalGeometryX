#define _USE_MATH_DEFINES
#include <iostream>
#include <utility>
#include <list>
#include <string>
#include "Point.h"
#include <vector>
#include <array>
#include <random>
#include <iterator>
#include <cmath>
#include <tuple>
#include <sstream>
#include <fstream>
#include <numeric>
#include <iomanip>
#include "boost/rational.hpp"
#include "boost/multiprecision/cpp_int.hpp"


struct Subdivision {
	struct QuadEdge;
	using QuadEdgeList = std::list<QuadEdge>;
	using QuadEdgeRef = QuadEdgeList::iterator;

	struct EdgeRef {
		QuadEdgeRef ptr;
		int n;

		Point** data() {
			return &(ptr->edges[n].second);
		}

		Point& Org() {
			return *ptr->edges[n].second;
		}
		Point& Dest() {
			return *ptr->edges[(n + 2) % 4].second;
		}

		EdgeRef Onext() {
			return ptr->edges[n].first;
		}
		EdgeRef Rot() {
			EdgeRef copy = *this;
			copy.n = (copy.n + 1) % 4;
			return copy;
		}
		EdgeRef InvRot() {
			EdgeRef copy = *this;
			copy.n = (copy.n + 3) % 4;
			return copy;
		}
		EdgeRef Sym() {
			return (*this).Rot().Rot();
		}
		EdgeRef Lnext() {
			return this->InvRot().Onext().Rot();
		}
		EdgeRef Oprev() {
			return this->Rot().Onext().Rot();
		}
		EdgeRef Rprev() {
			return this->Sym().Onext();
		}
		bool operator==(EdgeRef const& other) {
			return this->ptr == other.ptr && this->n == other.n;
		}
		bool operator!=(EdgeRef const& other) {
			return !(*this == other);
		}
		
	}; 
	
	struct QuadEdge { std::pair<EdgeRef, Point*> edges[4]; };

	QuadEdgeList quadEdges;
	EdgeRef makeEdge() {
		quadEdges.push_back(QuadEdge{});
		QuadEdgeRef qref = std::prev(quadEdges.end());
		qref->edges[0] = {{qref,0},nullptr};
		qref->edges[1] = {{qref,3},nullptr};
		qref->edges[2] = {{qref,2},nullptr};
		qref->edges[3] = {{qref,1},nullptr};
		return EdgeRef{qref,0};
	}
};

void splice(Subdivision::EdgeRef a, Subdivision::EdgeRef b)
{
	auto alpha = a.Onext().Rot();
	auto beta = b.Onext().Rot();

	std::swap(a.ptr->edges[a.n].first, b.ptr->edges[b.n].first);
	std::swap(alpha.ptr->edges[alpha.n].first, beta.ptr->edges[beta.n].first);
}

// Delaunay operators
Subdivision::EdgeRef connect(Subdivision& s, Subdivision::EdgeRef a, Subdivision::EdgeRef b)
{
	auto e = s.makeEdge();
	*e.data() = &a.Dest();
	*e.Sym().data() = &b.Org();

	splice(e, a.Lnext());
	splice(e.Sym(), b);
	return e;
}

void deleteEdge(Subdivision& s, Subdivision::EdgeRef e)
{
	splice(e, e.Oprev());
	splice(e.Sym(), e.Sym().Oprev());
	s.quadEdges.erase(e.ptr); // delete e
}

void swap(Subdivision::EdgeRef e)
{
	auto a = e.Oprev();
	auto b = e.Sym().Oprev();
	splice(e, a); splice(e.Sym(), b);
	splice(e, a.Lnext()); splice(e.Sym(), b.Lnext());
	*e.data() = *a.Sym().data();
	*e.Sym().data() = *b.Sym().data();
}

using cpp_int = boost::multiprecision::cpp_int;
using int_type = cpp_int;
using rational = boost::rational<int_type>;


rational from_double(double x)
{
	double significand{};
	int exp{};
	significand = frexp(x, &exp);
	rational r(int_type(significand*pow(2.0, 53)), int_type(pow(2.0, 53)));
	rational expr;
	if (exp >= 0)
		expr = rational(int_type(pow(2.0, exp)));
	else
		expr = rational(1, int_type(pow(2.0, -exp)));
	rational res = r*expr;

	//double err = double(res.numerator()) / double(res.denominator()) - x;
	//if (err != 0.0) {
	//	std::cout << "error: " << err << '\n';
	//	//std::exit(1);
	//}
	return res;
}

template <class T>
T det_3(std::array<T,9> const& m)
{
	return  + m[0] * (m[4] * m[8] - m[7] * m[5])
			- m[1] * (m[3] * m[8] - m[6] * m[5])
			+ m[2] * (m[3] * m[7] - m[6] * m[4]);
}

template <class T>
T det_4(std::array<T, 16> const& m)
{
	T d1 = m[10] * m[15] - m[14] * m[11];
	T d2 = m[9] * m[15] - m[13] * m[11];
	T d3 = m[9] * m[14] - m[13] * m[10];
	T d4 = m[8] * m[15] - m[12] * m[11];
	T d5 = m[8] * m[14] - m[12] * m[10];
	T d6 = m[8] * m[13] - m[12] * m[9];

	T M11 = m[5] * d1 - m[6] * d2 + m[7] * d3;
	
	T M12 = m[4] * d1 - m[6] * d4 + m[7] * d5;

	T M13 = m[4] * d2 - m[5] * d4 + m[7] * d6;

	T M14 = m[4] * d3 - m[5] * d5 + m[6] * d6;

	return m[0] * M11 - m[1] * M12 + m[2] * M13 - m[3] * M14;
}

bool inCircle(Point const& a, Point const& b, Point const& c, Point const& d)
{
	auto x = [](const Point& p) { return p.x; };
	auto y = [](const Point& p) { return p.y; };

	std::array<double, 16> md{
		x(a), y(a), sqNorm(a), 1.0,
		x(b), y(b), sqNorm(b), 1.0,
		x(c), y(c), sqNorm(c), 1.0,
		x(d), y(d), sqNorm(d), 1.0
	};
	std::array<rational, 16> m;
	std::transform(md.begin(), md.end(), m.begin(), [](double x) -> rational {
		return from_double(x);
	});
	rational det = det_4(m);
	//std::cout << a << ' ' << b << ' ' << c << ' ' << d << ' ' << det << '\n';

	return det > rational(0);
}

bool ccw(Point const& a, Point const& b, Point const& c)
{
	auto x = [](const Point& p) { return p.x; };
	auto y = [](const Point& p) { return p.y; };

	std::array<double, 9> md{
		x(a), y(a), 1.0,
		x(b), y(b), 1.0,
		x(c), y(c), 1.0
	};
	std::array<rational, 9> m;
	std::transform(md.begin(), md.end(), m.begin(), [](double x) -> rational {
		return from_double(x); 
	});

	rational det = det_3(m);
	//std::cout /*<< a << ' ' << b << ' ' << c << ' '*/ << det << '\n';
	return det > rational(0);
}

bool rightOf(Point const& p, Subdivision::EdgeRef& e)
{
	return ccw(p, e.Dest(), e.Org());
}

bool leftOf(Point const& p, Subdivision::EdgeRef& e)
{
	return ccw(p, e.Org(), e.Dest());
}

std::tuple<Subdivision, Subdivision::EdgeRef, Subdivision::EdgeRef>
delaunay(std::vector<Point>::iterator b, std::vector<Point>::iterator e)
{
	if (std::distance(b, e) == 2) {
		Subdivision s;
		auto a = s.makeEdge();
		*a.data() = &*b;
		*a.Sym().data() = &*(b + 1);
		return std::make_tuple(std::move(s),a,a.Sym());
	}
	else if (std::distance(b, e) == 3) {
		Subdivision s;
		auto s1 = b;
		auto s2 = b+1;
		auto s3 = b+2;

		auto a = s.makeEdge(); auto b = s.makeEdge();
		splice(a.Sym(), b);
		*a.data() = &*s1; 
		*a.Sym().data() = *b.data() = &*s2; 
		*b.Sym().data() = &*s3;
		
		if (ccw(*s1, *s2, *s3)) {
			connect(s, b, a);
			return std::make_tuple(std::move(s), a, b.Sym());
		}
		else if (ccw(*s3, *s2, *s1)) {
			auto c = connect(s, b, a);
			return std::make_tuple(std::move(s), c.Sym(), c);
		}
		else 
			return std::make_tuple(std::move(s), a, b.Sym());
	}
	// else if S.size() >= 4
	auto mid{b}; std::advance(mid, std::distance(b, e) / 2);
	Subdivision L, R;
	Subdivision::EdgeRef ldo, ldi, rdi, rdo;
	std::tie(L, ldo, ldi) = delaunay(b, mid);
	std::tie(R, rdi, rdo) = delaunay(mid, e);

	L.quadEdges.splice(L.quadEdges.end(), R.quadEdges); // physical merge first
	// then semantic merge:
	do {
		if (leftOf(rdi.Org(), ldi)) ldi = ldi.Lnext();
		else if (rightOf(ldi.Org(), rdi)) rdi = rdi.Rprev();
		else break;
	} while (true);

	auto basel = connect(L, rdi.Sym(), ldi);
	if (*ldi.data() == *ldo.data()) ldo = basel.Sym();
	if (*rdi.data() == *rdo.data()) rdo = basel;

	auto valid = [&basel](Subdivision::EdgeRef e) {return rightOf(e.Dest(), basel); };
	do // merge loop
	{
		auto lcand = basel.Sym().Onext();
		if (valid(lcand)) {
			while ( inCircle(basel.Dest(), basel.Org(), lcand.Dest(), lcand.Onext().Dest()) )
			{
				lcand = lcand.Onext();
				deleteEdge(L, lcand.Oprev());
			}
		}
		auto rcand = basel.Oprev();
		if (valid(rcand)) {
			while( inCircle(basel.Dest(), basel.Org(),rcand.Dest(),rcand.Oprev().Dest()) )
			{
				rcand = rcand.Oprev();
				deleteEdge(L, rcand.Onext());
			}
		}
		if (!valid(lcand) && !valid(rcand)) break; // we've reached the top

		if (!valid(lcand) ||
			valid(rcand) && inCircle(lcand.Dest(), lcand.Org(), rcand.Org(), rcand.Dest()))
			basel = connect(L, rcand, basel.Sym());
		else
			basel = connect(L, basel.Sym(), lcand.Sym());
	} while (true);
	return std::make_tuple(std::move(L), ldo, rdo);
}

struct Graphics
{
	double scale{30};
	double wid, height;
	double X0, Y0;
	std::string main_code;
	std::string ending;
	double pointRad{3.0};
	double lineWidth{2.0};
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


int main()
{
	Rect rect{{-5,-5},{13,7}};
	std::vector<Point> S = rectHull(rect, 13*2, 7*2);
	
	std::sort(S.begin(), S.end());
	auto le = std::unique(S.begin(), S.end(), [](Point p, Point q) {
		return std::abs(p.x - q.x) < 1.0e-10 && std::abs(p.y - q.y) < 1.0e-10;
	});
	S.erase(le, S.end());


	Subdivision s;
	Subdivision::EdgeRef l, r;
	
	std::tie(s,l,r) = delaunay(S.begin(), S.end());
	
	std::cout << S.size() << ' ' << s.quadEdges.size() << '\n';
	
	Graphics ga;
	for (Point const& p : S)
		ga.addPoint(p);
	for (Subdivision::QuadEdge& q : s.quadEdges)
		ga.addLine(*q.edges[0].second, *q.edges[2].second);
	std::ofstream outstr{"trian.xml"};
	ga.output(outstr);
}
