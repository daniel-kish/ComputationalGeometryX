#define _USE_MATH_DEFINES
#include <iostream>
#include <utility>
#include <list>
#include <vector>
#include "QuadEdge.h"
#include "boost/variant.hpp"

struct Point { double x, y; };
std::ostream& operator<<(std::ostream& os, Point p)
{
	return os << '(' << p.x << ' ' << p.y << ')';
}

int main()
{
	
	QuadEdgeList<boost::variant<Point*, std::string*> > list;
	std::vector<Point> pts{{0,0},{1,0},{0,1}};
	std::vector<std::string> faces{{"happy triangle"}};


	auto e1 = list.makeEdge();
	auto e2 = list.makeEdge();
	auto e3 = list.makeEdge();
	splice(e1,e2);
	splice(e1.Sym(), e3);
	splice(e3.Sym(), e2.Sym());

	
	e1.data() = e2.data() = &pts[0];
	e1.Sym().data() = e3.data() = &pts[1];
	e2.Sym().data() = e3.Sym().data() = &pts[2];
	
	e1.InvRot().data() = &faces[0];
	e1.Lnext().InvRot().data() = &faces[0];
	e1.Lnext().Lnext().InvRot().data() = &faces[0];
	
	auto ei = e1;
	do {
		std::cout << *boost::get<std::string*>(ei.InvRot().data()) << '\n';
		ei = ei.Lnext();
	} while (ei != e1);


	for (QuadEdgeList<boost::variant<Point*, std::string*> >::QuadEdge const& q : list) {
		std::cout << *boost::get<Point*>(q.edges[0].second)
			<< ' ' << *boost::get<Point*>(q.edges[2].second) << '\n';
	}
}
