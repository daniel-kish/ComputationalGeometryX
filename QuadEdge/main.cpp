#define _USE_MATH_DEFINES
#include <iostream>
#include "Header.h"

//
//struct Point { double x, y; };
//
//std::ostream& operator<<(std::ostream& os, Point p)
//{
//	return os << '(' << p.x << ' ' << p.y << ')';
//}
//
//using Data = boost::variant<Point*, std::string*>;
//

int main()
{
	QuadEdgeList<int> list;

	QuadEdgeList<int>::EdgeRef e;
	e = list.makeEdge();
	list.makeEdge();


	std::cout << list.size() << '\n';

	for (QuadEdgeList<int>::QuadEdgeRef i = list.begin(); i != list.end(); ++i)
	{
		QuadEdgeList<int>::EdgeRef e(i, 0);
		std::cout << e.data() << '\n';
	}
	for (QuadEdgeList<int>::QuadEdge const& q : list)
	{
		std::cout << q.records()[0].data << '\n';
	}
}
