#define _USE_MATH_DEFINES
#include <iostream>
#include "QuadEdge.h"
#include "boost\variant.hpp"


int main()
{
	QuadEdgeList<int> list;

	list.makeEdge();
	auto e = list.makeEdge();

	list.deleteEdge(e);

	std::cout << list.size() << '\n';

	for (QuadEdgeList<int>::QuadEdgeRef i = list.begin(); i != list.end(); ++i)
	{
		QuadEdgeList<int>::EdgeRef e(i, 0);
		std::cout << e.data() << ' ' << e.Sym().data() << '\n';
	}
	for (QuadEdgeList<int>::QuadEdge const& q : list)
	{
		std::cout << q.records()[0].data << '\n';
	}
}
