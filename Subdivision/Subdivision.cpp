#include "Subdivision.h"

Subdivision::VertexRef& Org(Edge e) {
	Subdivision::EdgeData& data = e.data();
	return boost::get<std::list<Subdivision::Vertex>::iterator>(data.var);
}

Subdivision::VertexRef& Dest(Edge e) {
	Subdivision::EdgeData& data = e.Sym().data();
	return boost::get<std::list<Subdivision::Vertex>::iterator>(data.var);
}

Subdivision::FaceRef& Left(Edge e) {
	Subdivision::EdgeData& data = e.InvRot().data();
	return boost::get<std::list<Subdivision::Face>::iterator>(data.var);
}

Subdivision::FaceRef& Right(Edge e) {
	Subdivision::EdgeData& data = e.Rot().data();
	return boost::get<std::list<Subdivision::Face>::iterator>(data.var);
}

Subdivision::Subdivision()
{
}

Subdivision::Subdivision(Point p1, Point p2)
{
	vertices.push_back(Vertex{p1});
	vertices.push_back(Vertex{p2});
	
	auto e = edges.makeEdge();

	auto vert = vertices.begin();
	e.data() = EdgeData{{},vert++};
	e.Sym().data() = EdgeData{{},vert};
}


void Subdivision::merge(Subdivision &other)
{
	edges.merge(other.edges);
	vertices.splice(vertices.end(), other.vertices);
	faces.splice(faces.end(), other.faces);
}

void Subdivision::deleteEdge(Edge e)
{
	edges.deleteEdge(e);
}

Edge Subdivision::add_vertex(Edge e, Point p)
{
	vertices.push_back(Vertex{p});
	auto vert = std::prev(vertices.end());

	auto a = edges.makeEdge();
	a.data() = EdgeData{{},Dest(e)};
	a.Sym().data() = EdgeData{{},vert};

	splice(e.Lnext(), a);

	return a;
}

Edge Subdivision::connect(Edge a, Edge b)
{
	auto e = edges.makeEdge();
	e.data() = EdgeData{{}, Dest(a)};
	e.Sym().data() = EdgeData{{},Org(b)};

	splice(e, a.Lnext());
	splice(e.Sym(), b);

	return e;
}
