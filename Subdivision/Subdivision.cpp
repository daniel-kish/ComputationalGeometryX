#include "Subdivision.h"

Subdivision::VertexRef& Org(EdgeRef e) {
	Subdivision::EdgeData& data = e.data();
	return boost::get<std::list<Subdivision::Vertex>::iterator>(data.var);
}

Subdivision::VertexRef& Dest(EdgeRef e) {
	Subdivision::EdgeData& data = e.Sym().data();
	return boost::get<std::list<Subdivision::Vertex>::iterator>(data.var);
}

Subdivision::FaceRef& Left(EdgeRef e) {
	Subdivision::EdgeData& data = e.InvRot().data();
	return boost::get<std::list<Subdivision::Face>::iterator>(data.var);
}

Subdivision::FaceRef& Right(EdgeRef e) {
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
	
	e.data() = EdgeData{{},vert};
	vert->leaves = e;
	++vert;
	e.Sym().data() = EdgeData{{},vert};
	vert->leaves = e.Sym();

	//for (VertexRef v = this->vertices.begin(); v != this->vertices.end(); ++v)
	//{
	//	if (!(Org(v->leaves) == v)) {
	//		std::cerr << "Houston, we've another problem!\n";
	//		std::exit(1);
	//	}
	//}
}


void Subdivision::merge(Subdivision &other)
{
	edges.merge(other.edges);
	vertices.splice(vertices.end(), other.vertices);
	faces.splice(faces.end(), other.faces);

	//for (VertexRef v = this->vertices.begin(); v != this->vertices.end(); ++v)
	//{
	//	if (!(Org(v->leaves) == v)) {
	//		std::cerr << "Houston, we've another problem!\n";
	//		std::exit(1);
	//	}
	//}
}

void Subdivision::deleteEdge(EdgeRef e)
{
	VertexRef org = Org(e);
	if (org->leaves == e) {
		if (e == e.Onext()) {
			std::cerr << "Houston, we've a problem!\n";
			std::exit(1);
		}
		org->leaves = e.Onext();
	}

	VertexRef dest = Dest(e);
	if (dest->leaves == e.Sym()) {
		if (e.Sym() == e.Sym().Onext()) {
			std::cerr << "Houston, we've a problem!\n";
			std::exit(1);
		}
		dest->leaves = e.Sym().Onext();
	}

	edges.deleteEdge(e);

	//for (VertexRef v = this->vertices.begin(); v != this->vertices.end(); ++v)
	//{
	//	if (!(Org(v->leaves) == v)) {
	//		std::cerr << "Houston, we've another problem!\n";
	//		std::exit(1);
	//	}
	//}
}

EdgeRef Subdivision::add_vertex(EdgeRef e, Point p)
{
	vertices.push_back(Vertex{p});
	auto vert = std::prev(vertices.end());

	auto a = edges.makeEdge();
	a.data() = EdgeData{{},Dest(e)};
	a.Sym().data() = EdgeData{{},vert};
	vert->leaves = a.Sym();

	splice(e.Lnext(), a);

	//for (VertexRef v = this->vertices.begin(); v != this->vertices.end(); ++v)
	//{
	//	if (!(Org(v->leaves) == v)) {
	//		std::cerr << "Houston, we've another problem!\n";
	//		std::exit(1);
	//	}
	//}

	return a;
}

EdgeRef Subdivision::connect(EdgeRef a, EdgeRef b)
{
	auto e = edges.makeEdge();
	e.data() = EdgeData{{}, Dest(a)};
	e.Sym().data() = EdgeData{{},Org(b)};

	splice(e, a.Lnext());
	splice(e.Sym(), b);

	//for (VertexRef v = this->vertices.begin(); v != this->vertices.end(); ++v)
	//{
	//	if (!(Org(v->leaves) == v)) {
	//		std::cerr << "Houston, we've another problem!\n";
	//		std::exit(1);
	//	}
	//}
	return e;
}
