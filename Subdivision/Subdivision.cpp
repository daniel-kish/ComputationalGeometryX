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

void printEdge(EdgeRef e)
{
	std::cout << Org(e)->point << ' ' << Dest(e)->point << '\n';
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
	
	e.data() = EdgeData{{},{},vert};
	vert->leaves = e;
	++vert;
	e.Sym().data() = EdgeData{{},{},vert};
	vert->leaves = e.Sym();

	//for (VertexRef v = this->vertices.begin(); v != this->vertices.end(); ++v)
	//{
	//	if (!(Org(v->leaves) == v)) {
	//		std::cerr << "ctor: v->leaves\n";
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
	//		std::cerr << "merge: v->leaves\n";
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
	//		std::cerr << "deleteEdge: v->leaves\n";
	//		std::exit(1);
	//	}
	//}
}

EdgeRef Subdivision::add_vertex(EdgeRef e, Point p)
{
	VertexRef eo = Org(e);
	VertexRef ed = Dest(e);

	vertices.push_back(Vertex{p});
	auto vert = std::prev(vertices.end());

	auto a = edges.makeEdge();
	a.data() = EdgeData{{},{},Dest(e)};
	a.Sym().data() = EdgeData{{},{},vert};
	vert->leaves = a.Sym();

	splice(e.Lnext(), a);

	//for (VertexRef v = this->vertices.begin(); v != this->vertices.end(); ++v)
	//{
	//	if (!(Org(v->leaves) == v)) {
	//		std::cerr << "add_vertex: v->leaves\n";
	//		std::exit(1);
	//	}
	//}
	//std::cerr << "add_vertex: ok\n";

	return a;
}

EdgeRef Subdivision::connect(EdgeRef a, EdgeRef b)
{
	auto e = edges.makeEdge();
	e.data() = EdgeData{{},{}, Dest(a)};
	e.Sym().data() = EdgeData{{},{},Org(b)};

	splice(e, a.Lnext());
	splice(e.Sym(), b);

	//for (VertexRef v = this->vertices.begin(); v != this->vertices.end(); ++v)
	//{
	//	if (!(Org(v->leaves) == v)) {
	//		std::cerr << "connect: v->leaves\n";
	//		std::exit(1);
	//	}
	//}
	return e;
}

// with faces initialized 

EdgeRef Subdivision::splitVertex(EdgeRef a, EdgeRef b, Point const& p)
{
	vertices.push_back(Vertex{p});
	VertexRef vnew = std::prev(vertices.end());

	EdgeRef e = edges.makeEdge().Rot();

	splice(a, e);
	splice(b, e.Sym());

	e.data().var = Org(b);
	e.InvRot().data().var = Left(a);
	e.Rot().data().var = Left(b);
	Org(b)->leaves = e;

	auto eiter = e.Sym();
	auto end{eiter};
	do {
		eiter.data().var = vnew;
		eiter = eiter.Onext();
	} while (eiter != end);
	vnew->leaves = e.Sym();

	return e;
}

void Subdivision::joinVertex(EdgeRef e)
{
	VertexRef todelete = Dest(e);
	auto a = e.Lnext();
	auto b = e.Oprev();

	VertexRef org_b = Org(b);

	if (a == e.Sym())
		a = b;

	splice(b, e.Sym());
	splice(a, e);

	auto eiter = a;
	auto end = eiter;
	do {
		eiter.data().var = org_b;
		eiter = eiter.Onext();
	} while (eiter != end);

	Org(b)->leaves = b;
	Left(b)->bounds = b;
	Left(a)->bounds = a;

	vertices.erase(todelete);
	edges.deleteEdge(e);
}

EdgeRef Subdivision::splitFace(EdgeRef a, EdgeRef b)
{
	assert(Left(a) == Left(b));

	FaceRef face = Left(a);
	VertexRef org = Org(a), dest = Org(b);

	faces.push_back(Face{});
	FaceRef fnew = std::prev(faces.end());

	EdgeRef e = edges.makeEdge();

	splice(e, a);
	splice(e.Sym(), b);

	e.data().var = org;
	e.Sym().data().var = dest;
	e.InvRot().data().var = face;
	face->bounds = e;

	auto eiter = e.Rot();
	auto end = eiter;

	do {
		eiter.data().var = fnew;
		eiter = eiter.Onext();
	} while (eiter != end);
	fnew->bounds = e.Sym();

	return e;
}

void Subdivision::joinFace(EdgeRef e)
{
	EdgeRef a = e.Oprev(), b = e.Lnext();

	if (a == e.Sym())
		a = b;

	splice(b, e.Sym());
	splice(a, e);
	FaceRef left_b = Left(b);

	auto eiter = a.InvRot();
	auto end = eiter;
	do {
		eiter.data().var = left_b;
		eiter = eiter.Onext();
	} while (eiter != end);

	Org(a)->leaves = a;
	Org(b)->leaves = b;
	Left(b)->bounds = b;

	faces.erase(Right(e));
	edges.deleteEdge(e);
}