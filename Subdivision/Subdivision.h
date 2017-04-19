#pragma once
#include "Point.h"
#include "boost/variant.hpp"
#include <list>
#include <vector>
#include "QuadEdge.h"

struct Subdivision 
{
	struct Vertex;
	struct Face;
	struct EdgeData {
		bool fixed; 
		bool boundary;
		boost::variant<std::list<Vertex>::iterator, std::list<Face>::iterator> var;
	};

	using Edges = QuadEdgeList<EdgeData>;
	struct Vertex {
		Point point;
		Edges::EdgeRef leaves;
	};
	struct Face {
		int mark{-1};
		Edges::EdgeRef bounds;
	};
	using VertexRef = std::list<Subdivision::Vertex>::iterator;
	using FaceRef = std::list<Subdivision::Face>::iterator;


	std::list<Face> faces;
	std::list<Vertex> vertices;
	Edges edges;
	FaceRef outer_face;

	Subdivision();
	Subdivision(Point p1, Point p2);
	Subdivision::Edges::EdgeRef connect(Subdivision::Edges::EdgeRef a, Subdivision::Edges::EdgeRef b);
	Subdivision::Edges::EdgeRef add_vertex(Subdivision::Edges::EdgeRef, Point);
	void deleteEdge(Subdivision::Edges::EdgeRef);
	void merge(Subdivision&);
	Subdivision::Edges::EdgeRef 
		splitVertex(Subdivision::Edges::EdgeRef a, Subdivision::Edges::EdgeRef b, Point const&);
	void joinVertex(Subdivision::Edges::EdgeRef);
	Subdivision::Edges::EdgeRef
		splitFace(Subdivision::Edges::EdgeRef a, Subdivision::Edges::EdgeRef b);
	void joinFace(Subdivision::Edges::EdgeRef a);
};

using EdgeRef = Subdivision::Edges::EdgeRef;
using VertexRef = Subdivision::VertexRef;
using FaceRef = Subdivision::FaceRef;

Subdivision::VertexRef& Org(EdgeRef e);
Subdivision::VertexRef& Dest(EdgeRef e);
Subdivision::FaceRef& Left(EdgeRef e);
Subdivision::FaceRef& Right(EdgeRef e);

void printEdge(EdgeRef);