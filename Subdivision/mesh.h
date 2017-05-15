#pragma once
#include "Subdivision.h"

void init_faces(Subdivision&);
void mark_outer_faces(Subdivision& s, EdgeRef outEdge);
bool encroaches(EdgeRef e, VertexRef v);
bool encroaches(EdgeRef e, Point p);
Point midpoint(EdgeRef e);
EdgeRef splitBoundaryEdge(Subdivision& s, EdgeRef e);
EdgeRef splitRegularEdge(Subdivision& s, EdgeRef e);

void splitEdges(Subdivision& s);

bool inside(Subdivision& s, Point x);
VertexRef insertMeshSite(Subdivision& s, Point x);
void eliminate_worst_triangle(Subdivision& dt);
void eliminate_triangle(Subdivision& dt, FaceRef face);
std::tuple<FaceRef, double> find_worst(Subdivision & dt);