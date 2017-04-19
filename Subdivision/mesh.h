#pragma once
#include "Subdivision.h"

void init_faces(Subdivision&);
void mark_outer_faces(Subdivision& s, EdgeRef outEdge);
bool encroaches(EdgeRef e, VertexRef v);
Point midpoint(EdgeRef e);
EdgeRef splitBoundaryEdge(Subdivision& s, EdgeRef e);