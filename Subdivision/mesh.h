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
void eliminate_triangle(Subdivision& dt, FaceRef face);

bool eliminate_worst_triangle(Subdivision& dt, double min_ratio);
bool eliminate_worst_triangle(Subdivision& dt, double min_ratio, double min_area);

bool eliminate_bad_triangle(Subdivision& dt, double min_ratio);
bool eliminate_bad_triangle(Subdivision& dt, double min_ratio, double min_area);


std::tuple<FaceRef, double> find_worst(Subdivision & dt);
std::tuple<FaceRef, double> find_biggest(Subdivision & dt);
std::tuple<FaceRef, double> find_smallest(Subdivision & dt);

FaceRef find_bad(Subdivision & dt, double ratio, double area);
FaceRef find_bad(Subdivision & dt, double ratio);

void ruppert_refinement(Subdivision & dt, double min_ratio, int max_iters=1000);
void ruppert_refinement(Subdivision & dt, double min_ratio, double min_area, int max_iters=1000);


void insertClosedLoop(Subdivision & dt, std::vector<Point> const& hole);

void deleteSite_wf(Subdivision& dt, VertexRef v);

bool chew_2nd_eliminate_worst(Subdivision& dt, double min_ratio);

void chew_2nd_refinement(Subdivision& dt, double min_ratio, int max_iters=1000);
void chew_2nd_refinement(Subdivision& dt, double min_ratio, double min_area, int max_iters = 1000);

void off_center_correction(Subdivision& dt, VertexRef site, double min_angle);