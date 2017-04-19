#pragma once
#include "Subdivision.h"
#include <tuple>
#include <vector>
#include <utility>

void swap(EdgeRef e);
bool rightOf(VertexRef v, EdgeRef e);
bool leftOf(VertexRef v, EdgeRef e);

bool leftOf(VertexRef v, std::pair<VertexRef, VertexRef>);
bool rightOf(VertexRef v, std::pair<VertexRef, VertexRef>);

bool leftOf(Point p, EdgeRef e);
bool rightOf(Point p, EdgeRef e);

bool incircle(VertexRef a, VertexRef b, VertexRef c, VertexRef d);

std::tuple<Subdivision, EdgeRef, EdgeRef>
delaunay_dnc(std::vector<Point>::iterator b, std::vector<Point>::iterator e);

EdgeRef locate(Subdivision& s, Point x);
EdgeRef locate(Subdivision & s, Point x, EdgeRef e);
VertexRef insertSite(Subdivision& s, Point x);
VertexRef insertSite(Subdivision& s, Point x, EdgeRef start);
void insertSiteSequence(Subdivision& s, std::vector<Point> seq);

EdgeRef insertEdge(Subdivision&, VertexRef, VertexRef);