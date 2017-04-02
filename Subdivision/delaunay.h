#pragma once
#include "Subdivision.h"
#include <tuple>
#include <vector>

void swap(Edge e);
bool rightOf(Vertex v, Edge e);
bool leftOf(Vertex v, Edge e);
double incircle(Vertex a, Vertex b, Vertex c, Vertex d);

std::tuple<Subdivision, Edge, Edge>
delaunay_dnc(std::vector<Point>::iterator b, std::vector<Point>::iterator e);
