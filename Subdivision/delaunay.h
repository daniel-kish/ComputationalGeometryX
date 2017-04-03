#pragma once
#include "Subdivision.h"
#include <tuple>
#include <vector>

void swap(Edge e);
bool rightOf(Vertex v, Edge e);
bool rightOf(Point p, Edge e);
bool leftOf(Vertex v, Edge e);
bool leftOf(Point p, Edge e);
bool incircle(Vertex a, Vertex b, Vertex c, Vertex d);

std::tuple<Subdivision, Edge, Edge>
delaunay_dnc(std::vector<Point>::iterator b, std::vector<Point>::iterator e);

std::tuple<Edge,bool> locate(Subdivision& s, Point x);
std::tuple<Edge, bool> locate(Subdivision & s, Point x, Edge e);
std::tuple<Edge, bool> insertSite(Subdivision& s, Point x);
std::tuple<Edge, bool> insertSite(Subdivision& s, Point x, Edge start);
void insertSiteSequence(Subdivision& s, std::vector<Point> seq);