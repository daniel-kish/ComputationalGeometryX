#include "mesh.h"
#include <vector>
#include "Subdivision.h"
#include "delaunay.h"
#include "predicates.h"
#include "geom.h"
#include "boost/math/constants/constants.hpp"

void init_faces(Subdivision& s)
{
	for (auto qref = s.edges.begin(); qref != s.edges.end(); ++qref)
	{
		EdgeRef e(qref, 1);
		//if (e.data().var.which() == 0) { // if references VertexRef
		//	e = e.Rot();
		//}

		e.data().var = s.faces.end();
		e.Sym().data().var = s.faces.end();
	}

	EdgeRef start{s.edges.begin(),1};
	//if (start.data().var.which() == 0) { // if references VertexRef
	//	start = start.Rot();
	//}

	std::vector<EdgeRef> stack{start};

	auto face = [](EdgeRef e) {
		return boost::get<Subdivision::FaceRef>(e.data().var);
	};

	while (!stack.empty())
	{
		EdgeRef e = stack.back(); stack.pop_back();
		if (face(e) != s.faces.end())
			continue;
		s.faces.push_back({});
		auto f = std::prev(s.faces.end());
		auto end = e;
		do {
			e.data().var = f;
			if (face(e.Sym()) == s.faces.end())
				stack.push_back(e.Sym());
			e = e.Onext();
		} while (e != end);
		f->bounds = e.Rot();
	}
}

void mark_outer_faces(Subdivision& s, EdgeRef outEdge)
{
	// 0 - out
	// 1 - in
	Left(outEdge)->mark = 0;
	std::vector<EdgeRef> stack{outEdge};
	while (!stack.empty())
	{
		EdgeRef e = stack.back(); stack.pop_back();
		int current_mark{Left(e)->mark};

		auto end = e;
		do {
			if (Right(e)->mark == -1) {
				Right(e)->mark = (current_mark + e.data().boundary) % 2;
				stack.push_back(e.Sym());
			}
			e = e.Lnext();
		} while (e != end);
	}
	s.outer_face = Left(outEdge);
}

bool encroaches(EdgeRef e, VertexRef v)
{
	double halflen = dist(Org(e)->point, Dest(e)->point) / 2.0;
	Point c = (Org(e)->point + Dest(e)->point)*0.5;
	if (dist(c, v->point) <= halflen)
		return true;
	return false;
}

bool encroaches(EdgeRef e, Point p)
{
	double halflen = dist(Org(e)->point, Dest(e)->point) / 2.0;
	Point c = (Org(e)->point + Dest(e)->point)*0.5;
	if (dist(c, p) <= halflen)
		return true;
	return false;
}

Point midpoint(EdgeRef e)
{
	return (Org(e)->point + Dest(e)->point) * 0.5;
}

EdgeRef swap_wf(Subdivision &s, EdgeRef e)
{
	assert(Left(e)->mark == Right(e)->mark);

	auto a = e.Oprev().Lnext();
	auto b = e.Lnext().Lnext();

	auto mark = Left(e)->mark;
	s.joinFace(e);
	e = s.splitFace(a, b);
	Left(e)->mark = mark;
	Right(e)->mark = mark;

	return e;
}

EdgeRef splitBoundaryEdge(Subdivision& s, EdgeRef e)
{
	assert(e.data().boundary);
	assert(e.Sym().data().boundary);
	assert(e.data().fixed);
	assert(e.Sym().data().fixed);
	assert(Left(e)->mark != Right(e)->mark);

	if (Left(e)->mark)
		e = e.Sym();

	auto e1 = s.splitVertex(e, e.Oprev(), midpoint(e));
	e1.data().boundary = true;
	e1.Sym().data().boundary = true;
	e1.data().fixed = true;
	e1.Sym().data().fixed = true;

	e1 = e1.Lnext();
	auto e2 = e1.Onext();

	auto mark1 = Left(e1)->mark;
	auto div1 = s.splitFace(e1, e1.Lnext().Lnext());
	Right(div1)->mark = Left(div1)->mark = mark1;

	e2 = div1.Onext();
	auto mark2 = Left(e2)->mark;
	auto div2 = s.splitFace(e2, e2.Lnext().Lnext());
	Right(div2)->mark = Left(div2)->mark = mark2;


	auto X = Org(div1);
	e = div2.Lnext();
	auto first = Dest(e);
	do {
		auto t = e.Oprev();
		if (!e.data().fixed && rightOf(Dest(t), e) && incircle(Org(e), Dest(t), Dest(e), X)) {
			e = swap_wf(s, e);
			assert(X == Dest(e));
			e = e.Oprev();
		}
		else if (Org(e) == first)
			break;
		else
			e = e.Onext().Lprev();
	} while (true);


	return e1;
}

EdgeRef splitRegularEdge(Subdivision& s, EdgeRef e)
{
	assert(e.data().fixed);
	assert(e.Sym().data().fixed);
	assert(!e.data().boundary);
	assert(!e.Sym().data().boundary);

	assert(Left(e)->mark);
	assert(Right(e)->mark);

	auto mark = Left(e)->mark;

	auto e1 = s.splitVertex(e, e.Oprev(), midpoint(e));
	e1.data().fixed = true;
	e1.Sym().data().fixed = true;


	e1 = e1.Lnext();

	auto div1 = s.splitFace(e1, e1.Lnext().Lnext());
	Left(div1)->mark = mark;
	Right(div1)->mark = mark;

	auto e2 = div1.Onext();

	auto div2 = s.splitFace(e2, e2.Lnext().Lnext());
	Left(div2)->mark = mark;
	Right(div2)->mark = mark;

	VertexRef X = Org(div2);
	e = div2.Lnext();
	auto first = Dest(e);
	do {
		auto t = e.Oprev();
		if (!e.data().fixed && rightOf(Dest(t), e) && incircle(Org(e), Dest(t), Dest(e), X)) {
			e = swap_wf(s, e);
			assert(X == Dest(e));
			e = e.Oprev();
		}
		else if (Org(e) == first)
			break;
		else
			e = e.Onext().Lprev();
	} while (true);

	return e1;
}

void splitEdges(Subdivision & dt)
{
	bool all_done = false;
	while (!all_done)
	{
		all_done = true;
		for (auto qref = dt.edges.begin(); qref != dt.edges.end(); ++qref)
		{
			EdgeRef e(qref);
			if (e.data().var.which() != 0)
				e = e.Rot();

			if (!e.data().fixed) continue;
			if (encroaches(e, Dest(e.Onext())) || encroaches(e, Dest(e.Oprev())))
			{
				if (e.data().boundary)
					splitBoundaryEdge(dt, e);
				else {
					assert(e.data().fixed);
					splitRegularEdge(dt, e);
				}
				all_done = false;
				break;
			}
		}
	}
}

bool inside(Subdivision & s, Point x)
{
	auto f = s.outer_face;
	auto e = f->bounds;

	if (rightOf(x, e)
		&& rightOf(x, e.Lnext())
		&& rightOf(x, e.Lnext().Lnext())
		)
		return true;
	return false;
}

VertexRef insertMeshSite(Subdivision& s, Point x)
{
	if (!inside(s, x)) {
		std::cerr << "fatal error: circumcenter is outside subdivision\n";
		return s.vertices.end();
	}

	EdgeRef e = locate(s, x);
	if (!e)
		return s.vertices.end();
	if (x == Org(e)->point || x == Dest(e)->point) // ignore
		return s.vertices.end();

	bool on_edge = onEdge(x, e);
	if (on_edge) {
		if (e.data().fixed) // 'x' encroaches 'e'
		{
			//std::cout << "on fixed edge " << e.data().boundary << '\n';
			if (e.data().boundary)
				splitBoundaryEdge(s, e);
			else {
				assert(e.data().fixed);
				splitRegularEdge(s, e);
			}
			return s.vertices.end();
		}
	}

	// check for conflicts
	for (auto qref = s.edges.begin(); qref != s.edges.end(); ++qref)
	{
		EdgeRef e(qref);
		if (e.data().var.which() != 0)
			e = e.Rot();

		if (!e.data().fixed) continue;
		if (encroaches(e, x))
		{
			//std::cout << "conflict with edge " << e.data().boundary << '\n';
			if (e.data().boundary)
				splitBoundaryEdge(s, e);
			else {
				assert(e.data().fixed);
				splitRegularEdge(s, e);
			}
			return s.vertices.end();
		}
	}

	//std::cout << "normal insertion\n";
	if (on_edge) {
		assert(Left(e)->mark == 1);
		assert(Left(e)->mark == Right(e)->mark);
		e = e.Oprev();
		s.joinFace(e.Onext());
		Left(e)->mark = 1;
	}

	VertexRef first = Org(e);
	EdgeRef base = s.splitVertex(e, e, x);
	VertexRef X = Dest(base);
	do {
		base = s.splitFace(e.Lnext(), base.Sym());
		Left(base)->mark = Right(base)->mark = 1;
		e = base.Oprev();
	} while (Dest(e) != first);

	assert(Dest(e) == first);
	assert(Dest(e.Onext()) == X);
	do {
		auto t = e.Oprev();
		if (!e.data().fixed
			&& rightOf(Dest(t), e)
			&& incircle(Org(e), Dest(t), Dest(e), X)) {
			e = swap_wf(s, e);
			assert(X == Dest(e));
			e = e.Oprev();
		}
		else if (Org(e) == first)
			break;
		else
			e = e.Onext().Lprev();
	} while (true);

	X->circumcenter = true;
	return X;
}


bool eliminate_worst_triangle(Subdivision& dt, double min_ratio)
{
	FaceRef face; double ratio;
	std::tie(face, ratio) = find_worst(dt);
	if (ratio <= min_ratio) return false;
	auto e = face->bounds;
	Point cc = circumCenter(Org(e)->point, Dest(e)->point, Dest(e.Onext())->point);
	insertMeshSite(dt, cc);
	return true;
}

bool eliminate_worst_triangle(Subdivision& dt, double min_ratio, double min_area)
{
	FaceRef face; double ratio;
	std::tie(face, ratio) = find_worst(dt);
	if (ratio <= min_ratio) {
		double area;
		std::tie(face, area) = find_biggest(dt);
		if (area <= min_area)
			return false;
	}
	auto e = face->bounds;
	Point cc = circumCenter(Org(e)->point, Dest(e)->point, Dest(e.Onext())->point);
	insertMeshSite(dt, cc);
	return true;
}

bool eliminate_bad_triangle(Subdivision& dt, double min_ratio)
{
	auto face = find_bad(dt, min_ratio);
	if (face == dt.faces.end()) return false;
	auto e = face->bounds;
	Point cc = circumCenter(Org(e)->point, Dest(e)->point, Dest(e.Onext())->point);
	insertMeshSite(dt, cc);
	return true;
}

bool eliminate_bad_triangle(Subdivision& dt, double min_ratio, double min_area)
{
	auto face = find_bad(dt, min_ratio, min_area);
	if (face == dt.faces.end()) return false;
	auto e = face->bounds;
	Point cc = circumCenter(Org(e)->point, Dest(e)->point, Dest(e.Onext())->point);
	insertMeshSite(dt, cc);
	return true;
}


void eliminate_triangle(Subdivision& dt, FaceRef face)
{
	auto e = face->bounds;
	Point cc = circumCenter(Org(e)->point, Dest(e)->point, Dest(e.Onext())->point);
	insertMeshSite(dt, cc);
}

std::tuple<FaceRef, double> find_worst(Subdivision & dt)
{
	FaceRef worst_face = dt.faces.end();
	double max_ratio = -1.0;

	for (auto face = dt.faces.begin(); face != dt.faces.end(); ++face)
	{
		if (!face->mark || face == dt.outer_face) continue;
		auto e = face->bounds;
		double ratio = quality_measure(Org(e)->point, Dest(e)->point, Dest(e.Onext())->point);
		if (ratio > max_ratio) {
			max_ratio = ratio;
			worst_face = face;
		}
	}
	return{worst_face,max_ratio};
}

std::tuple<FaceRef, double> find_biggest(Subdivision & dt)
{
	FaceRef biggest_face = dt.faces.end();
	double max_area = -1.0;

	for (auto face = dt.faces.begin(); face != dt.faces.end(); ++face)
	{
		if (!face->mark || face == dt.outer_face) continue;
		auto e = face->bounds;
		double area = triangleArea(Org(e)->point, Dest(e)->point, Dest(e.Onext())->point);
		if (area > max_area) {
			max_area = area;
			biggest_face = face;
		}
	}
	return{biggest_face,max_area};
}

std::tuple<FaceRef, double> find_smallest(Subdivision & dt)
{
	FaceRef smallest_face = dt.faces.end();
	double min_area = std::numeric_limits<double>::max();

	for (auto face = dt.faces.begin(); face != dt.faces.end(); ++face)
	{
		if (!face->mark || face == dt.outer_face) continue;
		auto e = face->bounds;
		double area = triangleArea(Org(e)->point, Dest(e)->point, Dest(e.Onext())->point);
		if (area < min_area) {
			min_area = area;
			smallest_face = face;
		}
	}
	return{smallest_face,min_area};
}

FaceRef find_bad(Subdivision & dt, double min_ratio, double min_area)
{
	FaceRef bad_face = dt.faces.end();

	for (auto face = dt.faces.begin(); face != dt.faces.end(); ++face)
	{
		if (!face->mark || face == dt.outer_face) continue;
		auto e = face->bounds;
		double ratio = quality_measure(Org(e)->point, Dest(e)->point, Dest(e.Onext())->point);
		double area = triangleArea(Org(e)->point, Dest(e)->point, Dest(e.Onext())->point);
		if (ratio > min_ratio) {
			min_ratio = ratio;
			bad_face = face;
			break;
		}
		else if (area > min_area) {
			min_area = area;
			bad_face = face;
			break;
		}
	}
	return bad_face;
}

FaceRef find_bad(Subdivision & dt, double max_ratio)
{
	FaceRef worst_face = dt.faces.end();

	for (auto face = dt.faces.begin(); face != dt.faces.end(); ++face)
	{
		if (!face->mark || face == dt.outer_face) continue;
		auto e = face->bounds;
		double ratio = quality_measure(Org(e)->point, Dest(e)->point, Dest(e.Onext())->point);
		if (ratio > max_ratio) {
			worst_face = face;
			break;
		}
	}
	return worst_face;
}

void ruppert_refinement(Subdivision & dt, double min_ratio, int max_iters)
{
	splitEdges(dt);
	int iters = 0;
	while (iters < max_iters && eliminate_worst_triangle(dt, min_ratio))
		++iters;
	std::cout << iters << '\n';
}

void ruppert_refinement(Subdivision & dt, double min_ratio, double min_area, int max_iters)
{
	splitEdges(dt);
	int iters = 0;
	while (iters < max_iters && eliminate_worst_triangle(dt, min_ratio, min_area))
		++iters;
	std::cout << iters << '\n';
}

void insertClosedLoop(Subdivision& dt, std::vector<Point> const & hole)
{
	std::vector<VertexRef> inserted;
	inserted.reserve(hole.size());
	for (Point const& p : hole)
		inserted.push_back(insertSite(dt, p));
	for (unsigned i = 0; i < inserted.size(); ++i) {
		auto ins = insertEdge(dt, inserted[i], inserted[(i + 1) % inserted.size()]);
		ins.data().boundary = true;
		ins.Sym().data().boundary = true;
	}
}

EdgeRef clip_delaunay_ear(Subdivision& dt, EdgeRef e)
{
	auto mark = Left(e)->mark;
	if (Dest(e.Lnext().Lnext()) == Org(e))
		return e;
	if (!leftOf(Dest(e.Lnext()), e))
		return e.Lnext();

	VertexRef a = Org(e), b = Dest(e), c = Dest(e.Lnext());
	auto ei = e.Lnext().Lnext();
	do {
		if (incircle(a, b, c, Dest(ei)))
			return e.Lnext();
		ei = ei.Lnext();
	} while (Dest(ei) != Org(e));

	auto new_e = dt.splitFace(e, e.Lnext().Lnext());
	Left(new_e)->mark = Right(new_e)->mark = mark;
	return new_e;
}

void deleteSite_wf(Subdivision & dt, VertexRef v)
{
	assert(v->circumcenter);
	auto e = v->leaves;
	auto b = e.Lnext();

	do {
		e = e.Onext();
		dt.joinFace(e.Oprev());
		Right(e)->mark = 1;
	} while (e != e.Onext());
	dt.joinVertex(e.Sym());

	auto res = b;
	do {
		b = res;
		res = clip_delaunay_ear(dt, b);
	} while (res != b);
}

VertexRef insertSite_wf(Subdivision& s, Point x, EdgeRef e)
{
	if (x == Org(e)->point || x == Dest(e)->point) // ignore
		return s.vertices.end();
	else if (onEdge(x, e)) {
		e = e.Oprev();
		s.joinFace(e.Onext());
	}

	// connect
	VertexRef first = Org(e);
	EdgeRef base = s.splitVertex(e, e, x);
	VertexRef X = Dest(base);
	X->circumcenter = true;
	do {
		base = s.splitFace(e.Lnext(), base.Sym());
		Left(base)->mark = Right(base)->mark = 1;
		e = base.Oprev();
	} while (Dest(e) != first);

	// inspect edges
	do {
		auto t = e.Oprev();
		if (!e.data().fixed && rightOf(Dest(t), e) && incircle(Org(e), Dest(t), Dest(e), X)) {
			e = swap_wf(s, e);
			assert(X == Dest(e));
			e = e.Oprev();
		}
		else if (Org(e) == first)
			break;
		else
			e = e.Onext().Lprev();
	} while (true);
	return X;
}

bool chew_2nd_eliminate_worst(Subdivision& dt, double min_ratio)
{
	FaceRef face; double ratio;

	std::tie(face, ratio) = find_worst(dt);
	if (ratio <= min_ratio)
		return false;
	auto e = face->bounds;
	Point c = circumCenter(Org(e)->point, Dest(e)->point, Dest(e.Onext())->point);

	bool found_edge{true};
	do
	{
		if (rightOf(c, e))
			e = e.Sym();
		else if (!rightOf(c, e.Onext()))
			e = e.Onext();
		else if (!rightOf(c, e.Dprev()))
			e = e.Dprev();
		else {
			found_edge = false;
			break;
		}
	} while (!e.data().fixed);

	// here c can be on edge or inside a triangle
	if (e.data().fixed && onEdge(c, e)) // and if 'c' is on a fixed edge, then work with edge
		found_edge = true;

	if (!found_edge) {
		auto v = insertSite_wf(dt, c, e);
		v->circumcenter = true;
		return true;
	}

	// found encroached edge
	// delete all encroaching vertices
	bool all_done = false;
	while (!all_done)
	{
		all_done = true;
		for (auto v = dt.vertices.begin(); v != dt.vertices.end(); ++v)
		{
			if (v->circumcenter && encroaches(e, v->point)) {
				deleteSite_wf(dt, v);
				all_done = false;
				break;
			}
		}
	}
	// split e
	if (e.data().boundary)
		splitBoundaryEdge(dt, e);
	else
		splitRegularEdge(dt, e);

	return true;
}

bool chew_2nd_eliminate_worst(Subdivision& dt, double min_ratio, double min_area)
{
	FaceRef face; double ratio;

	std::tie(face, ratio) = find_worst(dt);
	if (ratio <= min_ratio)
	{
		double area;
		std::tie(face, area) = find_biggest(dt);
		if (area <= min_area)
			return false;
	}
	
	auto e = face->bounds;
	Point c = circumCenter(Org(e)->point, Dest(e)->point, Dest(e.Onext())->point);

	bool found_edge{true};
	do
	{
		if (rightOf(c, e))
			e = e.Sym();
		else if (!rightOf(c, e.Onext()))
			e = e.Onext();
		else if (!rightOf(c, e.Dprev()))
			e = e.Dprev();
		else {
			found_edge = false;
			break;
		}
	} while (!e.data().fixed);

	// here c can be on edge or inside a triangle
	if (e.data().fixed && onEdge(c, e)) // and if 'c' is on a fixed edge, then work with edge
		found_edge = true;

	if (!found_edge) {
		auto v = insertSite_wf(dt, c, e);
		v->circumcenter = true;
		return true;
	}
	// found encroached edge
	// delete all encroaching vertices
	bool all_done = false;
	while (!all_done)
	{
		all_done = true;
		for (auto v = dt.vertices.begin(); v != dt.vertices.end(); ++v)
		{
			if (v->circumcenter && encroaches(e, v->point)) {
				deleteSite_wf(dt, v);
				all_done = false;
				break;
			}
		}
	}
	// split e
	if (e.data().boundary) {
		splitBoundaryEdge(dt, e);
	}
	else
		splitRegularEdge(dt, e);

	return true;
}

void chew_2nd_refinement(Subdivision& dt, double min_ratio, int iters)
{
	int i = 0;
	while (i++ < iters && chew_2nd_eliminate_worst(dt, min_ratio))
		;
	std::cout << "iters: " << i << '\n';
}

void chew_2nd_refinement(Subdivision& dt, double min_ratio, double min_area, int iters)
{
	int i = 0;
	while (i++ < iters && chew_2nd_eliminate_worst(dt, min_ratio, min_area))
		;
	std::cout << "iters: " << i << '\n';
}

bool off_center_correction(Subdivision& dt, VertexRef v, double min_angle, double q)
{
	assert(v->circumcenter);
	using boost::math::double_constants::pi;

	auto e = v->leaves;
	auto end = e;
	do {
		Point& A = Org(e)->point;
		Point& B = Dest(e)->point;
		Point& C = Dest(e.Onext())->point;
		Point AB = B - A, AC = C - A;
		double angle = acos(AB*AC / (norm(AB) * norm(AC))) * 180.0 / pi;
		if (angle < min_angle)
		{
			Point BA = A - B, BC = C - B;
			double t_P = (BA*BC) / (BC*BC);
			Point BP = BC*t_P;
			Point AP = BP - BA;
			A = A + AP*q;
			return false;
		}

		e = e.Onext();
	} while (e != end);
	return true;
}

double ratio_to_angle(double ratio) 
{
	return asin(1.0 / (2.0*ratio)) / boost::math::double_constants::pi * 180.0;
}

bool chew_2nd_eliminate_worst_correction(Subdivision& dt, double min_ratio, double q)
{
	FaceRef face; double ratio;

	std::tie(face, ratio) = find_worst(dt);
	if (ratio <= min_ratio)
		return false;
	auto e = face->bounds;
	Point c = circumCenter(Org(e)->point, Dest(e)->point, Dest(e.Onext())->point);

	bool found_edge{true};
	do
	{
		if (rightOf(c, e))
			e = e.Sym();
		else if (!rightOf(c, e.Onext()))
			e = e.Onext();
		else if (!rightOf(c, e.Dprev()))
			e = e.Dprev();
		else {
			found_edge = false;
			break;
		}
	} while (!e.data().fixed);

	// here c can be on edge or inside a triangle
	if (e.data().fixed && onEdge(c, e)) // and if 'c' is on a fixed edge, then work with edge
		found_edge = true;

	if (!found_edge) {
		auto v = insertSite_wf(dt, c, e);
		v->circumcenter = true;
		off_center_correction(dt, v, ratio_to_angle(min_ratio), q);
		return true;
	}

	// found encroached edge
	// delete all encroaching vertices
	bool all_done = false;
	while (!all_done)
	{
		all_done = true;
		for (auto v = dt.vertices.begin(); v != dt.vertices.end(); ++v)
		{
			if (v->circumcenter && encroaches(e, v->point)) {
				deleteSite_wf(dt, v);
				all_done = false;
				break;
			}
		}
	}
	// split e
	if (e.data().boundary)
		splitBoundaryEdge(dt, e);
	else
		splitRegularEdge(dt, e);

	return true;
}

bool chew_2nd_eliminate_worst_correction(Subdivision& dt, double min_ratio,
	double min_area, double q)
{
	FaceRef face; double ratio;

	std::tie(face, ratio) = find_worst(dt);
	if (ratio <= min_ratio)
	{
		double area;
		std::tie(face, area) = find_biggest(dt);
		if (area <= min_area)
			return false;
	}

	auto e = face->bounds;
	Point c = circumCenter(Org(e)->point, Dest(e)->point, Dest(e.Onext())->point);

	bool found_edge{true};
	do
	{
		if (rightOf(c, e))
			e = e.Sym();
		else if (!rightOf(c, e.Onext()))
			e = e.Onext();
		else if (!rightOf(c, e.Dprev()))
			e = e.Dprev();
		else {
			found_edge = false;
			break;
		}
	} while (!e.data().fixed);

	// here c can be on edge or inside a triangle
	if (e.data().fixed && onEdge(c, e)) // and if 'c' is on a fixed edge, then work with edge
		found_edge = true;

	if (!found_edge) {
		auto v = insertSite_wf(dt, c, e);
		v->circumcenter = true;
		off_center_correction(dt, v, ratio_to_angle(min_ratio), q);
		return true;
	}

	// found encroached edge
	// delete all encroaching vertices
	bool all_done = false;
	while (!all_done)
	{
		all_done = true;
		for (auto v = dt.vertices.begin(); v != dt.vertices.end(); ++v)
		{
			if (v->circumcenter && encroaches(e, v->point)) {
				deleteSite_wf(dt, v);
				all_done = false;
				break;
			}
		}
	}
	// split e
	if (e.data().boundary)
		splitBoundaryEdge(dt, e);
	else
		splitRegularEdge(dt, e);

	return true;
}

void chew_2nd_refinement_alper(Subdivision& dt, double q, double min_ratio, int iters)
{
	int i = 0;
	while (i++ < iters && chew_2nd_eliminate_worst_correction(dt, min_ratio, q))
		;
	//std::cout << "iters: " << i << '\n';
}

void chew_2nd_refinement_alper(Subdivision& dt, double q, double min_ratio, double min_area, int iters)
{
	int i = 0;
	while (i++ < iters && chew_2nd_eliminate_worst_correction(dt, min_ratio, min_area, q))
		;
	//std::cout << "iters: " << i << '\n';
}