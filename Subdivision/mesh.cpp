#include "mesh.h"
#include <vector>
#include "Subdivision.h"
#include "delaunay.h"
#include "geom.h"

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

	auto div1 = s.splitFace(e1, e1.Lnext().Lnext());
	Right(div1)->mark = Left(div1)->mark;


	e2 = div1.Onext();
	auto div2 = s.splitFace(e2, e2.Lnext().Lnext());
	Right(div2)->mark = Left(div2)->mark;

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
			e = swap_wf(s,e);
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

void eliminate_worst_triangle(Subdivision& dt)
{
	FaceRef face; double ratio;
	std::tie(face,ratio) = find_worst(dt);
	auto e = face->bounds;
	Point cc = circumCenter(Org(e)->point, Dest(e)->point, Dest(e.Onext())->point);
	insertMeshSite(dt, cc);
}


void eliminate_triangle(Subdivision& dt, FaceRef face)
{
	auto e = face->bounds;
	Point cc = circumCenter(Org(e)->point, Dest(e)->point, Dest(e.Onext())->point);
	insertMeshSite(dt, cc);
}

std::tuple<FaceRef,double> find_worst(Subdivision & dt)
{
	double max_ratio = -1.0;
	FaceRef worst_face = dt.faces.end();

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
	return {worst_face,max_ratio};
}
