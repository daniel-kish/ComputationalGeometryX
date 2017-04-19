#include "mesh.h"
#include <vector>
#include "Subdivision.h"
#include "delaunay.h"

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

Point midpoint(EdgeRef e)
{
	return (Org(e)->point + Dest(e)->point) * 0.5;
}

EdgeRef swap_wf(Subdivision &s, EdgeRef e)
{
	auto a = e.Oprev().Lnext();
	auto b = e.Lnext().Lnext();

	assert(Left(e)->mark == Right(e)->mark);
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
	assert(e.data().fixed);
	assert(Left(e)->mark != Right(e)->mark);

	if (Left(e)->mark)
		e = e.Sym();

	auto e1 = s.splitVertex(e, e.Oprev(), midpoint(e));
	e1.data().boundary = true;
	e1.data().fixed = true;
	
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


	return EdgeRef{};
}