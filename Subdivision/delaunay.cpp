#include "delaunay.h"
#include "predicates.h"


void swap(EdgeRef e)
{
	auto a = e.Oprev();
	auto b = e.Sym().Oprev();

	if (Org(e)->leaves == e)
		Org(e)->leaves = e.Onext();
	if (Dest(e)->leaves == e.Sym())
		Dest(e)->leaves = e.Sym().Onext();

	splice(e, a); splice(e.Sym(), b);
	splice(e, a.Lnext()); splice(e.Sym(), b.Lnext());
	
	e.data().var = a.Sym().data().var;
	e.Sym().data().var = b.Sym().data().var;
}

bool rightOf(Point p, EdgeRef e)
{
	return orient2d(p, Dest(e)->point, Org(e)->point) > 0.0;
}

bool leftOf(Point p, EdgeRef e)
{
	return orient2d(p, Org(e)->point, Dest(e)->point) > 0.0;
}

bool rightOf(VertexRef v, EdgeRef e)
{
	return orient2d(v->point, Dest(e)->point, Org(e)->point) > 0.0;
}

bool leftOf(VertexRef v, EdgeRef e)
{
	return orient2d(v->point, Org(e)->point, Dest(e)->point) > 0.0;
}

bool leftOf(VertexRef v, std::pair<VertexRef, VertexRef> e)
{
	return orient2d(e.first->point, e.second->point, v->point) > 0.0;
}

bool rightOf(VertexRef v, std::pair<VertexRef, VertexRef> e)
{
	return orient2d(e.first->point, v->point, e.second->point) > 0.0;
}

bool incircle(VertexRef a, VertexRef b, VertexRef c, VertexRef d) {
	double det = incircle((double*)&(a->point), (double*)&(b->point),
		(double*)&(c->point), (double*)&(d->point));
	return det > 0.0;
}

bool onEdge(Point c, EdgeRef e) {

	Point const& a = Org(e)->point;
	Point const& b = Dest(e)->point;

	double det = orient2d(a, b, c);
	if (det != 0.0)
		return false;

	Point ab = b - a;
	Point ac = c - a;
	double t_p = scp(ac, ab) / scp(ab, ab);

	if (t_p > 0.0 && t_p < 1.0)
		return true;
	
	return false;
}

std::tuple<Subdivision, EdgeRef, EdgeRef>
delaunay_dnc(std::vector<Point>::iterator b, std::vector<Point>::iterator e)
{
	if (std::distance(b, e) == 2) {
		Subdivision s(*b, *(b + 1));
		EdgeRef e(s.edges.begin());

		return std::make_tuple(std::move(s), e, e.Sym());
	}
	else if (std::distance(b, e) == 3) {
		auto s1 = b;
		auto s2 = b + 1;
		auto s3 = b + 2;

		Subdivision s(*b, *(b + 1));
		EdgeRef e1(s.edges.begin());
		EdgeRef e2 = s.add_vertex(e1, *(b + 2));

		double det = orient2d(*s1, *s2, *s3);
		if (det > 0.0) {
			s.connect(e2, e1);
			return std::make_tuple(std::move(s), e1, e2.Sym());
		}
		else if (det < 0.0) {
			EdgeRef e3 = s.connect(e2, e1);
			return std::make_tuple(std::move(s), e3.Sym(), e3);
		}
		else
			return std::make_tuple(std::move(s), e1, e2.Sym());
	}
	// else if S.size() >= 4 :

	auto mid{b}; std::advance(mid, std::distance(b, e) / 2);
	Subdivision L, R;
	EdgeRef ldo, ldi, rdi, rdo;
	std::tie(L, ldo, ldi) = delaunay_dnc(b, mid);
	std::tie(R, rdi, rdo) = delaunay_dnc(mid, e);

	L.merge(R); // physical merge

	// logical merge:
	do {
		if (leftOf(Org(rdi), ldi)) ldi = ldi.Lnext();
		else if (rightOf(Org(ldi), rdi)) rdi = rdi.Rprev();
		else break;
	} while (true);

	auto basel = L.connect(rdi.Sym(), ldi);
	if (Org(ldi) == Org(ldo)) ldo = basel.Sym();
	if (Org(rdi) == Org(rdo)) rdo = basel;

	auto valid = [&basel](EdgeRef e) {return rightOf(Dest(e), basel); };
	do // merge loop
	{
		auto lcand = basel.Sym().Onext();
		if (valid(lcand)) {
			while (incircle(Dest(basel), Org(basel), Dest(lcand), Dest(lcand.Onext())))
			{
				lcand = lcand.Onext();
				L.deleteEdge(lcand.Oprev());
			}
		}
		auto rcand = basel.Oprev();
		if (valid(rcand)) {
			while (incircle(Dest(basel), Org(basel), Dest(rcand), Dest(rcand.Oprev())))
			{
				rcand = rcand.Oprev();
				L.deleteEdge(rcand.Onext());
			}
		}
		if (!valid(lcand) && !valid(rcand)) break; // we've reached the top

		if (!valid(lcand) ||
			valid(rcand) && incircle(Dest(lcand), Org(lcand), Org(rcand), Dest(rcand)))
			basel = L.connect(rcand, basel.Sym());
		else
			basel = L.connect(basel.Sym(), lcand.Sym());
	} while (true);
	return std::make_tuple(std::move(L), ldo, rdo);
}

EdgeRef locate(Subdivision & s, Point x)
{
	std::size_t N = s.edges.size(); // ?
	EdgeRef e(s.edges.begin());

	do {
		N--;
		if (rightOf(x, e))
			e = e.Sym();
		else if (!rightOf(x, e.Onext()))
			e = e.Onext();
		else if (!rightOf(x, e.Dprev()))
			e = e.Dprev();
		else {
			//std::cout << s.edges.size() - N << '\n';
			return e;
		}
	} while (N);
	return EdgeRef{};
}

EdgeRef locate(Subdivision& s, Point x, EdgeRef e)
{
	std::size_t N = s.edges.size(); // ?
	do {
		N--;
		if (rightOf(x, e))
			e = e.Sym();
		else if (!rightOf(x, e.Onext()))
			e = e.Onext();
		else if (!rightOf(x, e.Dprev()))
			e = e.Dprev();
		else {
			//std::cout << s.edges.size() - N << '\n';
			return e;
		}
	} while (N);
	return EdgeRef{};
}

VertexRef insertSite(Subdivision& s, Point x)
{
	EdgeRef e = locate(s, x);
	if (!e) 
		return s.vertices.end();
	
	if (x == Org(e)->point || x == Dest(e)->point) // ignore
		return s.vertices.end();
	else if (onEdge(x, e)) {
		e = e.Oprev();
		s.deleteEdge(e.Onext());
	}
	
	// connect
	VertexRef first = Org(e);
	EdgeRef base = s.add_vertex(e.Lprev(), x);
	VertexRef X = Dest(base);
	
	do {
		base = s.connect(e, base.Sym());
		e = base.Oprev();
	} while (Dest(e) != first);

	assert(Dest(e) == first);
	assert(Dest(e.Onext()) == X);
	do {
		auto t = e.Oprev();
		if (!e.data().fixed && rightOf(Dest(t), e) && incircle(Org(e), Dest(t), Dest(e), X)) {
			swap(e);
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

VertexRef insertSite(Subdivision& s, Point x, EdgeRef start)
{
	EdgeRef e = locate(s, x, start);
	if (!e)
		return s.vertices.end();

	if (x == Org(e)->point || x == Dest(e)->point) // ignore
		return s.vertices.end();
	else if (onEdge(x, e)) {
		e = e.Oprev();
		s.deleteEdge(e.Onext());
	}

	// connect
	VertexRef first = Org(e);
	EdgeRef base = s.add_vertex(e.Lprev(), x);
	VertexRef X = Dest(base);

	do {
		base = s.connect(e, base.Sym());
		e = base.Oprev();
	} while (Dest(e) != first);

	do {
		auto t = e.Oprev();
		if (!e.data().fixed && rightOf(Dest(t), e) && incircle(Org(e), Dest(t), Dest(e), X)) {
			swap(e);
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

void insertSiteSequence(Subdivision & s, std::vector<Point> seq)
{
	if (seq.empty())
		return;

	VertexRef v = insertSite(s, seq[0]);
	if (v == s.vertices.end())
		return;

	EdgeRef close = v->leaves;
	for (unsigned int i = 1; i < seq.size(); ++i) {
		v = insertSite(s, seq[i], close);
		close = v->leaves;
	}
}

void triangulatePseudoPolygon(Subdivision& s, EdgeRef c)
{
	VertexRef a{Org(c)}, b{Dest(c)};

	if (Dest(c.Lnext().Lnext()) == a)
		return;

	auto begin = c.Lnext();
	auto end = c.Lprev();
	auto last = end.Lprev();
	auto e = begin;
	do {
		auto p = begin;
		bool fnd = true;
		do {
			if (incircle(a, b, Dest(e), Dest(p))) {
				fnd = false;
				break;
			}
			p = p.Lnext();
		} while (p != end);

		if (fnd) {
			bool done = false;
			if (e != begin) {
				auto r = s.connect(e, begin);
				triangulatePseudoPolygon(s, r);
				done = true;
			}
			if (e != last) {
				auto l = s.connect( done? e.Dprev():e , c);
				l = l.Sym();
				triangulatePseudoPolygon(s, l);
			}
			break;
		}
		e = e.Lnext();
	} while (e != end);

}

EdgeRef insertEdge(Subdivision& s, VertexRef a, VertexRef b)
{
	EdgeRef e = a->leaves;

	EdgeRef p = e;
	do {
		if (Dest(e) == b) {
			e.data().fixed = true;
			e.Sym().data().fixed = true;
			return e;
		}
		e = e.Onext();
	} while (e != p);

	p = e;
	do {
		if (rightOf(Dest(e), {a,b}) && leftOf(Dest(e.Onext()), {a,b}))
			break;
		e = e.Onext();
	} while (e != p);
	
	assert(Dest(e.Onext()) == Dest(e.Lnext()));

	EdgeRef first = e.Onext().Sym();
	assert(Dest(first) == a);

	e = e.Lnext();
	p = e;
	do {
		if (Dest(e.Oprev()) == b) {
			s.deleteEdge(e);
			break;
		}
		VertexRef n = Dest(e.Oprev());
		double det = orient2d(a->point, b->point, n->point);
		if (det > 0.0) // leftOf
			p = e.Oprev();
		else if (det < 0.0) // rightOf
			p = e.Dnext();
		s.deleteEdge(e);
		e = p;
	} while (true);

	auto ei = first;
	do
		ei = ei.Lnext();
	while (Dest(ei) != b);

	EdgeRef last = ei;
	
	auto c = s.connect(first,last.Lnext());

	c.data().fixed = true;
	c.Sym().data().fixed = true;

	triangulatePseudoPolygon(s, c);
	triangulatePseudoPolygon(s, c.Sym());

	return c;
}
