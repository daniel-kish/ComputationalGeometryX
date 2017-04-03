#include "delaunay.h"
#include "predicates.h"


void swap(Edge e)
{
	auto a = e.Oprev();
	auto b = e.Sym().Oprev();
	splice(e, a); splice(e.Sym(), b);
	splice(e, a.Lnext()); splice(e.Sym(), b.Lnext());
	e.data() = a.Sym().data();
	e.Sym().data() = b.Sym().data();
}

bool rightOf(Point p, Edge e)
{
	return orient2d(p, Dest(e)->point, Org(e)->point) > 0.0;
}

bool leftOf(Point p, Edge e)
{
	return orient2d(p, Org(e)->point, Dest(e)->point) > 0.0;
}

bool rightOf(Vertex v, Edge e)
{
	return orient2d(v->point, Dest(e)->point, Org(e)->point) > 0.0;
}

bool leftOf(Vertex v, Edge e)
{
	return orient2d(v->point, Org(e)->point, Dest(e)->point) > 0.0;
}

bool incircle(Vertex a, Vertex b, Vertex c, Vertex d) {
	double det = incircle((double*)&(a->point), (double*)&(b->point),
		(double*)&(c->point), (double*)&(d->point));
	return det > 0.0;
}

bool onEdge(Point c, Edge e) {

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

std::tuple<Subdivision, Edge, Edge>
delaunay_dnc(std::vector<Point>::iterator b, std::vector<Point>::iterator e)
{
	if (std::distance(b, e) == 2) {
		Subdivision s(*b, *(b + 1));
		Edge e(s.edges.begin());

		return std::make_tuple(std::move(s), e, e.Sym());
	}
	else if (std::distance(b, e) == 3) {
		auto s1 = b;
		auto s2 = b + 1;
		auto s3 = b + 2;

		Subdivision s(*b, *(b + 1));
		Edge e1(s.edges.begin());
		Edge e2 = s.add_vertex(e1, *(b + 2));


		double det = orient2d(*s1, *s2, *s3);
		if (det > 0.0) {
			s.connect(e2, e1); //connect(s, b, a);
			return std::make_tuple(std::move(s), e1, e2.Sym());
		}
		else if (det < 0.0) {
			Edge e3 = s.connect(e2, e1);
			return std::make_tuple(std::move(s), e3.Sym(), e3);
		}
		else
			return std::make_tuple(std::move(s), e1, e2.Sym());
	}
	// else if S.size() >= 4 :

	auto mid{b}; std::advance(mid, std::distance(b, e) / 2);
	Subdivision L, R;
	Edge ldo, ldi, rdi, rdo;
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

	auto valid = [&basel](Edge e) {return rightOf(Dest(e), basel); };
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

std::tuple<Edge,bool> locate(Subdivision & s, Point x)
{
	std::size_t N = s.edges.size(); // ?
	Edge e(s.edges.begin());

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
			return{e,true};
		}
	} while (N);
	return {e,false};
}

std::tuple<Edge, bool> locate(Subdivision & s, Point x, Edge e)
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
			return{e,true};
		}
	} while (N);
	return{e,false};
}

std::tuple<Edge, bool> insertSite(Subdivision& s, Point x)
{
	Edge e; bool fnd;
	std::tie(e,fnd) = locate(s, x);
	if (!fnd) 
		return {e,fnd};
	
	if (x == Org(e)->point || x == Dest(e)->point) // ignore
		return {e,false};
	else if (onEdge(x, e)) {
		std::cout << "here\n";
		e = e.Oprev();
		s.deleteEdge(e.Onext());
	}
	
	// connect
	Vertex first = Org(e);
	Edge base = s.add_vertex(e.Lprev(), x);
	Vertex X = Dest(base);
	
	do {
		base = s.connect(e, base.Sym());
		e = base.Oprev();
	} while (Dest(e) != first);

	Edge close;
	do {
		auto t = e.Oprev();
		if (rightOf(Dest(t), e) && incircle(Org(e), Dest(t), Dest(e), X)) {
			swap(e);
			assert(X == Dest(e));
			e = e.Oprev();
		}
		else if (Org(e) == first)
			break;
		else {
			close = e;
			e = e.Onext().Lprev();
		}
	} while (true);
	close = close.Onext().Sym();
	return {close,true};
}

std::tuple<Edge, bool> insertSite(Subdivision& s, Point x, Edge start)
{
	Edge e; bool fnd;
	std::tie(e, fnd) = locate(s, x, start);
	if (!fnd)
		return{e,fnd};

	if (x == Org(e)->point || x == Dest(e)->point) // ignore
		return{e,false};
	else if (onEdge(x, e)) {
		std::cout << "here\n";
		e = e.Oprev();
		s.deleteEdge(e.Onext());
	}

	// connect
	Vertex first = Org(e);
	Edge base = s.add_vertex(e.Lprev(), x);
	Vertex X = Dest(base);

	do {
		base = s.connect(e, base.Sym());
		e = base.Oprev();
	} while (Dest(e) != first);

	Edge close;
	do {
		auto t = e.Oprev();
		if (rightOf(Dest(t), e) && incircle(Org(e), Dest(t), Dest(e), X)) {
			swap(e);
			assert(X == Dest(e));
			e = e.Oprev();
		}
		else if (Org(e) == first)
			break;
		else {
			close = e;
			e = e.Onext().Lprev();
		}
	} while (true);
	close = close.Onext().Sym();
	return{close,true};
}

void insertSiteSequence(Subdivision & s, std::vector<Point> seq)
{
	if (seq.empty())
		return;

	Edge close; bool succ;
	std::tie(close,succ) = insertSite(s, seq[0]);
	if (!succ)
		return;

	for (int i = 1; i < seq.size(); ++i) {
		std::tie(close,succ) = insertSite(s, seq[i], close);
	}
}
