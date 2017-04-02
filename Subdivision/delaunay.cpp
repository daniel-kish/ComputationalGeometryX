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

bool rightOf(Vertex v, Edge e)
{
	return orient2d(v->point, Dest(e)->point, Org(e)->point) > 0.0;
}

bool leftOf(Vertex v, Edge e)
{
	return orient2d(v->point, Org(e)->point, Dest(e)->point) > 0.0;
}

double incircle(Vertex a, Vertex b, Vertex c, Vertex d) {
	return incircle((double*)&(a->point), (double*)&(b->point),
		(double*)&(c->point), (double*)&(d->point));
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
			while (incircle(Dest(basel), Org(basel), Dest(lcand), Dest(lcand.Onext())) > 0.0)
			{
				lcand = lcand.Onext();
				L.deleteEdge(lcand.Oprev());
			}
		}
		auto rcand = basel.Oprev();
		if (valid(rcand)) {
			while (incircle(Dest(basel), Org(basel), Dest(rcand), Dest(rcand.Oprev())) > 0.0)
			{
				rcand = rcand.Oprev();
				L.deleteEdge(rcand.Onext());
			}
		}
		if (!valid(lcand) && !valid(rcand)) break; // we've reached the top

		if (!valid(lcand) ||
			valid(rcand) && incircle(Dest(lcand), Org(lcand), Org(rcand), Dest(rcand)) > 0.0)
			basel = L.connect(rcand, basel.Sym());
		else
			basel = L.connect(basel.Sym(), lcand.Sym());
	} while (true);
	return std::make_tuple(std::move(L), ldo, rdo);
}