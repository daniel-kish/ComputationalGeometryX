#pragma once
#include <list>

template <class T>
class QuadEdgeList {
public:
	struct QuadEdge;
private:
	using List = std::list<QuadEdge>;
	using QuadEdgeRef = typename List::iterator;

public:
	class EdgeRef {
		QuadEdgeRef ptr;
		int n;

		friend QuadEdgeList;
		friend void splice(EdgeRef a, EdgeRef b) {
			auto alpha = a.Onext().Rot();
			auto beta = b.Onext().Rot();

			std::swap(a.ptr->edges[a.n].next, b.ptr->edges[b.n].next);
			std::swap(alpha.ptr->edges[alpha.n].next, beta.ptr->edges[beta.n].next);
		}
	public:
		EdgeRef(QuadEdgeRef, int);
		EdgeRef();
		T& data();
		EdgeRef Onext();
		EdgeRef Rot();
		EdgeRef InvRot();
		EdgeRef Sym();
		EdgeRef Lnext();
		EdgeRef Oprev();
		EdgeRef Rprev();
		bool operator==(EdgeRef const& other);
		bool operator!=(EdgeRef const& other);
	};

	struct QuadEdge { 
		struct {
			EdgeRef next;
			T data;
		} edges[4];
	};

private:
	List quadEdges;
public:
	EdgeRef makeEdge();
	void deleteEdge(EdgeRef e);
	std::size_t size();

	typename List::const_iterator begin() {
		return quadEdges.cbegin();
	}
	typename List::const_iterator end() {
		return quadEdges.cend();
	}
};

template <class T>
typename QuadEdgeList<T>::EdgeRef QuadEdgeList<T>::makeEdge()
{
	quadEdges.push_back(QuadEdge{});
	QuadEdgeRef qref = std::prev(quadEdges.end());
	qref->edges[0] = {{qref,0},T{}};
	qref->edges[1] = {{qref,3},T{}};
	qref->edges[2] = {{qref,2},T{}};
	qref->edges[3] = {{qref,1},T{}};
	return EdgeRef(qref, 0);
}

template <class T>
void QuadEdgeList<T>::deleteEdge(typename QuadEdgeList<T>::EdgeRef e)
{
	splice(e, e.Oprev());
	splice(e.Sym(), e.Sym().Oprev());
	quadEdges.erase(e.ptr);
}

template <class T>
std::size_t QuadEdgeList<T>::size()
{
	return quadEdges.size();
}

template <class T>
QuadEdgeList<T>::EdgeRef::EdgeRef(QuadEdgeRef qref, int n_) : ptr{qref}, n{n_}
{}

template <class T>
QuadEdgeList<T>::EdgeRef::EdgeRef()
{}

template <class T>
T& QuadEdgeList<T>::EdgeRef::data()
{
	return ptr->edges[n].data;
}


// EdgeRef ifc
template <class T>
typename QuadEdgeList<T>::EdgeRef QuadEdgeList<T>::EdgeRef::Onext() {
	return ptr->edges[n].next;
}

template <class T>
typename QuadEdgeList<T>::EdgeRef QuadEdgeList<T>::EdgeRef::Rot() {
	EdgeRef copy = *this;
	copy.n = (copy.n + 1) % 4;
	return copy;
}

template <class T>
typename QuadEdgeList<T>::EdgeRef QuadEdgeList<T>::EdgeRef::InvRot() {
	EdgeRef copy = *this;
	copy.n = (copy.n + 3) % 4;
	return copy;
}

template <class T>
typename QuadEdgeList<T>::EdgeRef QuadEdgeList<T>::EdgeRef::Sym() {
	return (*this).Rot().Rot();
}

template <class T>
typename QuadEdgeList<T>::EdgeRef QuadEdgeList<T>::EdgeRef::Lnext() {
	return this->InvRot().Onext().Rot();
}

template <class T>
typename QuadEdgeList<T>::EdgeRef QuadEdgeList<T>::EdgeRef::Oprev() {
	return this->Rot().Onext().Rot();
}

template <class T>
typename QuadEdgeList<T>::EdgeRef QuadEdgeList<T>::EdgeRef::Rprev() {
	return this->Sym().Onext();
}

template <class T>
bool QuadEdgeList<T>::EdgeRef::operator==(typename QuadEdgeList<T>::EdgeRef const& other) {
	return this->ptr == other.ptr && this->n == other.n;
}

template <class T>
bool QuadEdgeList<T>::EdgeRef::operator!=(typename QuadEdgeList<T>::EdgeRef const& other) {
	return !(*this == other);
}