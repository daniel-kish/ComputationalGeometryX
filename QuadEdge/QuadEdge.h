#pragma once
#include <list>
#include <array>

template <class T>
class QuadEdgeList {
public:
	class QuadEdge;
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

			a.ptr->recs[a.n];
			std::swap(a.ptr->recs[a.n].next, b.ptr->recs[b.n].next);
			std::swap(alpha.ptr->recs[alpha.n].next, beta.ptr->recs[beta.n].next);
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

	class QuadEdge { 
	public:
		struct Record {
			EdgeRef next;
			T data;
		};
		friend void splice(EdgeRef a, EdgeRef b);
	private:
		std::array<Record,4> recs[4];
	public:
		QuadEdge();
		std::array<Record, 4> const& records() const {
			return records;
		}
	};
	
private:
	List quadEdges;
public:
	EdgeRef makeEdge();
	void deleteEdge(EdgeRef e);
	std::size_t size();

	typename List::iterator begin() {
		return quadEdges.begin();
	}
	typename List::iterator end() {
		return quadEdges.end();
	}
};

template <class T>
typename QuadEdgeList<T>::EdgeRef QuadEdgeList<T>::makeEdge()
{
	quadEdges.push_back(QuadEdge{});
	QuadEdgeRef qref = std::prev(quadEdges.end());
	qref->recs[0] = {{qref,0},T{}};
	qref->recs[1] = {{qref,3},T{}};
	qref->recs[2] = {{qref,2},T{}};
	qref->recs[3] = {{qref,1},T{}};
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

template<class T>
QuadEdgeList<T>::QuadEdge::QuadEdge()
{}

template <class T>
QuadEdgeList<T>::EdgeRef::EdgeRef(QuadEdgeRef qref, int n_) : ptr{qref}, n{n_}
{}

template <class T>
QuadEdgeList<T>::EdgeRef::EdgeRef()
{}

template <class T>
T& QuadEdgeList<T>::EdgeRef::data()
{
	return ptr->recs[n].data;
}


// EdgeRef ifc
template <class T>
typename QuadEdgeList<T>::EdgeRef QuadEdgeList<T>::EdgeRef::Onext() {
	return ptr->recs[n].next;
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