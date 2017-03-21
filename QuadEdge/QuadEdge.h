#pragma once
#include <list>
#include <array>

template <class T>
class QuadEdgeList 
{
public:
	class QuadEdge;
	class EdgeRef;
	using List = std::list<QuadEdge>;
	using QuadEdgeRef = typename List::iterator;

	EdgeRef makeEdge();
	void deleteEdge(EdgeRef);
	std::size_t size() const;
	QuadEdgeRef begin();
	QuadEdgeRef end();
private:
	List quadEdges;
};

template <class T>
class QuadEdgeList<T>::EdgeRef
{
	typename QuadEdgeList<T>::QuadEdgeRef qref;
	int n;
public:
	EdgeRef(typename QuadEdgeList<T>::QuadEdgeRef, int);
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
	
	friend QuadEdgeList;
	friend void splice(typename QuadEdgeList<T>::EdgeRef a, typename QuadEdgeList<T>::EdgeRef b)
	{
		auto alpha = a.Onext().Rot();
		auto beta = b.Onext().Rot();
		std::swap(a.qref->recs[a.n].next, b.qref->recs[b.n].next);
		std::swap(alpha.qref->recs[alpha.n].next, beta.qref->recs[beta.n].next);
	}
};

template <class T>
class QuadEdgeList<T>::QuadEdge
{
public:
	struct Record {
		typename QuadEdgeList<T>::EdgeRef next;
		T data;
	};
	friend QuadEdgeList<T>;
	friend QuadEdgeList<T>::EdgeRef;
	friend void splice(typename QuadEdgeList<T>::EdgeRef a, typename QuadEdgeList<T>::EdgeRef b);
private:
	std::array<Record, 4> recs;
public:
	QuadEdge();
	std::array<Record, 4> const& records() const;
};


// QuadEdge
template <class T>
QuadEdgeList<T>::QuadEdge::QuadEdge()
{
}

template <class T>
std::array<typename QuadEdgeList<T>::QuadEdge::Record, 4> const&
QuadEdgeList<T>::QuadEdge::records() const
{
	return recs;
}

// EdgeRef 
template <class T>
QuadEdgeList<T>::EdgeRef::EdgeRef(typename QuadEdgeList<T>::QuadEdgeRef qref_, int n_)
	: qref{qref_}, n{n_}
{}

template <class T>
QuadEdgeList<T>::EdgeRef::EdgeRef()
{}

template <class T>
T& QuadEdgeList<T>::EdgeRef::data()
{
	return qref->recs[n].data;
}

template <class T>
typename QuadEdgeList<T>::EdgeRef QuadEdgeList<T>::EdgeRef::Onext()
{
	return qref->recs[n].next;
}

template <class T>
typename QuadEdgeList<T>::EdgeRef QuadEdgeList<T>::EdgeRef::Rot()
{
	typename QuadEdgeList<T>::EdgeRef copy = *this;
	copy.n = (copy.n + 1) % 4;
	return copy;
}

template <class T>
typename QuadEdgeList<T>::EdgeRef QuadEdgeList<T>::EdgeRef::InvRot()
{
	typename QuadEdgeList<T>::EdgeRef copy = *this;
	copy.n = (copy.n + 3) % 4;
	return copy;
}

template <class T>
typename QuadEdgeList<T>::EdgeRef QuadEdgeList<T>::EdgeRef::Sym()
{
	typename QuadEdgeList<T>::EdgeRef copy = *this;
	copy.n = (copy.n + 2) % 4;
	return copy;
}

template <class T>
typename QuadEdgeList<T>::EdgeRef QuadEdgeList<T>::EdgeRef::Lnext()
{
	return this->InvRot().Onext().Rot();
}

template <class T>
typename QuadEdgeList<T>::EdgeRef QuadEdgeList<T>::EdgeRef::Oprev()
{
	return this->Rot().Onext().Rot();
}

template <class T>
typename QuadEdgeList<T>::EdgeRef QuadEdgeList<T>::EdgeRef::Rprev()
{
	return this->Sym().Onext();
}

template <class T>
bool QuadEdgeList<T>::EdgeRef::operator==(typename QuadEdgeList<T>::EdgeRef const& other) {
	return this->qref == other.qref && this->n == other.n;
}

template <class T>
bool QuadEdgeList<T>::EdgeRef::operator!=(typename QuadEdgeList<T>::EdgeRef const& other) {
	return !(*this == other);
}

// QuadEdgeList
template <class T>
typename QuadEdgeList<T>::EdgeRef QuadEdgeList<T>::makeEdge()
{
	quadEdges.push_back(QuadEdge());
	QuadEdgeRef qref = std::prev(quadEdges.end());
	qref->recs[0] = {EdgeRef(qref,0),T{}};
	qref->recs[1] = {EdgeRef(qref,3),T{}};
	qref->recs[2] = {EdgeRef(qref,2),T{}};
	qref->recs[3] = {EdgeRef(qref,1),T{}};
	return EdgeRef(qref, 0);
}

template<class T>
void QuadEdgeList<T>::deleteEdge(typename QuadEdgeList<T>::EdgeRef e)
{
	splice(e, e.Oprev());
	splice(e.Sym(), e.Sym().Oprev());
	quadEdges.erase(e.qref);
}

template<class T>
std::size_t QuadEdgeList<T>::size() const
{
	return quadEdges.size();
}

template<class T>
typename QuadEdgeList<T>::QuadEdgeRef QuadEdgeList<T>::begin()
{
	return quadEdges.begin();
}

template<class T>
typename QuadEdgeList<T>::QuadEdgeRef QuadEdgeList<T>::end()
{
	return quadEdges.end();
}