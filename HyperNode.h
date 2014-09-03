#ifndef _HYPER_NODE_H
#define _HYPER_NODE_H 1

#include <vector>
#include <iostream>
#include <string>


#include "HyperEdge.h"


template<class T, class A>
class HyperNode
{
  public:
    T data;
    int id;
    std::vector<HyperEdge<A> > edges;

    HyperNode(const T& d, int i)
    {
        data = d;
        id = i;
    }

    std::string toString() const;

    void AddEdge(const HyperEdge<A>& edge) { edges.push_back(edge); }
};

template<class T, class A>
std::string HyperNode<T, A>::toString() const
{
    std::ostringstream oss;

    oss << data.toString() << + "\t\t\t\t= { ";

    oss << "(" << id <<") Edges = { ";
    for (int e = 0; e < edges.size(); e++)
    {
        oss << edges[e].toString();
        if (e+1 < edges.size()) oss << ", ";
    }
    oss << " }" << std::endl;

    return oss.str();
}


template<class T, class A>
std::ostream& operator<< (std::ostream& os, HyperNode<T, A>& node)
{
    os << node.toString();

    return os;
}
#endif
