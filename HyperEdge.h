#ifndef _HYPEREDGE_GUARD
#define _HYPEREDGE_GUARD 1

#include <vector>
#include <algorithm>
#include <iostream>
#include <string>
#include <sstream>
#include <string>


#include "EdgeAnnotation.h"


template<class A>
class HyperEdge
{
  public:
  
    int targetNode;
    std::vector<int> sourceNodes;
    A annotation;
    
    HyperEdge(const std::vector<int>& src, int target, const A& annot);
    ~HyperEdge() {}

    bool DefinesEdge(const std::vector<int>& thatAnte, int thatConseq) const;
    
    bool operator==(const HyperEdge<A>& that) const
    {
        return this->DefinesEdge(that.sourceNodes, that.targetNode);
    }

    std::string toString() const;

    template<class AS>  
    friend std::ostream& operator<< (std::ostream& os, HyperEdge<AS>& edge);
};

template<class A>
HyperEdge<A>::HyperEdge(const std::vector<int>& src, int target, const A& annot)
{
   sourceNodes = src;
   targetNode = target;
   annotation = annot;
}

/*************************************************************************/

template<class A>
bool HyperEdge<A>::DefinesEdge(const std::vector<int>& thatAnte, int thatConseq) const
{
    if (this->targetNode != thatConseq) return false;

    if (this->sourceNodes.size() != thatAnte.size()) return false;
	
    // Only need to check containment in one direction.
    for (int a = 0; a < thatAnte.size(); a++)
    {
        if (std::find(sourceNodes.begin(), sourceNodes.end(), thatAnte[a]) == sourceNodes.end())
        {
            return false;
        }
    }

	return true;
}

template<class A>
std::string HyperEdge<A>::toString() const
{
    std::ostringstream oss;

    oss << " { ";

    for (int n = 0; n < sourceNodes.size(); n++)
    {
        oss << sourceNodes[n];
        if (n+1 < sourceNodes.size()) oss << ", ";
    }
    oss << " } -> " << targetNode;

    return oss.str();
}

/*************************************************************************/

template<class A>
std::ostream& operator<< (std::ostream& os, HyperEdge<A>& edge)
{
    os << edge.toString();

    return os;
}

#endif
