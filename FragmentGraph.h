#ifndef _FRAGMENT_GRAPH_GUARD
#define _FRAGMENT_GRAPH_GUARD 1


#include <vector>
#include <iostream>
#include <string>
#include <map>


#include "FragmentGraphNode.h"


class Atom;


class FragmentGraph
{
  public:
    FragmentGraph();
    FragmentGraph* copy() const;

    // Returns the index of the new node
    unsigned int AddInitialNode(const Molecule* const mol);
    std::pair<unsigned int, unsigned int> AddEdgeAndNode(unsigned int fromConnID,
                                                         std::pair<unsigned int,
                                                         unsigned int> fromNodeIndex,
                                                         const Atom& to,
                                                         const Molecule& thatMol);

    bool IsIsomorphicTo(FragmentGraph* that) const;

    std::string toString() const;
    friend std::ostream& operator<< (std::ostream& os, const FragmentGraph& fg);

  private:
    // We order the nodes by the particular fragment used;
    // one list for each fragment 
    std::vector<FragmentGraphNode*>* orderedNodes;
    unsigned int numFragments;

    void printMolecules() const;
};
	
#endif
