#ifndef _FRAGMENT_GRAPH_NODE_GUARD
#define _FRAGMENT_GRAPH_NODE_GUARD 1


#include <vector>
#include <iostream>
#include <string>


#include "FragmentSubNode.h"


class Molecule;


class FragmentGraphNode
{
  public:
    FragmentGraphNode();
    FragmentGraphNode(const Molecule* const mol, unsigned int id);

    FragmentGraphNode* copy() const;

    // The degree of this node (cardinality of the connections)
    unsigned int degree();

    const Molecule* getMolecule() const { return theMolecule; } 
    FragmentSubNode* getSubNode(unsigned int id) const;
    unsigned int getNodeID() const { return graphID; }

    bool IsIsomorphicTo(FragmentGraphNode* that) const;

    std::string toString() const;
    friend std::ostream& operator<< (std::ostream& os, const FragmentGraphNode& node);

  private:
    const Molecule* theMolecule;
    unsigned int graphID;
    int theDegree;
    std::vector<FragmentSubNode*> subnodes;
};

#endif
