
#include <vector>
#include <iostream> 


#include "RigidFragmentSubNode.h"
#include "LinkerFragmentSubNode.h"
#include "FragmentGraphNode.h"
#include "Molecule.h"
#include "Utilities.h"


FragmentGraphNode::FragmentGraphNode() : theMolecule(),
                                         graphID(-1),
                                         theDegree(-1)
{
}

// ***********************************************************************
//
// Construction of a node for a linker / rigid
//
FragmentGraphNode::FragmentGraphNode(const Molecule* const mol, unsigned int id)
                                   : theMolecule(mol),
                                     graphID(id),
                                     theDegree(-1)
{
    subnodes.clear();

    // Create the subnodes for each possible connection
    std::vector<unsigned int> connectionIDs;
    mol->getConnectionIDs(connectionIDs);

    foreach_uints(conn_id_it, connectionIDs)
    {
        // If there is a connection, create a subnode (using the id).
        if (*conn_id_it != 0)
        {
            FragmentSubNode* subnode = 0;

            if (mol->IsLinker())
            {
                subnode = new LinkerFragmentSubNode(*conn_id_it, this);
            }
            else if (mol->IsRigid())
            {
                subnode = new RigidFragmentSubNode(*conn_id_it, this);
            }
            else
            {
                throw "Expected linker or rigid; not a complex molecule.";
            }

            subnodes.push_back(subnode);
        }
    }
}


// ***********************************************************************
//
// For copying (as part of copying of the entire graph).
//
FragmentGraphNode* FragmentGraphNode::copy() const
{
    FragmentGraphNode* theCopy = new FragmentGraphNode();

    theCopy->theMolecule = this->theMolecule;
    theCopy->graphID = this->graphID;
    
    foreach_subnodes(s_it, this->subnodes)
    {
        theCopy->subnodes.push_back((*s_it)->copy());
    }

    return theCopy;
}


// ***********************************************************************
//
// Acquire the proper subnode object by the given sub-node id.
//
FragmentSubNode* FragmentGraphNode::getSubNode(unsigned int id) const
{
    for (std::vector<FragmentSubNode*>::const_iterator it = this->subnodes.begin();
         it != this->subnodes.end();
         it++)
    {
        if ((*it)->getSubNodeID() == id) return *it;
    }
}

// ***********************************************************************
//
// We calculate the degrees only when the graph has been completely constructed.
//
bool FragmentGraphNode::IsIsomorphicTo(FragmentGraphNode* that) const
{
    if (this->theMolecule->getUniqueIndexID() != that->getMolecule()->getUniqueIndexID())
    {
        return false;
    }

    //
    // Match each subnode in that to this
    //
    std::vector<bool> marked;
    MakeBoolVector(marked, subnodes.size());
    foreach_subnodes(this_it, this->subnodes)
    {
        bool found = false;
        int counter = 0;
        foreach_subnodes(that_it, that->subnodes)
        {
            if (!marked[counter])
            {
                if ((*this_it)->IsIsomorphicTo(*that_it))
                {
                    marked[counter] = true;
                    found = true;
                    break;
                }
            }

            counter++;
        }

        if (!found) return false;
    }

    return !ContainsFalse(marked);
}

// ***********************************************************************
//
// We calculate the degrees only when the graph has been completely constructed.
// 
unsigned int FragmentGraphNode::degree()
{
    if (theDegree >= 0) return (unsigned int)(theDegree); 

    theDegree = 0;
    for (std::vector<FragmentSubNode*>::const_iterator it = subnodes.begin();
         it != subnodes.end();
         it++)
    {
        theDegree += (*it)->degree();
    }

    return theDegree;
}

// ***********************************************************************

std::string FragmentGraphNode::toString() const
{
    std::ostringstream oss;

    oss << "\tFragment Subnode (id=" << this->graphID << "):" << std::endl;

    foreach_subnodes(s_it, this->subnodes)
    {
        oss << **s_it;
    }

    return oss.str();
}

// ***********************************************************************

std::ostream& operator<< (std::ostream& os, const FragmentGraphNode& node)
{
    os << node.toString();

    return os;
}
