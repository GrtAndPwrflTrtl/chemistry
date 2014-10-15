
#include <string>
#include <sstream>


#include "Molecule.h"
#include "LinkerFragmentSubNode.h"
#include "FragmentGraphNode.h"
#include "Utilities.h"


// ************************************************************************************

LinkerFragmentSubNode::LinkerFragmentSubNode() : FragmentSubNode(-1, 0)
{
    connections.clear();
}

// ************************************************************************************

LinkerFragmentSubNode::LinkerFragmentSubNode(unsigned int id, FragmentGraphNode* parent)
                                            : FragmentSubNode(id, parent)
{
   connections.clear();
}

// ************************************************************************************

FragmentSubNode* LinkerFragmentSubNode::copy() const
{
    LinkerFragmentSubNode* theCopy = new LinkerFragmentSubNode();

    theCopy->uniqueSubnodeID = this->uniqueSubnodeID;
    theCopy->parentFragmentNode = this->parentFragmentNode;

    for (std::vector<FragmentSubNode*>::const_iterator it = this->connections.begin();
         it != this->connections.end();
         it++)
    {
        theCopy->connections.push_back(*it);
    }

    return theCopy;
}

// ************************************************************************************

void LinkerFragmentSubNode::addConnection(FragmentSubNode* connector)
{
/*
    if (parentFragmentNode->getMaxConnect() == connections.size())
    {
        throw "Attempt to add a connection when the maximum has been reached.";
    }
*/

    connections.push_back(connector);
}

// ************************************************************************************

bool LinkerFragmentSubNode::IsIsomorphicTo(FragmentSubNode* that)
{
    if (this->uniqueSubnodeID != that->getSubNodeID()) return false;

    if (! instanceOf<LinkerFragmentSubNode>(that)) return false;

    LinkerFragmentSubNode* thatL = dynamic_cast<LinkerFragmentSubNode*>(that);

    if (this->connections.size() != thatL->connections.size()) return false;

    //
    // Match each subnode in that to this
    //
    std::vector<bool> marked;
    MakeBoolVector(marked, this->connections.size());
    foreach_subnodes(this_it, this->connections)
    {
        bool found = false;
        int counter = 0;
        foreach_subnodes(that_it, thatL->connections)
        {
            if (!marked[counter])
            {
                if ((*this_it)->getSubNodeID() == (*that_it)->getSubNodeID())
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

// ************************************************************************************

std::string LinkerFragmentSubNode::toString() const
{
    std::ostringstream oss;

//    oss << uniqueSubnodeID << ": ";
//    oss << parentFragmentNode->getMolecule()->getFragmentId() << " ";

    if (this->connections.empty()) oss << "\t\t" << uniqueSubnodeID << " <--> empty"; 
    else
    {
        foreach_subnodes(n_it, this->connections)
        {
            if (uniqueSubnodeID != -1)
            {
                oss << "\t\t" << uniqueSubnodeID << " <--> ("
                    << (*n_it)->getParentNode()->getNodeID() << ", "
                    << (*n_it)->getSubNodeID()
                    << ")\t";
            }
        }
    }

    oss << std::endl;

    return oss.str();
}

// ************************************************************************************

