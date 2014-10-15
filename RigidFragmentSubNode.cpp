#include <sstream>


#include "Molecule.h"
#include "FragmentGraphNode.h"
#include "RigidFragmentSubNode.h"
#include "Utilities.h"


// ************************************************************************************

RigidFragmentSubNode::RigidFragmentSubNode() : connection(0), FragmentSubNode(-1, 0)
{
}

// ************************************************************************************

RigidFragmentSubNode::RigidFragmentSubNode(unsigned int id, FragmentGraphNode* parent) :
                                           connection(0), FragmentSubNode(id, parent)
{
}

// ************************************************************************************

FragmentSubNode* RigidFragmentSubNode::copy() const
{
    RigidFragmentSubNode* theCopy = new RigidFragmentSubNode();

    theCopy->uniqueSubnodeID = this->uniqueSubnodeID;
    theCopy->parentFragmentNode = this->parentFragmentNode;

    theCopy->connection = this->connection;

    return theCopy;
}


// ************************************************************************************

void RigidFragmentSubNode::addConnection(FragmentSubNode* connector)
{
    if (connection != 0) throw "Attempt to add a connection when one exists";

    connection = connector;
}

// ************************************************************************************

bool RigidFragmentSubNode::IsIsomorphicTo(FragmentSubNode* that)
{
    if (this->uniqueSubnodeID != that->getSubNodeID()) return false;

    if (! instanceOf<RigidFragmentSubNode>(that)) return false;

    RigidFragmentSubNode* thatRigid = dynamic_cast<RigidFragmentSubNode*>(that);

    // If both connections are non-existent, then these nodes match.
    if (this->connection == 0 && thatRigid->connection == 0) return true;

    // If either connections are null (invalid), they are not isomorphic
    if (this->connection == 0 || thatRigid->connection == 0) return false;

    // Otherwise, both connections exist; check the target subnode.
    return this->connection->getSubNodeID() == thatRigid->connection->getSubNodeID();
}

// ************************************************************************************

std::string RigidFragmentSubNode::toString() const
{
    std::ostringstream oss;

    if (uniqueSubnodeID != -1)
    {
        if (connection)
        {
            oss << "\t\t" << uniqueSubnodeID << " <--> ("
                << connection->getParentNode()->getNodeID() << ", "
                << connection->getSubNodeID()
                << ")";
        }
        else oss << "\t\t" << uniqueSubnodeID << " <--> empty";

        oss << std::endl;
    }

    return oss.str();
}

// ************************************************************************************

