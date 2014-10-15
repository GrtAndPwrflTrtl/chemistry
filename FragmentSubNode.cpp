#include <string>
#include <sstream>


#include "FragmentGraphNode.h"
#include "Molecule.h"
#include "FragmentSubNode.h"
#include "RigidFragmentSubNode.h"
#include "LinkerFragmentSubNode.h"
#include "Utilities.h"


// **********************************************************************************

FragmentSubNode::FragmentSubNode() : uniqueSubnodeID(-1), parentFragmentNode(0)
{
}

// **********************************************************************************

FragmentSubNode::FragmentSubNode(unsigned int id, FragmentGraphNode* parent)
                               : uniqueSubnodeID(id), parentFragmentNode(parent)
{
}

// **********************************************************************************

std::string FragmentSubNode::toString() const
{
    std::ostringstream oss;
	  
    if (this->parentFragmentNode->getMolecule()->IsLinker())
    { 
        oss << ((LinkerFragmentSubNode*)this)->toString();
    }
    else if (this->parentFragmentNode->getMolecule()->IsRigid())
    {
        oss << ((RigidFragmentSubNode*)this)->toString();
    }
    else throw "Expected a linker or rigid; found something else.";

    return oss.str();
}


// **********************************************************************************

std::ostream& operator<< (std::ostream& os, const FragmentSubNode& node)
{
    os << node.toString();

    return os;    
}
