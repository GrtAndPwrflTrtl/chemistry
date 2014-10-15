#include <vector>


#include "FragmentGraph.h"
#include "Molecule.h"
#include "Utilities.h"
#include "Atom.h"


// ***********************************************************************

FragmentGraph::FragmentGraph() : numFragments(0)
{
    orderedNodes = new std::vector<FragmentGraphNode*>[Molecule::NUM_UNIQUE_FRAGMENTS];
}

// ***********************************************************************

FragmentGraph* FragmentGraph::copy() const
{
    FragmentGraph* newGraph = new FragmentGraph();

    for (int f = 0; f < Molecule::NUM_UNIQUE_FRAGMENTS; f++)
    {
        foreach_nodes(n_it, this->orderedNodes[f])
        {
            newGraph->orderedNodes[f].push_back((*n_it)->copy());
        }
    }

    newGraph->numFragments = this->numFragments;

    return newGraph;
}

// ***********************************************************************

unsigned int FragmentGraph::AddInitialNode(const Molecule* const mol)
{
    if (mol->IsComplex()) throw "Cannot construct a fragment graph with non-fragment.";

    // Create the graph node for this molecule.
    FragmentGraphNode* node = new FragmentGraphNode(mol, numFragments++);

    // Add the new node (with sub-nodes) to the graph.
    unsigned int uniqueMolID = mol->getUniqueIndexID();

    unsigned int newIndex = orderedNodes[uniqueMolID].size();

    orderedNodes[uniqueMolID].push_back(node);

    return newIndex;
}

// ***********************************************************************
// Return the index (indices) of this new target molecule so atoms can be updated accordingly.
//
std::pair<unsigned int, unsigned int>
FragmentGraph::AddEdgeAndNode(unsigned int fromConnId,
                              std::pair<unsigned int, unsigned int> fromGraphNodeIndex,
                              const Atom& to,
                              const Molecule& thatMol)
{
/*
    std::cout << "Printing molecules: " << std::endl;
    printMolecules();


std::cerr << "Graph: " << std::endl << *this << std::endl;

std::cerr << "+++++++++++" << std::endl << thatMol << std::endl;



std::cerr << "Adding: " << fromConnId << " ("
          << fromGraphNodeIndex.first << ", "
          << fromGraphNodeIndex.second << ") " << to.toString() << std::endl;
*/

    // Create a new node for the 'to' node
    FragmentGraphNode* toNode = new FragmentGraphNode(&thatMol, numFragments++);

    // Add the new node (with sub-nodes) to the graph.
    orderedNodes[thatMol.getUniqueIndexID()].push_back(toNode);

    //
    // Attach the 'from' node to the 'to' node via subnodes
    //
    // Acquire the sub-nodes
    FragmentSubNode* fromSubNode = this->orderedNodes[fromGraphNodeIndex.first][fromGraphNodeIndex.second]->getSubNode(fromConnId);
    FragmentSubNode* toSubNode = toNode->getSubNode(to.getConnectionID());

    if (fromSubNode == 0) throw "From subnode not found.";
    if (toSubNode == 0) throw "To subnode not found.";

// std::cout << "from, to: " << fromSubNode << " " << toSubNode << std::endl;

    fromSubNode->addConnection(toSubNode);
    toSubNode->addConnection(fromSubNode);

    // Return the indices of the new 'to' molecule in the graph.
    return std::make_pair(thatMol.getUniqueIndexID(),
                          orderedNodes[thatMol.getUniqueIndexID()].size() - 1);
}

// ***********************************************************************
//
// We check isomorphism between two graphs by a single linear pass over the 
//
//
bool FragmentGraph::IsIsomorphicTo(FragmentGraph* that) const
{
    //
    // We should only arrive at this point in the code by having the same number
    // of fragments (nodes) in each graph
    //
    if (this->numFragments != that->numFragments) return false;

    //
    // We perform a linear pass over the nodes to verify that each molecule graph uses
    // the same number of each fragment  
    //
    for (int f = 0; f < Molecule::NUM_UNIQUE_FRAGMENTS; f++)
    {
        if (this->orderedNodes[f].size() != that->orderedNodes[f].size()) return false;
    }

    //
    // For each specific fragment, we perform isomorphism check
    // This is an n^2 operation
    //
    for (int f = 0; f < Molecule::NUM_UNIQUE_FRAGMENTS; f++)
    {
        std::vector<FragmentGraphNode*> thisList = orderedNodes[f];
        std::vector<FragmentGraphNode*> thatList = that->orderedNodes[f];
       
        if (!thisList.empty() && !thatList.empty())
        {
            std::vector<bool> marked;
            MakeBoolVector(marked, thisList.size());
            foreach_nodes(this_it, thisList)
            {
                bool found = false;
                int counter = 0;
                foreach_nodes(that_it, thatList)
                {
                    if (!marked[counter])
                    {
                        if ((*this_it)->IsIsomorphicTo(*that_it))
                        {
                            found = true;
                            marked[counter] = true;
                            break;
                        }
                    }
                    counter++;
                }
                if (!found) return false;
            }
            if (ContainsFalse(marked)) return false;
        }
    }

    return true;
}

// *****************************************************************************

std::string FragmentGraph::toString() const
{
    std::ostringstream oss;

    oss << "Fingerprint (fragment graph)" << std::endl;

    for (unsigned int f = 0; f < Molecule::NUM_UNIQUE_FRAGMENTS; f++)
    {
        oss << f << ": " << std::endl;

        foreach_nodes(n_it, this->orderedNodes[f])
        {
            oss << **n_it;
        }
    }

    return oss.str();
}


// *****************************************************************************

std::ostream& operator<< (std::ostream& os, const FragmentGraph& fg)
{
    os << fg.toString() << std::endl;

    return os;
}


// ***********************************************************************

void FragmentGraph::printMolecules() const
{
    for (unsigned int f = 0; f < Molecule::NUM_UNIQUE_FRAGMENTS; f++)
    {
        std::cout << f << ": " << std::endl;

        foreach_nodes(n_it, this->orderedNodes[f])
        {
            std::cout << *(*n_it)->getMolecule();
        }
    }
}


// ***********************************************************************

