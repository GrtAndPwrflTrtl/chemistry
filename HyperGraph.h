#ifndef _HYPER_GRAPH_GUARD
#define _HYPER_GRAPH_GUARD 1

#include <vector>
#include <exception>
#include <algorithm>
#include <utility>
#include <sstream>
#include <string>

#include "PebblerHyperGraph.h"
#include "HyperNode.h"
#include "HyperEdge.h"
#include "Utilities.h"
#include "Constants.h"

//
// The goal is three-fold in this class.
//   (1) Provides functionality to create a hypergraph: adding nodes and edges
//   (2) Convert all clauses to an integer hypergraph representation.
//   (3) Provide functionality to explore the hypergraph.
//
template<class T, class A>
class HyperGraph
{
  public:

    // The main graph data structure
    std::vector<HyperNode<T, A> > vertices;

    HyperGraph(unsigned int numBuckets)
    {
        // Initialize the database of nodes that have the same size;
        // these are essentially buckets to speed searching
        buckets = new std::vector<int>[numBuckets]; 
    }
    ~HyperGraph() { }
    int size() { return vertices.size(); }

    //
    // Integer-based representation of the main hypergraph
    //
    //PebblerHyperGraph<T> GetPebblerHyperGraph();
    int GetNodeIndex(const T& inputData);
    T& GetNode(int id);
    bool HasNode(const T& inputData);
    T GetNode(const T& inputData);
    bool AddNode(const T& inputData);
    // Check if the graph contains an edge defined by a many to one clause mapping
    bool HasEdge(const std::vector<T>& antecedent, const T& consequent);
    void AddEdge(const std::vector<T>& antecedent, const T& consequent, const A& annotation);
    
    template<class TS, class AS>
    friend std::ostream& operator<< (std::ostream& os, HyperGraph<TS, AS>& graph);
    std::string toString() const;

    PebblerHyperGraph<T, A> GetPebblerHyperGraph() const;
    std::vector<T> CollectData() const;


  private:
    // Check if the graph contains this specific grounded clause
    int ConvertToLocalIntegerIndex(const T& inputData);
    // Is this edge in the graph (using local, integer-based information)
    bool HasLocalEdge(const std::vector<int>& antecedent, int consequent);
    // Convert information to local, integer-based representation
    std::pair<std::vector<int>, int> ConvertToLocal(const std::vector<T>& antecedent, const T& consequent);

    // A 'database' of nodes based on the size; the class T must implement a method called size
    // These are buckets to expedite searching for equivalent molecules.
    std::vector<int>* buckets;
};


        
//
// Create the Integer-based representation of the main hypergraph
//
template<class T, class A>
PebblerHyperGraph<T, A> HyperGraph<T, A>::GetPebblerHyperGraph() const
{
    //
    // Strictly create the nodes
    //
    std::vector<PebblerHyperNode<T, A> > pebblerNodes;
    for (int v = 0; v < vertices.size(); v++)
    {
        pebblerNodes.push_back(PebblerHyperNode<T, A>(vertices[v].data, vertices[v].id));
    }

    //
    // Create all hyperedges
    //
    for (int v = 0; v < vertices.size(); v++)
    {
        HyperNode<T, A> node = vertices[v];

        for (int e = 0; e < node.edges.size(); e++)
        {
            std::sort(node.edges[e].sourceNodes.begin(), node.edges[e].sourceNodes.end());

            pebblerNodes[v].edges.push_back(PebblerHyperEdge<A>(node.edges[e].sourceNodes,
                                            node.edges[e].targetNode, node.edges[e].annotation));
        }
    }

    return PebblerHyperGraph<T, A>(pebblerNodes);
}


// Check if the graph contains this specific grounded clause
//
template<class T, class A>
int HyperGraph<T, A>::ConvertToLocalIntegerIndex(const T& inputData)
{
    int sz = inputData.size();

    for (std::vector<int>::const_iterator it = buckets[sz].begin();
         it != buckets[sz].end();
         it++)
    {
        if (vertices[*it].data == inputData) return *it;
    }

    return -1;
}

//
// Return the index of the given node
//
template<class T, class A>
int HyperGraph<T, A>::GetNodeIndex(const T& inputData)
{
    return ConvertToLocalIntegerIndex(inputData);
}

//
// Return the stored node in the graph
//
template<class T, class A>
T& HyperGraph<T, A>::GetNode(int id)
{
    if (id < 0 || id > vertices.size())
    {
        throw MakeString("Unexpected id in hypergraph node access: ", id);
    }

    return vertices[id].data;
}

//
// Check if the graph contains this specific grounded clause
//
template<class T, class A>
bool HyperGraph<T, A>::HasNode(const T& inputData)
{
    // Only check the exact set of nodes that have the same 'size'.
    // The particular bucket contains indices.
    int sz = inputData.size();

    for (std::vector<int>::const_iterator it = buckets[sz].begin(); it != buckets[sz].end(); it++)
    {
        if (vertices[*it].data == inputData) return true;
    }

    return false;
}

//
// Check if the graph contains this specific grounded clause
//
template<class T, class A>
T HyperGraph<T, A>::GetNode(const T& inputData)
{
    // Only check the exact set of nodes that have the same 'size'.
    // The particular bucket contains indices.
    int sz = inputData.size();

    for (std::vector<int>::const_iterator it = buckets[sz].begin(); it != buckets[sz].end(); it++)
    {
        if (vertices[*it].data == inputData) return vertices[*it].data;
    }

    throw null;
}

//
// Check if the graph contains this specific grounded clause
//
template<class T, class A>
bool HyperGraph<T, A>::AddNode(const T& inputData)
{
    if (HasNode(inputData)) return false;

                                      // <data,   id>
    vertices.push_back(HyperNode<T, A>(inputData, vertices.size()));

    // Place the index of the newly added node in the proper bucket.
    buckets[inputData.size()].push_back(vertices.size() - 1);

    return true;
}

//
// Is this edge in the graph (using local, integer-based information)
//
template<class T, class A>
bool HyperGraph<T, A>::HasLocalEdge(const std::vector<int>& antecedent, int consequent)
{
    for (int v = 0; v < vertices.size(); v++) 
    {
        for (int e = 0; e < vertices[v].edges.size(); e++) 
        {
            if (vertices[v].edges[e].DefinesEdge(antecedent, consequent)) return true;
        }
    }

    return false;
}

//
// Check if the graph contains an edge defined by a many to one clause mapping
//
template<class T, class A>
bool HyperGraph<T, A>::HasEdge(const std::vector<T>& antecedent, const T& consequent)
{
    std::pair<std::vector<int>, int> local = ConvertToLocal(antecedent, consequent);

    return HasLocalEdge(local.first, local.second);
}

//
// Convert information to local, integer-based representation
//
template<class T, class A>
std::pair<std::vector<int>, int> HyperGraph<T, A>::ConvertToLocal(const std::vector<T>& antecedent,
                                                                  const T& consequent)
{
    std::vector<int> localAnte;

    for (int a = 0; a < antecedent.size(); a++)
    {
        // This speeds things up since we don't have to look up each node
        // (avoid isomorphism check).
        int index = antecedent[a].getUniqueIndexID();
/*
        int index = ConvertToLocalIntegerIndex(antecedent[a]);
*/
        if (index == -1)
        {
            std::string err = MakeString("Source node not found as a hypergraph node: \n",
                                         antecedent[a].toString());
            std::cerr <<  err << std::endl;
            throw err;
        }

        localAnte.push_back(index);
    }

    int localConsequent = ConvertToLocalIntegerIndex(consequent);

    if (localConsequent == -1)
    {
        throw MakeString("Target value referenced not found as a hypergraph node", consequent.toString());
    }

    return std::make_pair(localAnte, localConsequent);
}

//
// Adding an edge to the graph
//
template<class T, class A>
void HyperGraph<T, A>::AddEdge(const std::vector<T>& antecedent, const T& consequent, const A& annotation)
{
    // Add a local representaiton of this edge to each node in which it is applicable
    if (HasEdge(antecedent, consequent)) return;

    std::pair<std::vector<int>, int> local = ConvertToLocal(antecedent, consequent);

    HyperEdge<A> edge(local.first, local.second, annotation);

//System.Diagnostics.Debug.WriteLine("Adding edge: " + edge.ToString());

    for (int s = 0; s < local.first.size(); s++) 
    {
        vertices[s].AddEdge(edge);
    }
}
template<class T, class A>
std::string HyperGraph<T, A>::toString() const
{
    std::ostringstream oss;

    for (int v = 0; v < vertices.size(); v++)
    {
        oss << v << ": " << vertices[v].toString() << std::endl;
    }     


    return oss.str();
}

template<class T, class A>
std::ostream& operator<< (std::ostream& os, HyperGraph<T, A>& graph)
{
    os << graph.toString();
    
    return os;
}

template<class T, class A>
std::vector<T> HyperGraph<T, A>::CollectData() const
{
    std::vector<T> theData;

    for (int v = 0; v < vertices.size(); v++)
    {
        theData.push_back(vertices[v].data);
    }

    return theData;
}


#endif
