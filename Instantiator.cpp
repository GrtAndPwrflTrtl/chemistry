#include <vector>
#include <queue>
#include <iostream>
#include <memory>
#include <time.h>
#include <pthread.h>


#include "Molecule.h"
#include "Rigid.h"
#include "Linker.h"


#include "HyperGraph.h"
#include "HyperEdge.h"
#include "EdgeAggregator.h"


#include "Instantiator.h"
#include "Utilities.h"
#include "IdFactory.h"
#include "Constants.h"


Instantiator::Instantiator(std::ostream& out) : ds(out) 
{
    // The hypergraph lock
    pthread_mutex_init(&graph_lock, NULL);

    // The threads and locks for the producer-consumer containers.
    queue_locks = new pthread_mutex_t[HIERARCHICAL_LEVEL_BOUND]; 
    queue_threads = new pthread_t[HIERARCHICAL_LEVEL_BOUND];
    completed_level = new bool[HIERARCHICAL_LEVEL_BOUND];
    level_queues = new std::queue<Molecule*>[HIERARCHICAL_LEVEL_BOUND];

    for (int m = 0; m < HIERARCHICAL_LEVEL_BOUND; m++)
    {
        // Initialize those container locks.
        pthread_mutex_init(&queue_locks[m], NULL);
	
        // Initialize the level-based threads.
        // ?

        // Initialize the fact that we have not computed this level. 
        completed_level[m] = false;

        // Initialize the level queues?
        // level_queues[m].clear();
    }
}



//
// Add the hyperedge to the hypergraph
//
void Instantiator::AddEdge(const std::vector<Molecule>& antecedent,
                           const Molecule& consequent,
                           const EdgeAnnotationT& annotation)
{
    pthread_mutex_lock(&graph_lock);

    graph.AddEdge(antecedent, consequent, annotation);

    pthread_mutex_unlock(&graph_lock);
}

//
// This function consumes from a level L queue and produces molecules in the level (L + 1) queue
//
void Instantiator::ProcessLevel(std::vector<Molecule*>& baseMols,
                                std::queue<Molecule*>& inSet,
                                std::queue<Molecule*>& outSet,
                                pthread_mutex_t& in_lock,
                                pthread_mutex_t& out_lock,
                                bool* previousLevelComplete,
                                bool* thisLevelComplete)
{
    struct timespec sleepTime;
    struct timespec remTime;     // Remaining time
    sleepTime.tv_sec = 0;
    sleepTime.tv_nsec = 100000000L; // 0.1 seconds

    //
    // Keep consuming molecules as long as the previous level is incomplete or this
    // level queue contains molecules to process.
    //
    while (!(*previousLevelComplete) || !inSet.empty())
    {
        // Nothing to process, currently, but the level is incomplete.
        if (inSet.empty())
        {
            nanosleep(&sleepTime, &remTime);
        }

        //
        // Process a molecule in the queue
        //
        else
        {
            //
            // Acquire a molecule to process.
            //
            pthread_mutex_lock(&in_lock);
            Molecule* molToProcess = inSet.front();
            inSet.pop();
            pthread_mutex_lock(&in_lock);

            //
            // Process the molecule by composing it with all the base molecules.
            //
            for (std::vector<Molecule*>::iterator baseMol = baseMols.begin();
                 baseMol != baseMols.end(); baseMol++)
            {
                std::vector<EdgeAggregator*>* newEdges = molToProcess->Compose(**baseMol);

                HandleNewMolecules(outSet, &out_lock, *newEdges);

                for (int i = 0; i < newEdges->size(); i++)
                {
                    delete (*newEdges)[i];
                }

                delete newEdges;
            }
        }
    }
    
    // Indicate this level is complete.
    *thisLevelComplete = true;
}

//
// Threaded construction of the hypergraph using a hierarchical list of threads and containers.
// We first construct the base case of 2-Molecules.
// Then, we inductively start constructing 3-Molecules, 4-Molecules, etc.
//
HyperGraph<Molecule, EdgeAnnotationT> Instantiator::ThreadedInstantiate(std::vector<Linker*>& linkers,
                                                                        std::vector<Rigid*>& rigids)
{
    // The set of all base molecules: linkers and rigids.
    std::vector<Molecule*> baseMolList;

    // Ids for the molecules in the graph.
    unsigned int graphId = 0;

    // Add all of the rigids to the worklist
    // Add the rigids to the hypergraph
    // Set the graph ID on the rigids AND the graph
    for (int r = 0; r < rigids.size(); r++)
    {
        baseMolList.push_back(rigids[r]);
        rigids[r]->setGraphID(graphId++);
        graph.AddNode(*rigids[r]);
    }

    // Add all of the linkers to the worklist AND the graph
    for (int ell = 0; ell < linkers.size(); ell++)
    {
        baseMolList.push_back(linkers[ell]);
        linkers[ell]->setGraphID(graphId++);
        graph.AddNode(*linkers[ell]);
    }

    //
    // Construct the set of 2-Molecules from the rigids and linkers.
    //
    const int BASE_LEVEL = 0;
    for (int m1 = 0; m1 < baseMolList.size() - 1; m1++)
    {
        for (int m2 = m1 + 1; m2 < baseMolList.size(); m2++)
        {
            std::vector<EdgeAggregator*>* newEdges = baseMolList[m1]->Compose(*baseMolList[m2]);

            HandleNewMolecules(level_queues[BASE_LEVEL+1], &queue_locks[BASE_LEVEL+1], *newEdges);

            for (int i = 0; i < newEdges->size(); i++)
            {
                delete (*newEdges)[i];
            }

            delete newEdges;
        }
    }

    // 1-Molecules and 2-Molecules have been processed.
    completed_level[BASE_LEVEL] = true;
    completed_level[BASE_LEVEL+1] = true;

    //
    // For each level, start a thread and compose the elements with the base set of molecules.
    //
    //                                           - 2 since 1-Molecules are complete
    for (int m = 1; m < HIERARCHICAL_LEVEL_BOUND - 2; m++)
    {
        //
        // Robert, spawn a thread with this function call
        //

        ProcessLevel(baseMolList, level_queues[m], level_queues[m+1],
                     queue_locks[m], queue_locks[m+1],
                     &completed_level[m-1], &completed_level[m]);
    }
}


//
// Main instantiation function for all figures stated in the given list;
// worklist technique to construct the graph
//
HyperGraph<Molecule, EdgeAnnotationT> Instantiator::Instantiate(std::vector<Linker*>& linkers,
                                                                std::vector<Rigid*>& rigids)
{
    // The worklist initialized to initial set of ground clauses from the figure
    std::queue<Molecule*> worklist;

    // Ids for the molecules in the graph.
    unsigned int graphId = 0;

    // Add all of the rigids to the worklist
    // Add the rigids to the hypergraph
    // Set the graph ID on the rigids AND the graph
    for (int r = 0; r < rigids.size(); r++)
    {
        worklist.push(rigids[r]);
        rigids[r]->setGraphID(graphId++);
        graph.AddNode(*rigids[r]);
    }
    
    // Add all of the linkers to the worklist AND the graph
    for (int ell = 0; ell < linkers.size(); ell++)
    {
        worklist.push(linkers[ell]);
        linkers[ell]->setGraphID(graphId++);
        graph.AddNode(*linkers[ell]);
    }

    //
    // Process all molecules until the worklist is empty: fixed point of combining of all linkers,
    // rigids, and new molecules.
    //
    while (!worklist.empty())
    {
        // Acquire the first molecule from the list for processing
        Molecule* currentMolecule = worklist.front();
        worklist.pop();

        if (HYPERGRAPH_CONSTR_DEBUG)
        {
            ds << "Working on: id(" << currentMolecule->getGraphID() << ")" << std::endl;
        }

        //
        // Take the current molecule and apply to all 'completed' molecules
        // (which are the nodes in the hypergraph)
        //
        int currentGraphSize = graph.size();
        for (int n = 0; n < currentGraphSize; n++)
        {
            std::vector<EdgeAggregator*>* newEdges;
            newEdges = currentMolecule->Compose(graph.vertices[n].data);

            HandleNewMolecules(worklist, 0, *newEdges);

            for (int i = 0; i < newEdges->size(); i++)
            {
                delete (*newEdges)[i];
            }
            delete newEdges;
        }
    }

    return graph;
}

//
// Add all new deduced clauses to the worklist if they have not been deduced before.
// If the given clause has been deduced before, update the hyperedges that were generated previously
//
// Forward Instantiation does not permit any cycles in the resultant graph.
//
void Instantiator::HandleNewMolecules(std::queue<Molecule*>& worklist,
                                      pthread_mutex_t* worklist_lock,
                                      std::vector<EdgeAggregator*>& newEdges)
{
    for (int e = 0; e < newEdges.size(); e++)
    {
        try
        {
            const Molecule& graphNode = graph.GetNode(*newEdges[e]->consequent);

            // Molecule is in the graph
            AddEdge(newEdges[e]->antecedent, graphNode, *newEdges[e]->annotation);
        }

        // The new consequent Molecule is not in the graph
        catch(unsigned int)
        {
            // Add the node to the graph and set its id
            if (graph.AddNode(*newEdges[e]->consequent))
            {
                newEdges[e]->consequent->setGraphID(graph.size() - 1);

                if (g_debug_output) std::cerr << "Added: " << newEdges[e]->consequent << std::endl; 
            }

            // Also add to the worklist; threaded for safety (if the lock is a valid pointer)
            if (worklist_lock) pthread_mutex_lock(worklist_lock);
            worklist.push(newEdges[e]->consequent);
            if (worklist_lock) pthread_mutex_unlock(worklist_lock);

            // Add the actual edge
            AddEdge(newEdges[e]->antecedent, *newEdges[e]->consequent, *newEdges[e]->annotation);
        }
    }
}
