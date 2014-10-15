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
#include "OBWriter.h"
#include "Utilities.h"
#include "IdFactory.h"
#include "Constants.h"
#include "OBWriter.h"


// Remove as debug
#include "FragmentGraph.h"


Instantiator::Instantiator(OBWriter*const obWriter, std::ostream& out) : writer(obWriter), ds(out)
{
    graph = new HyperGraph<Molecule, EdgeAnnotationT>(HIERARCHICAL_LEVEL_BOUND + 1);

    // The hypergraph lock
    pthread_mutex_init(&graph_lock, NULL);

    // The threads and locks for the producer-consumer containers.
    queue_locks = new pthread_mutex_t[HIERARCHICAL_LEVEL_BOUND+1]; 
    queue_threads = new pthread_t[HIERARCHICAL_LEVEL_BOUND+1];
    completed_level = new bool[HIERARCHICAL_LEVEL_BOUND+1];
    level_queues = new std::queue<Molecule*>[HIERARCHICAL_LEVEL_BOUND+1];
    arg_pointer = new Instantiator_ProcessLevel_Thread_Args[HIERARCHICAL_LEVEL_BOUND+1];
    moleculeLevelCount = new int[HIERARCHICAL_LEVEL_BOUND + 1];

    for (int m = 1; m <= HIERARCHICAL_LEVEL_BOUND; m++)
    {
        // Initialize those container locks.
        pthread_mutex_init(&queue_locks[m], NULL);
	
        // Initialize the level-based threads.
        // ?

        // Initialize the fact that we have not computed this level. 
        completed_level[m] = false;

        // Initialize the level queues?
        // level_queues[m].clear();

        // set up arg structs
        arg_pointer[m].m=m;
        arg_pointer[m].this_pointer=this;

        moleculeLevelCount[m] = 0;
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

    graph->AddEdge(antecedent, consequent, annotation);

    pthread_mutex_unlock(&graph_lock);
}

//
// Add the hypernode to the hypergraph
//
bool Instantiator::AddNode(Molecule& mol)
{
    bool added = false;

//std::cout << "edding: " << mol << std::endl;
//std::cout << "Adding: " << mol.getFingerprint()->toString() << std::endl;

    pthread_mutex_lock(&graph_lock);

    if (graph->AddNode(mol))
    {
        mol.setUniqueIndexID(graph->size() - 1);
        added = true;
    }

    pthread_mutex_unlock(&graph_lock);

    return added;
}

//void Instantiator::ProcessLevel(std::vector<Molecule*>& baseMols,
//                                std::queue<Molecule*>& inSet,
//                                std::queue<Molecule*>& outSet,
//                                pthread_mutex_t& in_lock,
//                                pthread_mutex_t& out_lock,
//                                bool* previousLevelComplete,
//                                bool* thisLevelComplete)

void *ProcessLevel(void *ptr_void)
{
    //  unpacking arguments structure into mutiple local pointers
    Instantiator_ProcessLevel_Thread_Args * args = (Instantiator_ProcessLevel_Thread_Args *)ptr_void;
    int m = args->m; // level number
    Instantiator * This=(Instantiator *)args->this_pointer; // this pointer of calling class (Instantiator)

    //
    //  recast variables for local use (from the spawned thread record we were passed)
    //
    std::vector<Molecule*> *baseMols = &(This->baseMolecules);
    std::queue<Molecule*> *inSet = &(This->level_queues[m-1]);
    std::queue<Molecule*> *outSet = &(This->level_queues[m]);
    pthread_mutex_t *in_lock = &(This->queue_locks[m-1]);
    pthread_mutex_t *out_lock = &(This->queue_locks[m]);
    bool* previousLevelComplete = &(This->completed_level[m-1]);
    bool* thisLevelComplete = &(This->completed_level[m]);

    //
    // A structure for sleeping for 0.1 seconds
    //
    struct timespec sleepTime;
    struct timespec remTime;     // Remaining time
    sleepTime.tv_sec = 0;
    sleepTime.tv_nsec = 100000000L; // 0.1 seconds

    //
    // Keep consuming molecules as long as the previous level is incomplete or this
    // level queue contains molecules to process.
    //
    while (!(*previousLevelComplete) || !inSet->empty())
    {
        // Nothing to process, currently, but the level is incomplete.
        if (inSet->empty())
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
            pthread_mutex_lock(in_lock);
            Molecule* molToProcess = inSet->front();
            inSet->pop();
            pthread_mutex_unlock(in_lock);
            This->moleculeLevelCount[m-1]++;
//std::cout << "Took molecule off level " << m-1 << " queue" << std::endl;

            //
            // Process the molecule by composing it with all the base molecules.
            //
            for (std::vector<Molecule*>::iterator baseMol = baseMols->begin();
                 baseMol != baseMols->end(); baseMol++)
            {
                std::vector<EdgeAggregator*>* newEdges = molToProcess->Compose(**baseMol);

                This->HandleNewMolecules(*outSet, out_lock, *newEdges);

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

    std::cerr << "Level " << (m-1) << " created "
              << This->moleculeLevelCount[m-1] << " molecules." << std::endl; 

    std::cerr << "Level " << m << " complete." << std::endl; 
}

//
// Threaded construction of the hypergraph using a hierarchical list of threads and containers.
// We first construct the base case of 2-Molecules.
// Then, we inductively start constructing 3-Molecules, 4-Molecules, etc.
//
HyperGraph<Molecule, EdgeAnnotationT>* Instantiator::ThreadedInstantiate(
                                                                std::vector<Linker*>& linkers,
                                                                std::vector<Rigid*>& rigids)
{
    InitializeBaseMolecules(rigids, linkers, baseMolecules);

    // Add  all the base molecules to the hypergraph
    foreach_molecules(m_it, baseMolecules)
    {
        graph->AddNode(**m_it);
    }

    //
    // Construct the set of 2-Molecules from the rigids and linkers.
    //
    for (int m1 = 0; m1 < baseMolecules.size(); m1++)
    {
        for (int m2 = m1; m2 < baseMolecules.size(); m2++)
        {
            std::vector<EdgeAggregator*>* newEdges =
                                          baseMolecules[m1]->Compose(*baseMolecules[m2]);

            HandleNewMolecules(level_queues[2], &queue_locks[2], *newEdges);

            for (int i = 0; i < newEdges->size(); i++)
            {
                delete (*newEdges)[i];
            }

            delete newEdges;
        }
    }

    // 1-Molecules and 2-Molecules have been processed.
    completed_level[0] = true;
    completed_level[1] = true;
    completed_level[2] = true;

    // Indicate size of 1-M and 2-M lists
    moleculeLevelCount[1] = baseMolecules.size();

    //
    // For each level, start a thread and compose the elements with the base set of molecules.
    //
    for (int m = 3; m <= HIERARCHICAL_LEVEL_BOUND; m++)
    {
        if (~pthread_create(&queue_threads[m], NULL, ProcessLevel, (void*)&arg_pointer[m]))
	    {if (g_debug_output) {std::cout << "Level " << m << " thread created" << std::endl;}}
	else
            {if (g_debug_output) {std::cout << "Level " << m << " creation failed" << std::endl;}}
    }
    for (int m = 3; m <= HIERARCHICAL_LEVEL_BOUND; m++)
    {
	(void) pthread_join(queue_threads[m], NULL);
	if (g_debug_output) std::cout << "Level " << m << " thread removed" << std::endl;
    }

    std::cout << "Level\t" << "# Molecules" << std::endl; 
    for (int m = 1; m <= HIERARCHICAL_LEVEL_BOUND; m++)
    {
       std::cout << m << "\t" << moleculeLevelCount[m] << std::endl; 
    }

    // Tell the output engine we have completed synthesis.
    // This function then spins until the thread pool is complete. 
    this->writer->IndicateSynthesisComplete();

    return graph;
}


//
// Main instantiation function for all figures stated in the given list;
// worklist technique to construct the graph
//
/*
HyperGraph<Molecule, EdgeAnnotationT>* Instantiator::Instantiate(std::vector<Linker*>& linkers,
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
        graph->AddNode(*rigids[r]);
    }
    
    // Add all of the linkers to the worklist AND the graph
    for (int ell = 0; ell < linkers.size(); ell++)
    {
        worklist.push(linkers[ell]);
        linkers[ell]->setGraphID(graphId++);
        graph->AddNode(*linkers[ell]);
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
        int currentGraphSize = graph->size();
        for (int n = 0; n < currentGraphSize; n++)
        {
            std::vector<EdgeAggregator*>* newEdges;
            newEdges = currentMolecule->Compose(graph->vertices[n].data);

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
*/

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
/*
        std::cout << "Considering: "
                  << *newEdges[e]->consequent->getFingerprint() << std::endl;
*/
        try
        {
            const Molecule& graphNode = graph->GetNode(*newEdges[e]->consequent);

/*
            std::cout << "Molecule is already in the graph..." << std::endl;

            std::cout << "Failed to Add as redundant: "
                      << *newEdges[e]->consequent->getFingerprint() << std::endl;
*/

            // Molecule is in the graph
            AddEdge(newEdges[e]->antecedent, graphNode, *newEdges[e]->annotation);
        }

        // The new consequent Molecule is not in the graph
        catch(unsigned int)
        {
            // Add a node to the graph and set its id
            AddNode(*newEdges[e]->consequent);

/*
            std::cout << "Added: "
                      << *newEdges[e]->consequent->getFingerprint() << std::endl;
*/

            // Only output a molecule (in its complete form) on the fly if
            // we have a threaded version (which we always do now).
            if (worklist_lock) this->writer->OutputMolecule(*newEdges[e]->consequent);

            //
            // Also add to the worklist; threaded for safety (if the lock is a valid pointer)
            //
            if (worklist_lock) pthread_mutex_lock(worklist_lock);

            worklist.push(newEdges[e]->consequent);

            if (worklist_lock) pthread_mutex_unlock(worklist_lock);

            //std::cout << "Added molecule to a queue" << std:: endl;

            // Add the actual edge
            AddEdge(newEdges[e]->antecedent, *newEdges[e]->consequent, *newEdges[e]->annotation);
        }
    }
}


//
// Initialize the linkers and rigids as required; the baseMolecules list will then be
// used as a reference container throughout synthesis.
//
void Instantiator::InitializeBaseMolecules(const std::vector<Rigid*>& rigids,
                                           const std::vector<Linker*>& linkers,
                                           std::vector<Molecule*>& baseMolecules)
{
    // Clear the list just in case.
    baseMolecules.clear();

    // Assign the linkers and rigids unique ids; these correspond EXACTLY to the indices of
    // the containers used for determing molecular (non)-isomorphism.
    foreach_rigids(r_it, rigids)
    {
        (*r_it)->setUniqueIndexID(moleculeIDFactory.getNextId());
        baseMolecules.push_back(*r_it);
    }

    foreach_linkers(l_it, linkers)
    {
        (*l_it)->setUniqueIndexID(moleculeIDFactory.getNextId());
        baseMolecules.push_back(*l_it);
    }

    // The set of base molecules is static in the synthesis process; therefore,
    // we set the (static) reference base set of molecules in the Molecule class
    // so the corresponding molecular fingerprint graph can be constructed and compared.
    Molecule::SetBaseMoleculeInfo(baseMolecules, rigids.size(), linkers.size());

    // Each molecule will contain a reference count of the number of each specific
    // linker / rigid in the particular molecule.
    foreach_molecules(m_it, baseMolecules)
    {
        (*m_it)->initFragmentDevices();
        (*m_it)->initGraphRepresentation();
    }
}
