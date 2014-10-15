#ifndef _INSTANTIATOR_GUARD
#define _INSTANTIATOR_GUARD 1


#include <vector>
#include <queue>
#include <iostream>
#include <memory>
#include <pthread.h>


#include "Molecule.h"
#include "Rigid.h"
#include "Linker.h"
#include "HyperGraph.h"
#include "EdgeAnnotation.h"
#include "IdFactory.h"
#include "OBWriter.h"


// threads require a struct to pass multiple arguments
struct Instantiator_ProcessLevel_Thread_Args
{
    int m; // level number
    void* this_pointer; // this pointer to calling class (Instantiator)
};

class Instantiator
{
  private:
    // Create necessary synthesis containers and init the linkers and rigids

    void InitializeBaseMolecules(const std::vector<Rigid*>& rigids,
                                 const std::vector<Linker*>& linkers,
                                 std::vector<Molecule*>& baseMolecules);

    // To generate unique molecular ids
    IdFactory moleculeIDFactory;

    // Contains all processed clauses and relationships amongst the clauses
    HyperGraph<Molecule, EdgeAnnotationT>* graph;

    // debug stream
    std::ostream& ds;

    void HandleNewMolecules(std::queue<Molecule*>& worklist,
                            pthread_mutex_t* wl_lock,
                            std::vector<EdgeAggregator*>& newEdges);

    void AddEdge(const std::vector<Molecule>& antecedent,
                 const Molecule& consequent,
                 const EdgeAnnotationT& annotation);

    bool AddNode(Molecule& mol);

    /*void ProcessLevel(std::vector<Molecule*>& baseMols,
                      std::queue<Molecule*>& inSet,
                      std::queue<Molecule*>& outSet,
                      pthread_mutex_t& in_lock,
                      pthread_mutex_t& out_lock,
                      bool* previousLevelComplete,
                      bool* thisLevelComplete);*/

    // Lock the hypergraph (for adding)
    pthread_mutex_t graph_lock;

    // A list of locks for all of the producer-consumer queues.
    pthread_mutex_t* queue_locks;

    // All of the hierarchical level threads.
    pthread_t* queue_threads;

    // Indicator that a level has completed processing.
    bool* completed_level;

    // The actual producer-consumer queue for each level.
    std::queue<Molecule*>* level_queues;

    // array of args for each level thread
    Instantiator_ProcessLevel_Thread_Args *arg_pointer;

    // set of linkers and rigids (1-molecules)
    std::vector<Molecule*> baseMolecules;

    // Molecules per level (count) for debug
    int* moleculeLevelCount;

    // For output of molecules on the fly.
    OBWriter* const writer;

  public:
    Instantiator(OBWriter*const obWriter, std::ostream& out = std::cout);

    // Main instantiation function for all linkers and rigidss; worklist technique to construct the graph

/*
    HyperGraph<Molecule, EdgeAnnotationT>* Instantiate(std::vector<Linker*>& linkers,
                                                      std::vector<Rigid*>& rigids);
*/

    // Main instantiation function for all linkers and rigidss; worklist technique to construct the graph
    HyperGraph<Molecule, EdgeAnnotationT>* ThreadedInstantiate(std::vector<Linker*>& linkers,
                                                              std::vector<Rigid*>& rigids);

    // thread must be implemented as friend class
    friend void *ProcessLevel(void * args); // worker thread
};

#endif
