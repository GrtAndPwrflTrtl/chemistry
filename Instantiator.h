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

class Instantiator
{
  private:
    // Contains all processed clauses and relationships amongst the clauses
    HyperGraph<Molecule, EdgeAnnotationT> graph;

    // debug stream
    std::ostream& ds;

    void HandleNewMolecules(std::queue<Molecule*>& worklist,
                            pthread_mutex_t* wl_lock,
                            std::vector<EdgeAggregator*>& newEdges);

    void AddEdge(const std::vector<Molecule>& antecedent,
                 const Molecule& consequent,
                 const EdgeAnnotationT& annotation);

    void ProcessLevel(std::vector<Molecule*>& baseMols,
                      std::queue<Molecule*>& inSet,
                      std::queue<Molecule*>& outSet,
                      pthread_mutex_t& in_lock,
                      pthread_mutex_t& out_lock,
                      bool* previousLevelComplete,
                      bool* thisLevelComplete);

    pthread_mutex_t graph_lock;
    pthread_mutex_t* queue_locks;
    pthread_t* queue_threads;
    bool* completed_level;
    std::queue<Molecule*>* level_queues;


  public:
    Instantiator(std::ostream& out = std::cout);

    // Main instantiation function for all linkers and rigidss; worklist technique to construct the graph
    HyperGraph<Molecule, EdgeAnnotationT> Instantiate(std::vector<Linker*>& linkers,
                                                      std::vector<Rigid*>& rigids);

    // Main instantiation function for all linkers and rigidss; worklist technique to construct the graph
    HyperGraph<Molecule, EdgeAnnotationT> ThreadedInstantiate(std::vector<Linker*>& linkers,
                                                              std::vector<Rigid*>& rigids);

};

#endif
