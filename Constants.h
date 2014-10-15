#ifndef _CONSTANTS_GUARD
#define _CONSTANTS_GUARD 1


const unsigned int null = 0;


// Debugging constants
const bool DEBUG = true;
const bool HYPERGRAPH_CONSTR_DEBUG = true;
const bool g_debug_output = false;


const int THREAD_POOL_SIZE = 10;


// skip the entire synthesis, just output lipinski descriptors for
//  the input fragments to "initial_fragments_logfile.txt" and exit
const bool g_calculate_lipinski_descriptors_for_input_fragments_only = false;

// upper bound for level-based threading synthesis
extern unsigned int HIERARCHICAL_LEVEL_BOUND;

// Limiting factor on molecule generation.
extern double MOLWT_UPPERBOUND;
extern double HBD_UPPERBOUND;
extern double HBA1_UPPERBOUND;
extern double LOGP_UPPERBOUND;

#endif
