#ifndef _CONSTANTS_GUARD
#define _CONSTANTS_GUARD 1


const unsigned int null = 0;

// Debugging constants
const bool DEBUG = true;
const bool HYPERGRAPH_CONSTR_DEBUG = true;

// Limiting factor on molecule generation.
const double MAX_DALTON_WEIGHT = 500.0;

const int THREAD_POOL_SIZE = 10;

// print debugging output to command line?
const bool g_debug_output=true;

// skip the entire synthesis, just output lipinski descriptors for
//  the input fragments to "initial_fragments_logfile.txt" and exit
const bool g_calculate_lipinski_descriptors_for_input_fragments_only=false;

#endif
