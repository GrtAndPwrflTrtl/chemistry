#ifndef _VALIDATOR_GUARD
#define _VALIDATOR_GUARD 1


#include <vector>


#include <openbabel/mol.h>


#include "HyperGraph.h"


//
// A class to perform validation on the synthesis:
//    Given (1) a hypergraph (formed from linkers and rigids)
//          (2) a set of Open Babel molecules,
//    Verify that the set of molecules are vertices in the hypergraph.
//
class Validator
{
  public:
    Validator(const HyperGraph<Molecule, EdgeAnnotationT>& g) : graph(g) {}
    bool Validate(OpenBabel::OBMol&);
    void Validate(const std::string& fileName);
    void Validate(std::vector<OpenBabel::OBMol>&);

  private:
    HyperGraph<Molecule, EdgeAnnotationT> graph; 
};

#endif
