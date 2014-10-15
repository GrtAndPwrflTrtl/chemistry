#ifndef _RIGID_GUARD
#define _RIGID_GUARD 1

#include <vector>

#include <openbabel/mol.h>

#include "EdgeAggregator.h"
#include "Molecule.h"

class Rigid : public Molecule
{
  public:
    Rigid(OpenBabel::OBMol* obmol, const std::string& name);
    Rigid() : uniqueFragmentID(-1) {}
    ~Rigid() {}

    bool operator==(const Rigid& that) const
    {
        return this->getUniqueIndexID() == that.getUniqueIndexID();
    }

    unsigned int getFragmentId() const { return this->uniqueFragmentID; }  

  protected:  
    virtual void parseAppendix(std::string& suffix);
    
  private:
    unsigned int uniqueFragmentID;
};

#endif
