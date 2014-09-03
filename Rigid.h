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
    Rigid() {}
    ~Rigid() {}

    bool operator==(const Rigid& that) const
    {
        // Check type and then call base class equality.
        if (that.type != RIGID) return false;

        Molecule thisMol = *this;
        Molecule thatMol = that;

       return thisMol == thatMol;
    }

  protected:  
    virtual void parseAppendix(std::string& suffix);
    
  private:
};

#endif
