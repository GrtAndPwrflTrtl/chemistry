#ifndef _LINKER_GUARD
#define _LINKER_GUARD 1

#include <vector>

#include <openbabel/mol.h>

#include "EdgeAggregator.h"
#include "Molecule.h"

class Linker : public Molecule
{
  public:
    Linker(OpenBabel::OBMol*, const std::string& name);
    Linker() {}
    ~Linker() {}

    bool operator==(const Linker& that) const
    {
        // Check type and then call base class equality.
        if (that.type != LINKER) return false;

        Molecule thisMol = *this;
        Molecule thatMol = that;

        return thisMol == thatMol;
    }

  protected:
    virtual void parseAppendix(std::string& suffix);

  private:
};

#endif
