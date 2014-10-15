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
    Linker() : uniqueFragmentID(-1) {}

    ~Linker() {}

    bool operator==(const Linker& that) const
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
