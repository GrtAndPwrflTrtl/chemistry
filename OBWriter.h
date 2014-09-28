#ifndef _OB_WRITER_GUARD
#define _OB_WRITER_GUARD 1


#include <vector>
#include <iostream>


#include "Molecule.h"
#include "Thread_Pool.h"


//
// A class to dump all of the molecules to a file.
//
class OBWriter
{
  public:
    OBWriter(const char*);
    void write(std::vector<Molecule> molecules);

  private:
    std::ofstream out;
    unsigned int mCounter; 

    void ScrubAndExportSMI(std::vector<Molecule>& molecules);
    void CallsBeforeWriting(std::vector<Molecule>& molecules);
};

#endif
