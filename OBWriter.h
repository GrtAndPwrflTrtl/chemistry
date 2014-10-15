#ifndef _OB_WRITER_GUARD
#define _OB_WRITER_GUARD 1


#include <vector>
#include <string>
#include <iostream>
#include <queue>
#include <pthread.h>


#include <openbabel/mol.h>


#include "Molecule.h"
#include "Thread_Pool.h"
#include "IdFactory.h"


//
// A class to dump all of the molecules to a file.
//
class OBWriter
{
  public:
    OBWriter(unsigned int numThreads);
    ~OBWriter();

    // static void InitializeFile(const char* fileName);
    void OutputMolecule(Molecule&);
    static int OutputSingleMolecule(std::string smiMol);
    static std::vector<OpenBabel::OBMol*> compliantMols;

    void IndicateSynthesisComplete();
    void InitiateOutputThreadPool();

    static void InitializeFile(const std::string& outFile);

    void write(std::vector<Molecule> molecules);

  private:
    unsigned int mCounter; 
    unsigned int mFailCounter; 
    bool writing_complete;
    bool writing_started;

    static pthread_mutex_t output_file_lock;
    static pthread_mutex_t valid_molecule_lock;
    static pthread_mutex_t id_lock;
    static IdFactory molIDmaker;
    static std::ofstream out;

    Thread_Pool<std::string, int>* pool;  

    void Initialize();
    static std::string ScrubAndConvertToSMI(OpenBabel::OBMol& mol);

    void ScrubAndExportSMI(std::vector<Molecule>& molecules);
    void CallsBeforeWriting(std::vector<Molecule>& molecules);
};

#endif
