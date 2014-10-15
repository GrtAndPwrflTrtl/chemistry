#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include <pthread.h>
#include <unistd.h>


#include <openbabel/mol.h>
#include <openbabel/obconversion.h>


#include "OBWriter.h"
#include "Molecule.h"
#include "Rigid.h"
#include "Linker.h"
#include "obgen.h"
#include "Thread_Pool.h"
#include "Constants.h"
#include "Utilities.h"
#include "IdFactory.h"
#include "Options.h"


// Static Definitions
pthread_mutex_t OBWriter::valid_molecule_lock;
pthread_mutex_t OBWriter::output_file_lock;
pthread_mutex_t OBWriter::id_lock;
IdFactory OBWriter::molIDmaker(1000);
std::ofstream OBWriter::out;
std::vector<OpenBabel::OBMol*> OBWriter::compliantMols;

// ****************************************************************************

OBWriter::OBWriter(unsigned int threadCount) : mCounter(0),
                                               mFailCounter(0),
                                               writing_complete(false),
                                               writing_started(false) 
{
    pthread_mutex_init(&output_file_lock, NULL);
    pthread_mutex_init(&id_lock, NULL);

    // Create the thread pool
    pool = new Thread_Pool<std::string, int>(threadCount, OBWriter::OutputSingleMolecule);
}

// ****************************************************************************

OBWriter::~OBWriter()
{
    // Killing the thread pool to force all threads to join.
    delete pool;
}

// ****************************************************************************

void OBWriter::InitializeFile(const std::string& outFile)
{
    out.open(outFile.c_str());

    if (out.fail()) throw "Output stream opening failed.";
}

// ****************************************************************************

void OBWriter::Initialize()
{
    pthread_mutex_init(&OBWriter::valid_molecule_lock, NULL);
    pthread_mutex_init(&OBWriter::output_file_lock, NULL);
    pthread_mutex_init(&OBWriter::id_lock, NULL);
}

// ****************************************************************************

void OBWriter::IndicateSynthesisComplete()
{
    std::cerr << "Synthesis is complete; writing continues." << std::endl;

    // Spin until writing is complete.
    while (writing_started && mCounter > pool->out_q_size())
    {
        // Sleep 5 seconds; obgen takes a while.
        sleep(5);
    }
    std::cerr << "Writing of the molecules with obgen is complete." << std::endl;
}

// ****************************************************************************

void OBWriter::OutputMolecule(Molecule& mol)
{
    // If this is the first call to output, save the fact we are writing
    if (!writing_started)
    {
        OBWriter::Initialize();
        writing_started = true;
    }

    //
    // Output molecule
    //
    //std::cerr << "Current counter: fail(" << this->mFailCounter << "); pass ("
    //          << this->mCounter << ")" << std::endl;

    //
    // Process the molecule for output
    //
    // (1) Lock around open babel; this makes a copy with the copy constructor.
    pthread_mutex_lock(&Molecule::openbabel_lock);  
    OpenBabel::OBMol theMol = *(mol.getOpenBabelMol());

    // The molecule must be Lipinski compliant (using Open Babel)
    // We use the copy as not to disrupt the approximations we use during synthesis.
    if (!Molecule::isOpenBabelLipinskiCompliant(theMol))
    {
        this->mFailCounter++;
        pthread_mutex_unlock(&Molecule::openbabel_lock);
        return;
    }
    else this->mCounter++;

    // (2) Export to SMI
    std::string smiMol = OBWriter::ScrubAndConvertToSMI(theMol);
    pthread_mutex_unlock(&Molecule::openbabel_lock);

    //
    // (3) Add the molecule to the queue for processing.
    //
    pool->push(smiMol);
}

// ****************************************************************************
/*
void OBWriter::Initialize()
{
    pthread_mutex_init(&OBWriter::valid_molecule_lock, NULL);
    pthread_mutex_init(&OBWriter::output_file_lock, NULL);
    pthread_mutex_init(&OBWriter::id_lock, NULL);
}
*/
// ****************************************************************************

int OBWriter::OutputSingleMolecule(std::string smiMol)
{
    //
    // Generate a unique identification number for this molecule.
    //
    pthread_mutex_lock(&OBWriter::id_lock);
    unsigned int id = molIDmaker.getNextId();
    pthread_mutex_unlock(&OBWriter::id_lock);

    // 
    // Write SMI to a temp file
    //
    std::ostringstream smiFileOut;
    smiFileOut << "synth_log_smi_" << id << ".smi";
    std::string smiFile = smiFileOut.str();

    std::ofstream smiOut;
    smiOut.open(smiFile.c_str());
    smiOut << smiMol;
    smiOut.close();

    std::ostringstream molFileOut;
    molFileOut << "synth_log_temp_mol_" << id << ".sdf";
    std::string molFile = molFileOut.str();

    //
    // we suppress all stderr messages with 2> and keep stdout with output.
    //
    std::ostringstream obgenCall;
    obgenCall << "./obgen " << smiFileOut.str() << " 2> /dev/null " << " > " << molFile.c_str();

    std::cerr << "Calling: " << obgenCall.str() << std::endl;

    // Spawn obgen process to work on the temp file.
    system(obgenCall.str().c_str());

    //
    // Convert to SDF and append to the final output file; maintain molecule for validation.
    //
    ifstream in;
    in.open(molFile.c_str());

    //
    // We keep the SDF version of the synthesized molecule
    //
    // Begin open babel usage
    pthread_mutex_lock(&Molecule::openbabel_lock);

    OpenBabel::OBMol* mol = new OpenBabel::OBMol();
    OpenBabel::OBConversion SDF_conv;
    SDF_conv.SetInAndOutFormats("SDF", "SDF");
    SDF_conv.ReadFile(mol, molFile.c_str());

    //
    // Append output to a total output file.
    //
    pthread_mutex_lock(&output_file_lock);
    OBWriter::out << SDF_conv.WriteString(mol);
    pthread_mutex_unlock(&output_file_lock);

    // End open babel usage
    pthread_mutex_unlock(&Molecule::openbabel_lock);

    // Save the valid molecule
    pthread_mutex_lock(& OBWriter::valid_molecule_lock);
    OBWriter::compliantMols.push_back(mol);
    pthread_mutex_unlock(& OBWriter::valid_molecule_lock);

    // Close the input file MOL file.
    in.close();

    std::cerr << molFile << " completed." << std::endl;

    return 0;
}

// ****************************************************************************

std::string OBWriter::ScrubAndConvertToSMI(OpenBabel::OBMol& theMol)
{
    OpenBabel::OBConversion SMI_conv(&std::cin, &std::cout);

    // set conversion type(s) and verify it worked
    if(!SMI_conv.SetInAndOutFormats("SMI","SMI"))
    {
        throw "SetInAndOutFormats failed!";
    }

//    std::string s; // temporary buffer

/*
    std::ofstream logfile("synth_log_ScrubAndExportSMI_logfile.txt",
                           std::ofstream::out | std::ofstream::app); // append
*/

    // Pre-emptive extra run of OBGen, seems to stop segmentation fault
    OBGen::fast_obgen(&theMol);

    // Write to then read from SMI; should remove xyz coords
    //if (g_debug_output) std::cout << "ScrubAndExportSMI: WriteString:" << std::endl;

    return SMI_conv.WriteString(&theMol);

/*
    s = SMI_conv.WriteString(mol);

    SMI_conv.ReadString(mol, s);

    s = SDF_conv.WriteString(mol);
*/
}

// ****************************************************************************

void OBWriter::ScrubAndExportSMI(std::vector<Molecule>& molecules)
{
    OpenBabel::OBConversion SDF_conv(&std::cin, &std::cout); // conversion to/from SDF (has xyz coords)
    OpenBabel::OBConversion SMI_conv(&std::cin, &std::cout); // conversion to/from SMI (no xyz coords)

    // set conversion type(s) and verify it worked
    if(!SMI_conv.SetInAndOutFormats("SMI","SMI") || !SDF_conv.SetInAndOutFormats("SDF","SDF"))
    {
        std::cerr << "SetInAndOutFormats failed!" << std::endl;
        return;    
    }

    std::string s; // temporary buffer
    int i = 1;     // counter (for debugging output)

    std::ofstream logfile("ScrubAndExportSMI_logfile.txt", std::ofstream::out | std::ofstream::app); // append

    if (logfile.is_open())
    {
        std::cout  << "ScrubAndExportSMI: logfile.is_open";
    }
    else
    {
        std::cerr << "ScrubAndExportSMI: Unable to open logfile";
    }

    //
    // Iterate through all molecules.
    //
    for (std::vector<Molecule>::iterator it = molecules.begin(); it != molecules.end(); it++)
    {
        // Set up data struct and display basic information
        if (g_debug_output) std::cout << "ScrubAndExportSMI: getOpenBabelMol..." << std::endl;

        OpenBabel::OBMol* mol = it->getOpenBabelMol();  // OBMol mol;

        std::cout << "ScrubAndExportSMI: #" << i 
                  << ": NumAtoms=" << mol->NumAtoms()
                  << ", NumBonds=" << mol->NumBonds() << std::endl;
        logfile << "ScrubAndExportSMI: #" << i
                << ": NumAtoms=" << mol->NumAtoms()
                << ", NumBonds=" << mol->NumBonds() << std::endl;
        //std::cout << it->toString() << std::endl; // more detailed data

        // Pre-emptive extra run of OBGen, seems to stop segmentation fault
        if (g_debug_output) std::cout << "ScrubAndExportSMI: OBGen(fast)..." << std::endl;
        OBGen::fast_obgen(mol);

        // log "before" molecule
        //logfile << "molecule #" << i << ", SDF format, before:" << std::endl;
        //s=SDF_conv.WriteString(mol);
        //logfile << s;

        // Write to then read from SMI; should remove xyz coords
        if (g_debug_output) std::cout << "ScrubAndExportSMI: WriteString:" << std::endl;

        s = SMI_conv.WriteString(mol);

        logfile << "molecule #" << i << ", SMI: " << s; // log "SMI" version

        if (SMI_conv.ReadString(mol, s))
        {
            if (g_debug_output) std::cout << "ScrubAndExportSMI: ReadString - successful" << std::endl;

            // log "after" molecule
            if (g_debug_output) std::cout << "molecule #" << i << ", SDF format, after:" << std::endl;
                
            //logfile << "molecule #" << i << ", SDF format, after:" << std::endl;
            s = SDF_conv.WriteString(mol);

            if (g_debug_output) std::cout << s;
            //logfile << s;
        }
        else if (g_debug_output) std::cout << "ScrubAndExportSMI: ReadString - failed" << std::endl;

        if (g_debug_output)
        {
            std::cout << "-----------------------------------------------------------------" << std::endl;
        }

        logfile << "-----------------------------------------------------------------" << std::endl;
        i++;
    }

    logfile.close();
}

// ****************************************************************************

void OBWriter::CallsBeforeWriting(std::vector<Molecule>& molecules)
{
    int counter = 0;
    for (std::vector<Molecule>::iterator it = molecules.begin(); it != molecules.end(); it++)
    {
        OpenBabel::OBMol* obmol = it->getOpenBabelMol();

        // if (g_debug_output)
        std::cerr << "Converting molecule " << counter++ << " with obgen" << std::endl; 

        OBGen::obgen(obmol);
    }
}


// ****************************************************************************

void OBWriter::write(std::vector<Molecule> molecules)
{
    //
    // Refine the list of molecules to those that are lipinski compliant.
    //
    std::vector<Molecule> synthMolecules;
    for (std::vector<Molecule>::iterator it = molecules.begin(); it != molecules.end(); it++)
    {
        it->openBabelPredictLipinski();

        if (it->IsComplex() && it->isLipinskiCompliant())
        {
            synthMolecules.push_back(*it);
        }
    }

    //
    // Convert to SMI and then use obgen to regenerate 3D molecules.
    //
    std::cout << "OBWriter::write: ScrubAndExportSMI(synthMolecules)..." << std::endl;
    ScrubAndExportSMI(synthMolecules);

    std::cout << "OBWriter::write: CallsBeforeWriting(synthMolecules)..." << std::endl;
    CallsBeforeWriting(synthMolecules);

    // Converter to output the synthesized molecules
    OpenBabel::OBConversion toSDF(&std::cin, &this->out);
    toSDF.SetOutFormat("SDF");

    //
    // Process all refined molecules
    //
    for (std::vector<Molecule>::iterator it = synthMolecules.begin();
         it != synthMolecules.end();
         it++)
    {
        //
        // Print the molecule number
        //
        this->out << "#### ";
        this->out << mCounter++;
        this->out << " ####";

        //
        // Print the names of all the linkers / rigids.
        //
        std::vector<Rigid*> rigids;
        it->getRigids(rigids);
        std::vector<Linker*> linkers;
        it->getLinkers(linkers);

        foreach_rigids(r_it, rigids)
        {
            this->out << (*r_it)->getName() << std::endl;
        }

        foreach_linkers(l_it, linkers)
        {
            this->out << (*l_it)->getName() << std::endl;
        }
        
        //
        // Take the SMI and convert it back to SDF by populating the coordinates.
        //
        OpenBabel::OBMol* obmol = it->getOpenBabelMol();
        toSDF.Write(obmol);
    }
}

// ****************************************************************************
