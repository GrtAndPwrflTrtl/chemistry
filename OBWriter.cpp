#include <vector>
#include <iostream>
#include <sstream>

#include <openbabel/obconversion.h>


#include "OBWriter.h"
#include "Molecule.h"
#include "Rigid.h"
#include "Linker.h"
#include "obgen.h"
#include "Thread_Pool.h"
#include "Constants.h"

// Static allocation of the thread pool.
// Thread_Pool<OpenBabel::OBMol*, bool> OBWriter::pool(THREAD_POOL_SIZE, OBGen::obgen);


// ****************************************************************************

OBWriter::OBWriter(const char* outFile) : mCounter(0)
{
    out.open(outFile);

    if (out.fail()) throw "Output stream opening failed.";
}

// ****************************************************************************

// This function is no longer called.

int OBWriter::CountComplexLipinski(std::vector<Molecule>& molecules)
{
    int numComplex = 0;

    for (std::vector<Molecule>::iterator it = molecules.begin(); it != molecules.end(); it++)
    {
        // The output molecules must be lipinski compliant (and non-linker, non-rigid).
        if (it->IsComplex() && it->isLipinskiCompliant())
        {
            numComplex++;
        }
    }

    return numComplex;
}


// ****************************************************************************

void OBWriter::ScrubAndExportSMI(std::vector<Molecule>& molecules)
{
    // init
    OpenBabel::OBConversion SDF_conv(&std::cin, &std::cout); // conversion to/from SDF (has xyz coords)
    OpenBabel::OBConversion SMI_conv(&std::cin, &std::cout); // conversion to/from SMI (no xyz coords)

    if(SMI_conv.SetInAndOutFormats("SMI","SMI") && SDF_conv.SetInAndOutFormats("SDF","SDF")) // set conversion type(s) and verify it worked
    {
        // more init
        std::string s; // temporary buffer
        int i=1; // counter (for debugging output)
        std::ofstream logfile("ScrubAndExportSMI_logfile.txt", std::ofstream::out | std::ofstream::app); // append
        if (logfile.is_open())
            std::cout  << "ScrubAndExportSMI: logfile.is_open";
        else
            std::cerr << "ScrubAndExportSMI: Unable to open logfile";

        for (std::vector<Molecule>::iterator it = molecules.begin(); it != molecules.end(); it++) // iterate through molecules
        {
            if (it->IsComplex() && it->isLipinskiCompliant()) // otherwise we dont care about it anyway, so skip to save time/effort
            {
                // Set up data struct and display basic information
                if (g_debug_output) std::cout << "ScrubAndExportSMI: getOpenBabelMol..." << std::endl;
                OpenBabel::OBMol* mol = it->getOpenBabelMol();  // OBMol mol;
                std::cout << "ScrubAndExportSMI: #" << i << ": NumAtoms=" << mol->NumAtoms() << ", NumBonds=" << mol->NumBonds() << std::endl;
                logfile << "ScrubAndExportSMI: #" << i << ": NumAtoms=" << mol->NumAtoms() << ", NumBonds=" << mol->NumBonds() << std::endl;
		//std::cout << it->toString() << std::endl; // more detailed data

                // pre-emptive extra run of OBGen, seems to stop segmentation fault
                if (g_debug_output) std::cout << "ScrubAndExportSMI: OBGen(fast)..." << std::endl;
                OBGen::fast_obgen(mol);

                // log "before" molecule
                //logfile << "molecule #" << i << ", SDF format, before:" << std::endl;
                //s=SDF_conv.WriteString(mol);
                //logfile << s;

                // Write to then read from SMI; should remove xyz coords
                if (g_debug_output) std::cout << "ScrubAndExportSMI: WriteString:" << std::endl;
                s=SMI_conv.WriteString(mol);
                logfile << "molecule #" << i << ", SMI: " << s; // log "SMI" version
                if(SMI_conv.ReadString(mol, s))
                   {if (g_debug_output) {std::cout << "ScrubAndExportSMI: ReadString - successful" << std::endl;}}
                else
                   std::cout << "ScrubAndExportSMI: ReadString - failed" << std::endl;

                // log "after" molecule
                if (g_debug_output) std::cout << "molecule #" << i << ", SDF format, after:" << std::endl;
                //logfile << "molecule #" << i << ", SDF format, after:" << std::endl;
                s=SDF_conv.WriteString(mol);
                if (g_debug_output) std::cout << s;
                //logfile << s;
            }
            if (g_debug_output) std::cout << "-----------------------------------------------------------------" << std::endl;
            logfile << "-----------------------------------------------------------------" << std::endl;
        i++;
        }
        logfile.close();
    }
    else
        std::cerr << "SetInAndOutFormats failed!" << std::endl;
}

// ****************************************************************************

void OBWriter::CallsBeforeWriting(std::vector<Molecule>& molecules)
{
    int counter = 0;
    for (std::vector<Molecule>::iterator it = molecules.begin(); it != molecules.end(); it++)
    {
        // The output molecules must be lipinski compliant (and non-linker, non-rigid).
        if (it->IsComplex() && it->isLipinskiCompliant())
        {
            OpenBabel::OBMol* obmol = it->getOpenBabelMol();

            std::cerr << "Converting molecule " << counter++ << " with obgen" << std::endl; 

            // pool.push(obmol);
            OBGen::obgen(obmol);
        }
    }
}


// ****************************************************************************

void OBWriter::write(std::vector<Molecule> molecules)
{
    //
    // Refine the list of molecules to those that are lipinski compliant.
    for (std::vector<Molecule>::iterator it = molecules.begin(); it != molecules.end(); it++)
        {
            it->predictLipinski();
        }

    // int numToProcess = CountComplexLipinski(molecules);

    std::cout << "OBWriter::write: ScrubAndExportSMI(molecules)..." << std::endl;
    ScrubAndExportSMI(molecules);

    std::cout << "OBWriter::write: CallsBeforeWriting(molecules)..." << std::endl;
    CallsBeforeWriting(molecules);

    // while(pool.out_q_size() < numToProcess);

    std::cerr << "obgen processing complete." << std::endl;


    //
    // Open Babel conversion from SDF to SMI format
    //
/*
    std::ostringstream oss;
    OpenBabel::OBConversion SDFtoSMI(&std::cin, &oss);
    conv.SetInAndOutFormats("SDF","SMI");
*/
    OpenBabel::OBConversion toSDF(&std::cin, &this->out);
    toSDF.SetOutFormat("SDF");


    for (std::vector<Molecule>::iterator it = molecules.begin(); it != molecules.end(); it++)
    {
        // The output molecules must be lipinski compliant (and non-linker, non-rigid).
        // Change IsLipinski to some function call to recalculate the OPENBABEL lipinski comliance check
        if (it->IsComplex() && it->isLipinskiCompliant())
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
            const std::vector<Rigid*>& rigids = it->getRigids();
            const std::vector<Linker*>& linkers = it->getLinkers();

            for (std::vector<Rigid*>::const_iterator rit = rigids.begin(); rit != rigids.end(); rit++)
            {
                this->out << (*rit)->getName() << std::endl;
            }

            for (std::vector<Linker*>::const_iterator lit = linkers.begin(); lit != linkers.end(); lit++)
            {
                this->out << (*lit)->getName() << std::endl;
            }
/*
            //
            // (1) Convert SDF to SMI
            // (2) Call OBGEN to regenerate 3D coordinates
            //
            // We go through this process since we don't want to change any internals of the OBmol;
            // we are acting on a copy, per se.
            //
            // stream written to
            std::ostringstream oss;

            // Setting up open babel conversion.
            OpenBabel::OBConversion SDFtoSMI(&std::cin, &oss);
            SDFtoSMI.SetInAndOutFormats("SDF","SMI");

            // The molecule we are writing.
            OpenBabel::OBMol* obmol = it->getOpenBabelMol();

            // Write the molecule to the temp stream
            SDFtoSMI.Write(obmol);

            //
            // Read the molecule as SMI
            // 
            OpenBabel::OBConversion SMItoSDF(oss.str(), &this->out);
            conv.SetInAndOutFormats("SDF","SMI");

*/
            //
            // Take the SMI and convert it back to SDF by populating the coordinates.
            //
            OpenBabel::OBMol* obmol = it->getOpenBabelMol();
            // OBGen::obgen(obmol);
            toSDF.Write(obmol);
        }
    }
}

// ****************************************************************************
