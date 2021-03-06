#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <cstdio>
#include <sstream>
#include <cstdlib>


//
// Open Babel
//
#include <openbabel/obconversion.h>
#include <openbabel/mol.h>
#include <openbabel/generic.h>
#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/groupcontrib.h>


//
// This project molecular representation
//
#include "Atom.h"
#include "Bond.h"
#include "Molecule.h"
#include "Rigid.h"
#include "Linker.h"

//
// File processing in / out.
//
#include "OBWriter.h"
#include "Options.h"
#include "Validator.h"


//
// Synthesis-Based Functionality
//
#include "HyperGraph.h"
#include "EdgeAnnotation.h"
#include "Instantiator.h"

#include "PebblerHyperGraph.h"
#include "Utilities.h"
#include "IdFactory.h"
#include "Constants.h"


//
// Global set of linkers and rigids read from the input files.
//
std::vector<Linker*> linkers;
std::vector<Rigid*> rigids;

void Cleanup(std::vector<Linker*>& linkers, std::vector<Rigid*>& rigids);

bool splitMolecule(std::ifstream& infile, std::string& name,
                   std::string& prefix, std::string& suffix)
{
    prefix = "";
    suffix = "";

    std::string line = "";

    // Eat #### in large files (if it exists)
    eatWhiteLines(infile); 
    if (infile.peek() == '#')
    {
        getline(infile, line);
        name = line;
        eatWhiteLines(infile);
    }

    getline(infile, line);
    prefix += line + '\n';

    // Nothing left to read...
    if (infile.eof() || infile.fail()) return false;

    // Read the prefix (end indicated by END)
    while(line.find("END") == std::string::npos)
    {
        getline(infile, line);
        prefix += line + '\n';
    }

    // Add '$$$$' to the prefix.
    // prefix += "\n$$$$";

    // Set suffix equal to remainder of the file
    while (line.find("$$$$") == std::string::npos)
    {
        getline(infile, line);
        suffix += line + '\n';
    }

    return true;
}

Molecule* createLocalMolecule(OpenBabel::OBMol* mol, MoleculeT mType,
                              const std::string& name, std::string& suffix)
{
    //
    // Add the suffix as comment data to the actual OBMol object.  
    //
    OpenBabel::OBCommentData* cData = new OpenBabel::OBCommentData();
    cData->SetAttribute("Comment");
    cData->SetData(suffix);
    mol->SetData(cData);

    //
    // Create this particular molecule type based on the name of the file.
    //
    if (mType == LINKER)
    {
        return new Linker(mol, name);
    }
    else if (mType == RIGID)
    {
        return new Rigid(mol, name);
    }
    
    return 0;
}

void addMolecule(char type, Molecule* molecule)
{
    if (type == 'l')
    {
        linkers.push_back(static_cast<Linker*>(molecule));        
    }
    else if (type == 'r')
    {
        rigids.push_back(static_cast<Rigid*>(molecule));        
    }
}

void readMoleculeFile(const char* fileName)
{
    //
    // Input parser conversion functionality for Open babel
    //
    OpenBabel::OBConversion obConversion;
    obConversion.SetInFormat("SDF");

    //
    // Open the file, split the current molecule into Molecule Data (prefix)
    // and Our Data (Suffix)
    //
    std::ifstream infile;
    infile.open(fileName);

    std::string name = "UNKNOWN";
    std::string prefix = "";
    std::string suffix = "";
    
    while(splitMolecule(infile, name, prefix, suffix))
    {
        //
        // If the name of molecule is not given, overwrite it with the name of the file.
        //
        if (name == "UNKNOWN")
        {
           name = "####   ";
           name += fileName;
           name += "    ####";
        }

        if (g_debug_output) std::cerr << "Name: " << std::endl << name << std::endl;
        if (g_debug_output) std::cerr << "Prefix: " << std::endl << prefix << std::endl;
        if (g_debug_output) std::cerr << "Suffix: " << std::endl << suffix << std::endl;

        // Create and parse using Open Babel
        OpenBabel::OBMol* mol = new OpenBabel::OBMol();
        bool notAtEnd = obConversion.ReadString(mol, prefix);

        // Assign all needed data to the molecule (comment data)
        Molecule* local = createLocalMolecule(mol, fileName[0] == 'l' ? LINKER : RIGID,
                                              name, suffix);

        // calculate the molecular weight, H donors and acceptors and the plogp
        local->openBabelPredictLipinski();

        // add to logfile
        if (local->islipinskiPredicted())
        {
            std::ofstream logfile("synth_log_initial_fragments_logfile.txt",
                                  std::ofstream::out | std::ofstream::app); // append
            logfile << fileName << "\nMolWt = " << local->getMolWt() << "\n";
            logfile << "HBD = " << local->getHBD() << "\n";
            logfile << "HBA1 = " << local->getHBA1() << "\n";
            logfile << "logP = " << local->getlogP() << "\n";
            logfile << std::endl;
            logfile.close();
        }
        else std::cerr << "Main: predictLipinski failed somehow!" << endl;

        if (g_debug_output) std::cout << "Local: " << *local << "|" << std::endl;
    
        // Add to the linker or rigid list as needed.
        addMolecule(fileName[0], local); 
    }
}


//
// Parse each input data files
//
bool readInputFiles(const Options& options)
{
    for (std::vector<std::string>::const_iterator it = options.inFiles.begin();
         it != options.inFiles.end(); it++)
    {
        if ((*it)[0] != 'l' && (*it)[0] != 'r')
        {
            cerr << "Unexpected file prefix: \'" << (*it)[0]
                 << "\' with file " << *it << endl;
            return false;
        }

        readMoleculeFile((*it).c_str());
    }

    return true;
}


int main(int argc, char** argv)
{
    if (argc < 2)
    {
        std::cerr << "Usage: <program> [SDF-file-list] -o <output-file> -v <validation-file>"
                  << " -pool <#obgen-threads>" << std::endl;
        return 1;
    }

    //
    // Remove log files from a previous run.
    //
    std::string rmString = "rm synth_log*";
    system(rmString.c_str());

    //
    // Global options object.
    //
    Options options(argc, argv);
    if (!options.parseCommandLine())
    {
        std::cerr << "Command-line parsing failed; exiting." << std::endl;
        return 1;
    }

    // 
    // Output command-line option information
    //
    if (!Options::THREADED)
    {
        std::cout << "\'-threaded\' not specified, but will be used." << std::endl;
    }

    // Printing the specified Tanimoto value to the user.
    std::cerr << "Tanimoto Coefficient Threshold Specified: "
              << Options::TANIMOTO << std::endl;
    std::cerr << "OBGEN output thread pool size: "
              << Options::OBGEN_THREAD_POOL_SIZE << std::endl;


    if (!readInputFiles(options)) return 1;

    //
    // Bypass synthesis for acquiring information about the input fragments.
    //    
    if (g_calculate_lipinski_descriptors_for_input_fragments_only)
    {
        std::cout << "Calculated Lipinski Descriptors for input fragments, now exiting early."
                  << " (Flag set in Constants.h)" << std:: endl;
        return 0;
    }

    // Output object for the nodes of the hypergraph.
    OBWriter* writer = new OBWriter(Options::OBGEN_THREAD_POOL_SIZE);
    writer->InitializeFile(options.outFile);

    // The main object that performs synthesis.
    Instantiator instantiator(writer, cout);

    // Instantiation build the hypergraph; this is the main data structure for the
    // resultant molecules.
    // Also creates the hypergraph using threaded or non-threaded techniques.
    HyperGraph<Molecule, EdgeAnnotationT>* graph = instantiator.ThreadedInstantiate(linkers,
                                                                                    rigids);

    std::cout << "Hypergraph contains (" << graph->size() << ") nodes" << std::endl;
    std::cout << OBWriter::compliantMols.size()
              << " are Lipinski compliant molecules" << std::endl;

    //
    // Validate the molecules specified in the validation file (command-line -v)
    //
    Validator validator(OBWriter::compliantMols);
    validator.Validate(options.validationFile);

    // Deleting the writer will kill the thread pool.
    delete writer; 

    // For later: pebbling
    //PebblerHyperGraph<Molecule, EdgeAnnotationT> pebblerGraph = graph->GetPebblerHyperGraph();

    Cleanup(linkers, rigids);

std::cerr << "Exiting the main thread." << std::endl;

    return 0;
}

void Cleanup(std::vector<Linker*>& linkers, std::vector<Rigid*>& rigids)
{
    for (int ell = 0; ell < linkers.size(); ell++)
    {
        delete linkers[ell];
    }

    for (int r = 0; r < rigids.size(); r++)
    {
        delete rigids[r];
    }
}
