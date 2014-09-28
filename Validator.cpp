#include <vector>
#include <fstream>


#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/fingerprint.h>


#include "HyperGraph.h"
#include "Validator.h"
#include "Options.h"

//
// Validate a single molecule.
//
bool Validator::Validate(OpenBabel::OBMol& validationMol)
{
    std::vector<unsigned int> validationFP;
    bool returnval = false;

    if (g_debug_output) std::cerr << "Atoms: " << validationMol.NumAtoms() << std::endl;
    if (g_debug_output) std::cerr << "Bonds: " << validationMol.NumBonds() << std::endl;

    // Create the fingerprint for the validation molecule
    OpenBabel::OBFingerprint* fpType = OpenBabel::OBFingerprint::FindFingerprint("");

    if (g_debug_output) std::cerr << "Using Fingerprint: " << fpType->Description() << std::endl;

    // Acquire the fingerprint of the validation molecule so we can use it for Tanimoto comparison.
    fpType->GetFingerprint(&validationMol, validationFP);

    if (g_debug_output)
    {
        std::cerr << "Validation: " << std::endl;
        for (std::vector<unsigned int>::const_iterator it = validationFP.begin();
             it != validationFP.end(); it++)
        {
            std::cerr << *it << "|";
        }
        std::cerr << std::endl;
    }

    //
    // Check this validation fingerprint against all of fingerprints in the hypergraph
    // 
    // Tracks the largest tanimoto coefficient of the synthesized molecule (and index?)
    //
    double maxTanimoto = -1;
    int maxIndex = -1;
    for (int v = 0; v < graph.vertices.size(); v++)
    {
        // Acquire the fingerprint of the hypergraph molecule.
        std::vector<unsigned int> hgFP;
        graph.vertices[v].data.GetFingerprint(hgFP);

        // Debug output of the fingerprint.
        if (g_debug_output)
        {
            std::cerr << "Hypergraph: " << std::endl;
            for (std::vector<unsigned int>::const_iterator it = hgFP.begin(); it != hgFP.end(); it++)
            {
                std::cerr << *it << "|";
            }
            std::cerr << std::endl;
        }

        // Acquire the tanimoto coefficient for the validation molecule and
        // the synthesized molecules in the hypergraph.
        double tanimoto = OpenBabel::OBFingerprint::Tanimoto(validationFP, hgFP);

        if (tanimoto > maxTanimoto)
        {
            maxTanimoto = tanimoto;
            maxIndex = v;
        }

        if (g_debug_output) std::cerr << "Tanimoto: " << tanimoto << std::endl;

        if (tanimoto > (double)Options::TANIMOTO) 
        {
            returnval = true;
            break;
        }
    }


    std::ofstream logfile("Max_Tanimoto_logfile.txt", std::ofstream::out | std::ofstream::app); // append
    logfile << "Validation Molecule: " << validationMol.GetTitle() << " with ";
    logfile << "Synth Molecule: " << graph.vertices[maxIndex].data.getName() << "\n";
    logfile << maxIndex << ": maxTanimoto = " << maxTanimoto;
    // if (returnval) logfile << " - Validated\n";
    // else logfile << " - Failed to Validate\n";
    logfile << std::endl;
    logfile.close();

    return returnval;
}


//
// Validate a list of molecules.
//
void Validator::Validate(std::vector<OpenBabel::OBMol>& molecules)
{
    for (std::vector<OpenBabel::OBMol>::iterator it = molecules.begin(); it != molecules.end(); it++)
    {
        if (!Validate(*it))
        {
            std::cerr << "Failed to validate: " << it->GetTitle() << std::endl;
        }
    }
}


//
// For validation purposes, we use the original files in MOL format (before
// being stripped for linkers and rigids).
//
// (1) Parse the input files for all molecules to validate. 
// (2) Acquire all of the hypergraph molecules and perform obgen to re-acquire hydrogen bonds.
// (3) Compared the molecules using the Tanimoto similarity. (Get fingerprints of both molecules and compare).
//
void Validator::Validate(const std::string& fileName)
{
    //
    // Input parser conversion functionality for Open babel
    //
    OpenBabel::OBConversion obConversion;
    obConversion.SetInFormat("MOL2");

    // The molecules to validate.
    std::vector<OpenBabel::OBMol> molecules;
   
    std::cerr << "Reading Validation file " << fileName << std::endl;
 
    //
    // Read all of the OBMol objects using Open Babel; using their sample code style for reading.
    //
    OpenBabel::OBMol* mol = new OpenBabel::OBMol();
    bool notAtEnd = obConversion.ReadFile(mol, fileName.c_str());
    molecules.push_back(*mol);

    while(notAtEnd)
    {
        // Create and parse using Open Babel
        OpenBabel::OBMol* mol = new OpenBabel::OBMol();
        notAtEnd = obConversion.Read(mol);
        if (notAtEnd) molecules.push_back(*mol);
    }

    // Validate all molecules in the given file.
    Validate(molecules);

/*
    // Destroy all the newly-created molecules we validated.
    for (std::vector<OpenBabel::OBMol>::iterator it = molecules.begin(); it != molecules.end(); it++)
    {
        delete &(*it);
    }
*/
}


