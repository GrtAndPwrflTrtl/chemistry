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

std::cerr << "Atoms: " << validationMol.NumAtoms() << std::endl;
std::cerr << "Bonds: " << validationMol.NumBonds() << std::endl;

    // Create the fingerprint for the validation molecule
    OpenBabel::OBFingerprint* fpType = OpenBabel::OBFingerprint::FindFingerprint("");

std::cerr << "Using Fingerprint: " << fpType->Description() << std::endl;

    fpType->GetFingerprint(&validationMol, validationFP);

std::cerr << "Validation: " << std::endl;
for (std::vector<unsigned int>::const_iterator it = validationFP.begin(); it != validationFP.end(); it++)
{
    std::cerr << *it << "|";
}
std::cerr << std::endl;


    //
    // Check this validation fingerprint against all of fingerprints in the hypergraph
    // 
    for (int v = 0; v < graph.vertices.size(); v++)
    {
        std::vector<unsigned int> hgFP;
        graph.vertices[v].data.GetFingerprint(hgFP);

std::cerr << "Hypergraph: " << std::endl;
for (std::vector<unsigned int>::const_iterator it = hgFP.begin(); it != hgFP.end(); it++)
{
    std::cerr << *it << "|";
}
std::cerr << std::endl;


        double tanimoto = OpenBabel::OBFingerprint::Tanimoto(validationFP, hgFP);

        std::cerr << "Tanimoto: " << tanimoto << std::endl;

        if (tanimoto > (double)Options::TANIMOTO) return true;
    }

    return false;
}


//
// Validate all of the molecules.
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
    // Read all of the OBMol objects using Open Babel
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

    Validate(molecules);
/*
    // Destroy all the newly-created molecules we validated.
    for (std::vector<OpenBabel::OBMol>::iterator it = molecules.begin(); it != molecules.end(); it++)
    {
        delete &(*it);
    }
*/
}


