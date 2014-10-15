#include <vector>
#include <iostream>
#include <string>
#include <map>
#include <sstream>


#include <openbabel/mol.h>
#include <openbabel/bond.h>


#include "Linker.h"
#include "EdgeAggregator.h"
#include "AtomT.h"


Linker::Linker(OpenBabel::OBMol* obmol, const std::string& name) : uniqueFragmentID(-1),
                                                                   Molecule(obmol, name, LINKER)
{
    //
    // Acquire the comment data, make a copy, parse that comment.
    //
    OpenBabel::OBCommentData* comment = static_cast<OpenBabel::OBCommentData*>(obmol->GetData("Comment"));

    std::string commentStr = comment->GetData();

    parseAppendix(commentStr);
}

//
// Parse suffix to add max connection for each atom.
//
void Linker::parseAppendix(std::string& suffix)
{
    // std::cerr << "Linker::parseAppendix: " << suffix << "|" << std::endl;

    // Use a string stream instead of manipulatiing the string
    std::stringstream suffStream(suffix);

    //
    // Read until we get "> <"
    //
    std::string line = ""; 
    while(line.find("> <") == std::string::npos)
    {
        getline(suffStream, line);
    }

    //
    // Now, read the MAX Connections for each atom.
    //
    int maxConnections = -1;
    std::string atomType;

    for(int x = 0; x < this->getNumberOfAtoms(); x++)
    {
        suffStream >> maxConnections;
        suffStream >> atomType;

        // A linker can link to any atom.
        this->atoms[x].setCanConnectToAnyAtom();
        this->atoms[x].setMaxConnect(maxConnections);
        this->atoms[x].setAtomType(* new AtomT(atomType));
        this->atoms[x].setOwnerMolecule(this);
        this->atoms[x].setOwnerMoleculeType(LINKER);
    }
}
