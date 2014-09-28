#include <cstring>
#include <string>
#include <vector>
#include <sstream>
#include <cctype>


#include "Rigid.h"
#include "EdgeAggregator.h"
#include "Utilities.h"


Rigid::Rigid(OpenBabel::OBMol* obmol, const std::string& name) : Molecule(obmol, name, RIGID)
{
    //
    // Acquire the comment data, make a copy, parse that comment.
    //
    OpenBabel::OBCommentData* comment = static_cast<OpenBabel::OBCommentData*>(obmol->GetData("Comment"));

    std::string commentStr = comment->GetData();

    parseAppendix(commentStr);

    this->rigids.push_back(this);
}

void Rigid::parseAppendix(std::string& suffix)
{
    // Use a string stream instead of manipulatiing the string
    std::stringstream suffStream(suffix);

    //
    // Read until we get "> <"
    //
    std::string line = "";
    while(line.find("> <") == std::string::npos)
    {
        getline(suffStream, line);
        std::cout << "---------" << line << "**********" << std::endl;
    }

    //
    // Read the Atom Types
    //
    std::string atomType;
    for(int x = 0; x < this->getNumberOfAtoms(); x++)
    {
        suffStream >> atomType;

        this->atoms[x].setAtomType(* new AtomT(atomType));
        this->atoms[x].setOwnerMoleculeType(RIGID);
    }

    // Get the next line.
    getline(suffStream, line);

    //
    // Read until we get "> <"
    //
    while(line.find("> <") == std::string::npos)
    {
        getline(suffStream, line);
    }

    //
    // Read Branches
    //
    int atomId = -1;

    while (!isspace(suffStream.peek()))    
    {
        suffStream >> atomId;

        while (suffStream.peek() != '\n' && suffStream.peek() != '\r')
        {
            suffStream >> atomType;

            this->atoms[atomId - 1].setMaxConnect(1);
            this->atoms[atomId - 1].addConnectionType(* new AtomT(atomType));

            eatWhiteToNewLineOrChar(suffStream);
        }

        // Get the newline
        suffStream.get();
    }

    //
    // Read through the $$$$
    //
    while(line.find("$$$$") == std::string::npos)
    {
        getline(suffStream, line);
    }
}
