#include <string>
#include <vector>
#include <iostream>
#include <cstring>
#include <cstdlib>


#include "Options.h"


double Options::TANIMOTO;


Options::Options(int argCount, char** vals) : argc(argCount), argv(vals)
{
    // Default values
    outFile = "molecules.sdf";
    validationFile = "";

    Options::TANIMOTO = 0.95;
}

bool Options::parseCommandLine()
{
    for (int i = 1; i < argc; i++)
    {
        if (argv[i][0] == '-')
        {
            if (!handleOption(i)) return false;
        }
        else
        {
            inFiles.push_back(argv[i]);
        }
    }

    return true;
}

//
// Deal with the actual options specified on the command-line.
//
bool Options::handleOption(int& index)
{
    //
    // Check for an error: last command-line argument is an option (with no subsequent file)
    //
    if (index + 1 >= argc)
    {
        std::cerr << "Specified option " << argv[index]
                  << " not followed by a file name." << std::endl;
        return false;
    }

    if (strcmp(argv[index], "-o") == 0)
    {
        outFile = argv[++index];
    }
    if (strcmp(argv[index], "-v") == 0)
    {
        validationFile = argv[++index];
    }
    if (strncmp(argv[index], "-tc", 3) == 0)
    {
        // not directly following; e.g. -tc 0.95
        if (strcmp(argv[index], "-tc") == 0)
        {
            Options::TANIMOTO = atof(argv[++index]);
        }
        // -tc followed by a double value; e.g. -tc0.95
        else
        {
            Options::TANIMOTO = atof(&argv[index][3]);
        }
    }

    return true;
}

