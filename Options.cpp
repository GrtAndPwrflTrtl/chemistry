#include <string>
#include <vector>
#include <iostream>
#include <cstring>
#include <cstdlib>


#include "Options.h"
#include "Constants.h"


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

    if (strncmp(argv[index], "-mw", 3) == 0)
    {
        if (strcmp(argv[index], "-mw") == 0)
            MOLWT_UPPERBOUND = atof(argv[++index]);
        else
            MOLWT_UPPERBOUND = atof(&argv[index][3]);
    }
    if (strncmp(argv[index], "-hd", 3) == 0)
    {
        if (strcmp(argv[index], "-hd") == 0)
            HBD_UPPERBOUND = atof(argv[++index]);
        else
            HBD_UPPERBOUND = atof(&argv[index][3]);
    }
    if (strncmp(argv[index], "-ha", 3) == 0)
    {
        if (strcmp(argv[index], "-ha") == 0)
            HBA1_UPPERBOUND = atof(argv[++index]);
        else
            HBA1_UPPERBOUND = atof(&argv[index][3]);
    }
    if (strncmp(argv[index], "-lp", 3) == 0)
    {
        if (strcmp(argv[index], "-lp") == 0)
            LOGP_UPPERBOUND = atof(argv[++index]);
        else
            LOGP_UPPERBOUND = atof(&argv[index][3]);
    }
    if (strncmp(argv[index], "-hl", 3) == 0)
    {
        if (strcmp(argv[index], "-hl") == 0)
            HIERARCHICAL_LEVEL_BOUND = atoi(argv[++index]);
        else
            HIERARCHICAL_LEVEL_BOUND = atoi(&argv[index][3]);
    }



    return true;
}

