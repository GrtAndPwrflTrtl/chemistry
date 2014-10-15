#ifndef _OPTIONS_GUARD
#define _OPTIONS_GUARD 1

#include <string>
#include <vector>

//
// A aggregation class for options specified in the command-line
//
class Options
{
  public:
    Options(int, char**);
    bool parseCommandLine();

    std::string outFile;
    std::string validationFile;
    std::vector<std::string> inFiles;

    static double TANIMOTO;
    static bool THREADED;
    static unsigned int OBGEN_THREAD_POOL_SIZE;

  private:
    int argc;
    char** argv;

    bool handleOption(int& index);
};

#endif
