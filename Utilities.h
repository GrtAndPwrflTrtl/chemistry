#ifndef _UTILITIES_GUARD
#define _UTILITIES_GUARD 1

#include <cmath>
#include <vector>
#include <iostream>
#include <string>

using namespace std;

// log base 2 of the input value.
double log2(double value);

// How many bits in binary to represent this number?
unsigned int numBinaryBits(unsigned int value);

string MakeString(const char[], int);
std::string MakeString(const char s1[], const char s2[]);
std::string MakeString(const char s1[], std::string s2);

template<class T>
bool Contains(std::vector<T> list, const T& val)
{
    for (int i = 0; i < list.size(); i++)
    {
        if (list[i].Equals(val)) return true;
    }
    return false;
}

//
// To simulate instanceof in Java.
//
template<typename CheckType, typename InstanceType>
bool instanceOf(InstanceType &Instance)
{
  return dynamic_cast<CheckType*>(&Instance) != NULL;
}

void eatWhiteToNewLineOrChar(std::istream&);
void eatWhiteLines(std::istream&);

#endif
