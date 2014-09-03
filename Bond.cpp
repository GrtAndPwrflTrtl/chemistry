#include <iostream>
#include <cstring>
#include <sstream>
#include <string>
#include <vector>
#include <cstdio>

#include "Bond.h"

/*********************************************************************************************************/

Bond::Bond(int id, int origin, int target) : bondID(id),
                                             originAtomID(origin),
                                             targetAtomID(target)
{
}

/*****************************************************************************************/

std::string Bond::toString() const
{
    std::ostringstream oss;

    oss << "id = " << bondID
        << " from atom " << originAtomID
        << " to atom " << targetAtomID;

    return oss.str();
}


/*****************************************************************************************/

std::ostream& operator<< (std::ostream& os, Bond& bond)
{
    os << bond.toString();

    return os;
}

/*****************************************************************************************/
