#include <iostream>
#include <cstring>
#include <sstream>
#include <string>
#include <vector>
#include <cstdio>


#include "Atom.h"


/*********************************************************************************************************/

Atom::~Atom()
{
    // delete &atomType;
}

/*********************************************************************************************************/

Atom::Atom(int id, AtomT& aType)
{
    this->atomID = id;
    this->atomType = aType;
    this->maxConnect = 0;
    this->canConnectToAnyAtom = false;
    this->ownerType = COMPLEX; // This should not be allowed.
}

/*********************************************************************************************************/

Atom::Atom(int id, bool canConnectToAnyAtom)
{
    this->atomID = id;
    this->maxConnect = 0;
    this->canConnectToAnyAtom = canConnectToAnyAtom;
}

/*********************************************************************************************************/

Atom::Atom(int id, int maxConn, std::string connType, bool canConnectToAnyAtom)
{
    this->atomID = id;
    this->maxConnect = maxConn;
    this->canConnectToAnyAtom = canConnectToAnyAtom;
}

/*********************************************************************************************************/

void Atom::SetBasedOn(const Atom& that)
{
    this->atomType = that.atomType;
    this->canConnectToAnyAtom = that.canConnectToAnyAtom;
    this->maxConnect = that.maxConnect;
    this->ownerType = that.ownerType;

    for (int a = 0; a < that.allowableTypes.size(); a++)
    {
        this->allowableTypes.push_back(that.allowableTypes[a]);
    }

/* Open Babel handles all bonds.
    for (int c = 0; c < that.connections.size(); c++)
    {
        this->connections.push_back(that.connections[c]);
    }
*/

    for (int c = 0; c < that.extConnections.size(); c++)
    {
        this->extConnections.push_back(that.extConnections[c]);
    }
}


/*********************************************************************************************************/

void Atom::addConnection(int atomIndex)
{
    this->connections.push_back(atomIndex);
}

/*********************************************************************************************************/

void Atom::addExternalConnection(int atomIndex)
{
// std::cerr << "Adding External Index: " << atomIndex << std::endl;

    this->extConnections.push_back(atomIndex);
}


/*********************************************************************************************************/

void Atom::addConnectionType(const AtomT& aType)
{
   if (std::find(allowableTypes.begin(), allowableTypes.end(), aType) == allowableTypes.end())
   {
       allowableTypes.push_back(aType);
   }
   else
   {
       std::cerr << "Duplicate allowable bond-type: " << aType << std::endl;
   }
}

/*********************************************************************************************************/

bool Atom::CanConnectTo(const Atom& that) const
{
    //
    // Disallow Linker-Linker connections.
    //
    if (this->ownerType == LINKER && that.ownerType == LINKER) return false;

    //
    // Are there any allowable spots in the atoms to connect?
    //
    if(!this->SpaceToConnect()) return false;

    if(!that.SpaceToConnect()) return false;

//std::cerr << "\tSpace avail" << std::endl;

    //
    // Does this atom allow the connection to that?
    //
    bool canConnect = false;
    for (int t = 0; t < this->allowableTypes.size(); t++)
    {
        if (this->allowableTypes[t] == that.getAtomType())
        {
            canConnect = true;
            break;
        }
    }

    // If there are no allowable connections, check to see if it accepts all connections.
    if (!canConnect)
    {
        if (!this->CanConnectToAny()) return false;
    }

    //
    // Does that atom allow the connection to this?
    //
    canConnect = false;
    for (int t = 0; t < that.allowableTypes.size(); t++)
    {
        if (that.allowableTypes[t] == this->atomType)
        {
            canConnect = true;
            break;
        }
    }

    // If there are no allowable connections, check to see if it accepts all connections.    
    if (canConnect) return true;

    return that.CanConnectToAny();
}

/*********************************************************************************************************/

bool Atom::operator==(const Atom& that) const
{
//std::cerr << "This: " << this->toString() << std::endl;
//std::cerr << "That: " << that.toString() << std::endl;

    //
    // Atom Types must match
    //
    if (atomType != that.atomType) return false;

//std::cerr << "(a)" << std::endl;

    if (this->canConnectToAnyAtom != that.canConnectToAnyAtom) return false;

//std::cerr << "(b)" << std::endl;
    if (this->maxConnect != that.maxConnect) return false;

//std::cerr << "(c)" << std::endl;
    //
    // Check the allowable connection types, connections, and external connections
    //
    if (this->allowableTypes.size() != that.allowableTypes.size()) return false;

    if (this->connections.size() != that.connections.size()) return false;

    if (this->extConnections.size() != that.extConnections.size()) return false;

//std::cerr << "(d)" << std::endl;
    for (int t = 0; t < this->allowableTypes.size(); t++)
    {
        if (std::find(that.allowableTypes.begin(),
                      that.allowableTypes.end(), this->allowableTypes[t]) == that.allowableTypes.end())
        {
            return false;
        }
    }

//std::cerr << "(e)" << std::endl;
    for (int c = 0; c < this->connections.size(); c++)
    {
        if (std::find(that.connections.begin(),
                      that.connections.end(), this->connections[c]) == that.connections.end())
        {
            return false;
        }
    }

//std::cerr << "(f)" << std::endl;
    for (int c = 0; c < this->extConnections.size(); c++)
    {
        if (std::find(that.extConnections.begin(),
                      that.extConnections.end(), this->extConnections[c]) == that.extConnections.end())
        {
            return false;
        }
    }

//std::cerr << "(g)" << std::endl;
    return true;
}


/*********************************************************************************************************/

std::string Atom::toString() const
{
    std::ostringstream oss;

    oss << atomID << ": " << atomType.toString();
    oss << " Connections{ Max: " << maxConnect;

    oss << " Allow: ";

    if (CanConnectToAny()) oss << " All Conn";
    else
    {
        for (int a = 0; a < allowableTypes.size(); a++)
        {
            oss << allowableTypes[a].toString() << " ";
        }
    }   
    oss << "\tBonds: {";

    for (int c = 0; c < connections.size(); c++)
    {
        oss << connections[c] << " ";
    }
 
    oss << "}\t\tExt Bonds: { ";

    for (int c = 0; c < extConnections.size(); c++)
    {
        oss << extConnections[c] << " ";
    }

    oss << "} }";

    return oss.str();
}

/*********************************************************************************************************/

std::ostream& operator<< (std::ostream& os, Atom& atom)
{
    os << atom.toString();

    return os;
}

/*********************************************************************************************************/
