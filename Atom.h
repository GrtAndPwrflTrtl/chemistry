#ifndef _ATOM_GUARD
#define _ATOM_GUARD 1

#include <string>
#include <vector>
#include <algorithm>


#include "AtomT.h"


typedef enum
{
    RIGID,
    LINKER,
    COMPLEX
} MoleculeT;


class Atom
{
  private:
    int atomID;

    // This atom's type
    AtomT atomType;

    // The classification of the original linker / rigid that contained this atom
    MoleculeT ownerType;     

    // In the case of a linker atom, we can connect to anything
    bool canConnectToAnyAtom;

    // The maximum number of connections to this atom allowable.
    int maxConnect;

    // std::vector<int> bonds;
    
    // Allowable atom types for connections
    std::vector<AtomT> allowableTypes;

    // The actual (index) connections within the existent linker or rigid.
    std::vector<int> connections;

    // The actual (index) connections within the existent linker or rigid.
    std::vector<int> extConnections;


  public:
    //
    // Get functions
    //
    int getAtomID() const { return this->atomID; }
    AtomT getAtomType() const { return this->atomType; }
    int getMaxConnect() const { return this->maxConnect; }

    bool CanConnectToAny() const { return canConnectToAnyAtom; }
    bool SpaceToConnect() const
    {
//        std::cerr << "\t" << this->maxConnect << " " << extConnections.size() << std::endl;

        return this->maxConnect > extConnections.size();
    }
    bool CanConnectTo(const Atom& that) const;

    //
    // Set Functions
    //
    void setAtomID(int x) { this->atomID = x; }
    void setAtomType(AtomT x) { this->atomType = x; }

    void setCanConnectToAnyAtom() { this->canConnectToAnyAtom = true; }
    void setMaxConnect(int x) { this->maxConnect = x; }

    void SetBasedOn(const Atom& atom);
    void setOwnerMoleculeType(MoleculeT t) { ownerType = t; }

    Atom(int id, AtomT& atomType);
    Atom(int id, bool canConnectToAnyAtom);
    Atom(int id, int maxConn, std::string connType, bool canConnectToAnyAtom);
    ~Atom();

    void addConnectionType(const AtomT& aType);
    void addConnection(int);
    void addExternalConnection(int);

    std::string toString() const;
    friend std::ostream& operator<< (std::ostream& os, Atom& atom);
    bool operator==(const Atom& that) const;
};

#endif
