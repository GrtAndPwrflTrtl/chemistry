#include <cstring>
#include <vector>


#include<openbabel/descriptor.h>
#include<openbabel/fingerprint.h>


#include "Molecule.h"
#include "Bond.h"
#include "Atom.h"
#include "obgen.h"
#include "Thread_Pool.h"


#include "EdgeAggregator.h"
#include "EdgeAnnotation.h"
#include "IdFactory.h"
#include "Utilities.h"
#include "Constants.h"
#include "Options.h"



// Static allocation of the thread pool.
//Thread_Pool<OpenBabel::OBMol*, bool> Molecule::pool(THREAD_POOL_SIZE, OBGen::obgen);


Molecule::Molecule() {}
Molecule::~Molecule()
{
/*
    delete &obmol;

    for (int a = 0; a < this->atoms.size(); a++)
    {
        delete &this->atoms[a];
    }

    for (int b = 0; b < this->bonds.size(); b++)
    {
        delete &this->bonds[b];
    }
*/
}

Molecule::Molecule(OpenBabel::OBMol* mol, const std::string& n, MoleculeT t) : graphID(-1),
                                                                               moleculeID(-1),
                                                                               lipinski(true),
                                                                               obmol(mol),
                                                                               name(n),
                                                                               type(t) 
{
    // Create the initial atom / bond data based on obmol.
    localizeOBMol();

    //
    // Acquire the fingerprint of the new molecule.
    //
    OpenBabel::OBFingerprint* fpType = OpenBabel::OBFingerprint::FindFingerprint("");
    fpType->GetFingerprint(mol, this->fingerprint);
}

void Molecule::localizeOBMol()
{
    int numOfAtoms = this->obmol->NumAtoms();
    int numOfBonds = this->obmol->NumBonds();

    //
    // Translate the OB atoms into our local atoms.
    //
    for(int x = 0; x < numOfAtoms; x++)
    {
        this->addAtom(*new Atom(atomIdMaker.getNextId(), *new AtomT()));
    }

    //
    // Translate the OB Bonds into our local bonds.
    //
    for (int x = bondIdMaker.min(); x < numOfBonds; x++)
    {
        OpenBabel::OBBond* oneObBond = this->obmol->GetBondById(x);

        this->addBond((int)oneObBond->GetBeginAtom()->GetId(), (int)oneObBond->GetEndAtom()->GetId());
    }

//std::cerr << *this << std::endl;
}

bool Molecule::operator==(const Molecule& that) const
{
//std::cerr << "Checking: == " << this->getGraphID() << ":" << that.getGraphID() << std::endl;

    //
    // The foremost comparitor is the Tanimoto value
    //
    std::vector<unsigned int> thisFP = fingerprint;
    std::vector<unsigned int> thatFP;
    that.GetFingerprint(thatFP);
    double tanimoto = OpenBabel::OBFingerprint::Tanimoto(thisFP, thatFP);

    std::cerr << "Tanimoto: " << tanimoto << std::endl;

    if (tanimoto > Options::TANIMOTO) return true;

 

    //
    // Check type as well as sizes of atoms, bonds, linkers, and rigids.
    //
    if (this->type != that.type) return false;

//std::cerr << "1" << std::endl;

    if (this->atoms.size() != that.atoms.size()) return false;

//std::cerr << "2" << std::endl;
    if (this->bonds.size() != that.bonds.size()) return false;

//std::cerr << "3" << std::endl;
    if (this->rigids.size() != that.rigids.size()) return false;

//std::cerr << "4" << std::endl;
    if (this->linkers.size() != that.linkers.size()) return false;

//std::cerr << "5" << std::endl;
    //
    // Check the contents of the atoms, bonds, rigids, and linkers
    //
    for (int a = 0; a < this->atoms.size(); a++)
    {
        if (find(that.atoms.begin(), that.atoms.end(), this->atoms[a]) == that.atoms.end())
        {
//std::cerr << "Not Found: " << this->atoms[a].toString() << endl;

            return false;
        }
    }

//std::cerr << "6" << std::endl;
    for (int b = 0; b < this->bonds.size(); b++)
    {
        if (find(that.bonds.begin(), that.bonds.end(), this->bonds[b]) == that.bonds.end()) return false;
    }

//std::cerr << "7" << std::endl;
    for (int r = 0; r < this->rigids.size(); r++)
    {
        if (find(that.rigids.begin(), that.rigids.end(), this->rigids[r]) == that.rigids.end()) return false;
    }

//std::cerr << "8" << std::endl;
    for (int ell = 0; ell < this->linkers.size(); ell++)
    {
        if (find(that.linkers.begin(), that.linkers.end(), this->linkers[ell]) == that.linkers.end()) return false;
    }

//std::cerr << "9" << std::endl;
    return true;
}

// *****************************************************************************

//
// Exhaustively create any possible compositions between two Molecules
//
/*
    1. There should be no linker-to-linker connections.
       The only acceptable connections are rigid-to-rigid, and linker-to-rigid.

    2. The only loops that we should allow are loops that are inherent to the original linkers and rigids.
       In other words, if an original linker or rigid has a ring of atoms, then that loop is ok, but if you
       were to draw the molecule as a graph where each node is not an atom, but either a linker or a rigid,
       then that graph should not have any loops in it.

    3. The molecule should obey all 4 of Lipinski's rules, also known as the rule of 5.
       Here is a summary:
           a. A molecule should have a mass of less than 500 daltons, which can be calculated with
              the obmol function GetExactMass().
           b. No more than 5 hydrogen bond doners.
           c. No more than 10 hydrogen bond acceptors.
           d. An octanol-water partition coefficient log P not greater than 5.
*/

//
// (a) Check if the molecular weight is too heavy.
//
bool Molecule::exceedsMaxMolecularMass()
{
    return this->obmol->GetMolWt() > MAX_DALTON_WEIGHT;
}

bool Molecule::isLipinskiCompliant()
{
    //
    // (b) Hydrogen Bond donors
    //
    OpenBabel::OBDescriptor* pDescr = OpenBabel::OBDescriptor::FindType("HBD");
    if(pDescr)
    {
        if(pDescr->Predict(this->obmol) > 5.0) return false; 
    }
    else
    {
        std::cerr << "HBD not found." << std::endl;
    }

    //
    // (c) Hydrogen Bond Acceptors
    //
    pDescr = OpenBabel::OBDescriptor::FindType("HBA1");
    if(pDescr)
    {
        if (pDescr->Predict(this->obmol) > 10.0) return false;
    }
    else
    {
        std::cerr << "HBA1 not found." << std::endl;
    }

    //
    // Octanol-water partition coefficient log P not greater than 5
    //
    pDescr = OpenBabel::OBDescriptor::FindType("logP");
    if(pDescr)
    {
        if (pDescr->Predict(this->obmol) > 5.0) return false;
    }
    else
    {
        std::cerr << "logP not found." << std::endl;
    }

    return true;
}

//
// The method by which we create new molecules is such that loops with linkers and rigids will not occur.
//
bool Molecule::ContainsLoops() const
{
    return false;
}

//
// Overall molecular criteria satisfaction.
//
bool Molecule::satisfiesMoleculeSynthesisCriteria()
{
    //
    // The only eliminating criteria is for molecule mass to be too large.
    //
    if (exceedsMaxMolecularMass()) return false;

    // Do the linkers / rigids create a loop in the molecule?
    if (ContainsLoops()) return false;

    // Check Lipinski compliance
    lipinski = isLipinskiCompliant();

    return true;
}

std::vector<EdgeAggregator*>* Molecule::Compose(const Molecule& that) const
{
    std::vector<EdgeAggregator*>* newMolecules = new std::vector<EdgeAggregator*>();

    //
    // For each atom in this molecule, does it connect to an atom in that molecule?
    //
    for (unsigned int thisA = 0; thisA < atoms.size(); thisA++)
    {
        for (unsigned int thatA = 0; thatA < that.atoms.size(); thatA++)
        {

//std::cerr << thisA << ":" << thatA << std::endl;

            //
            // We've established the fact that these two particular atoms are connectable
            // Can we actually connect these two molecules at these two atoms?
            //
            if (atoms[thisA].CanConnectTo(that.atoms[thatA]))
            {

std::cerr << "Connection Possible: " << std::endl;
std::cerr << "\t" << atoms[thisA].toString() << std::endl;
std::cerr << "\t" << that.atoms[thatA].toString() << std::endl;

                // Create a new molecule;
                // The indices are the new indices when the atoms and bonds are combined together.
                Molecule* newMol = ComposeToNewMolecule(that, thisA + 1, thatA + this->atoms.size() + 1);

                // Check if the new Molecule satisfies all criteria.
                if (newMol->satisfiesMoleculeSynthesisCriteria())
                {
// std::cerr << "Criteria Satisfied." << std::endl;

// std::cerr << "Created: " << *newMol << std::endl;

                    // Antecedent
                    std::vector<Molecule> ante;
                    ante.push_back(*this);
                    ante.push_back(that);

                    // Add the new molecule / edge to the list of new molecules
                    newMolecules->push_back(new EdgeAggregator(ante, newMol, new EdgeAnnotationT()));

                    // Only if lipinski compliant do we run a molecule through obgen.
                    // To save time, spawn a thread.
                    if (newMol->lipinski)
                    {
                        // Robert, this is the exact call we need to place into a thread.
                        //pool.push(newMol->obmol);
                        //OBGen::obgen(newMol->obmol);
                    }
                }
            }
        }
    } 

    return newMolecules;
}

// *****************************************************************************
//
// Create the local informations:
//    (1) Using the same style invoked by OpenBabel, we modify the indices by adding their
//        number of atoms in this to the index of the atoms in that.
//    (2) Bonds based on indices will be updated accordingly.
//
Molecule* Molecule::ComposeToNewMolecule(const Molecule& that,
                                         int thisAtomIndex,
                                         int thatAtomIndex) const
{
    //
    // Combine the Open Babel representations.
    //
    OpenBabel::OBMol* newOBMol = new OpenBabel::OBMol(*this->obmol);

    *newOBMol += *that.obmol;

    // Add the new Open Babel bond.
    // This, along with the construction of the local molecule takes care of the new bond
    newOBMol->AddBond(thisAtomIndex, thatAtomIndex, 1); // Order of the bond is 1.

    // Remove the comment information as it is no longer relevant to this molecule.
    newOBMol->DeleteData("Comment");

    //
    //
    // Transfer the local data.
    //
    //
    // Create the new Molecule object; it will create the localized information.
    Molecule* newLocal = new Molecule(newOBMol, "complex", COMPLEX);

    //
    // Copy the local atom information
    //
    int newAtomCount = 0;
    for (int a = 0; a < this->atoms.size(); a++, newAtomCount++)
    {
        newLocal->atoms[newAtomCount].SetBasedOn(this->atoms[a]);
    }

    for (int a = 0; a < that.atoms.size(); a++, newAtomCount++)
    {
        newLocal->atoms[newAtomCount].SetBasedOn(that.atoms[a]);
    }

    //
    // Combine all the linkers and rigids
    //
    for (int r = 0; r < this->rigids.size(); r++)
    {
        newLocal->rigids.push_back(this->rigids[r]);
    }

    for (int r = 0; r < that.rigids.size(); r++)
    {
        newLocal->rigids.push_back(that.rigids[r]);
    }

    for (int ell = 0; ell < this->linkers.size(); ell++)
    {
        newLocal->linkers.push_back(this->linkers[ell]);
    }

    for (int ell = 0; ell < that.linkers.size(); ell++)
    {
        newLocal->linkers.push_back(that.linkers[ell]);
    }

    // Add local information to the new molecule.
    // Bonds in open babel start indexing at 1.
    newLocal->atoms[thisAtomIndex-1].addExternalConnection(thatAtomIndex-1);
    newLocal->atoms[thatAtomIndex-1].addExternalConnection(thisAtomIndex-1);

    return newLocal;
}

// *****************************************************************************

//
// On-Demand acquisition of the fingerprint one time.
//
void Molecule::GetFingerprint(std::vector<unsigned int>& returnFP) const
{
    if (this->fingerprint.empty())
    {
        std::cerr << "Expected pre-computed fingerprint." << std::endl;
    }

    returnFP = this->fingerprint;
}


// *****************************************************************************

int Molecule::getAtomIndex(int id) const
{
    for(int x = 0; x < this->atoms.size(); x++)
    {
        if(this->atoms[x].getAtomID() == id)
        {
            return x;
        }
    }
    
    return -1;
}

// *****************************************************************************

Atom Molecule::getAtom(int id) const
{
    int index = getAtomIndex(id);
    
    if (index == -1) throw "Bond id not found";
    
    return atoms[index];
}

// *****************************************************************************

bool Molecule::addBond(int xID, int yID) // , eTypeOfBondT bt, eStatusBitT s)
{
    int xIndex = getAtomIndex(xID);
    int yIndex = getAtomIndex(yID);

    if (xIndex == -1 || yIndex == -1) return false;

    atoms[xIndex].addConnection(yIndex);
    atoms[yIndex].addConnection(xIndex);

    this->bonds.push_back(Bond(this->bonds.size(), xID, yID));

    return true;
}

// *****************************************************************************

int Molecule::getBondIndex(int id) const
{
    for(int x = 0; x < this->bonds.size(); x++)
    {
        if(this->bonds[x].getBondID() == id)
        {
            return x;
        }
    }
    
    return -1;
}

// *****************************************************************************

int Molecule::getBondIndex(int xID, int yID) const
{
    for(int x = 0; x < this->bonds.size(); x++)
    {
        if(this->bonds[x].getOriginAtomID() == xID && this->bonds[x].getTargetAtomID() == yID)
        {
            return x;
        }
         
        if(this->bonds[x].getOriginAtomID() == yID && this->bonds[x].getTargetAtomID() == xID)
        {
            return x;
        }
    }

    return -1;
}

// *****************************************************************************

Bond Molecule::getBond(int id) const
{
    int index = getBondIndex(id);
    
    if (index == -1) throw "Bond id not found";
    
    return bonds[index];
}

// *****************************************************************************

Bond Molecule::getBond(int xID, int yID) const
{
    int index = getBondIndex(xID, yID);
    
    if (index == -1) throw "Bond id not found";
    
    return bonds[index];
}

// *****************************************************************************

void Molecule::addAtom(const Atom& a)
{
    this->atoms.push_back(a);

    // TODO: Add atom to this->obmol
}


// *****************************************************************************

std::string Molecule::toString() const
{
    std::ostringstream oss;

    oss << "Molecule: " << graphID << " ";

    if (IsLinker()) oss << " is a linker.";
    else if (IsRigid())  oss << " is a rigid.";
    else if (IsComplex())
    {
        oss << "There are ";
        oss << getNumberOfRigids() << " rigids, and ";
        oss << getNumberOfLinkers();
        oss << " linkers." << std::endl;
    }

    oss << "There are ";
    oss << getNumberOfAtoms();
    oss << " atoms, and ";
    oss << getNumberOfBonds();
    oss << " bonds." << std::endl;

    if (getNumberOfAtoms() > 0)
    {
        oss << "Atoms:" << std::endl;

        for(int x = 0; x < getNumberOfAtoms(); x++)
        {
            oss << "\t" << getAtom(x).toString() << std::endl;
        }
    }

    if (getNumberOfBonds() > 0)
    {
        oss << "Bonds:" << std::endl;

        for (int x = 0; x < getNumberOfBonds(); x++)
        {
            oss << "\t" << getBond(x).toString() << std::endl;
        }
    }

    return oss.str();
}


// *****************************************************************************

std::ostream& operator<< (std::ostream& os, Molecule& mol)
{
    os << mol.toString() << std::endl;

    return os;
}
