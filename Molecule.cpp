#include <cstring>
#include <vector>
#include <bitset>
#include <utility>


#include<openbabel/descriptor.h>
#include<openbabel/fingerprint.h>


#include "Molecule.h"
#include "Bond.h"
#include "Atom.h"
#include "obgen.h"
#include "Thread_Pool.h"
#include "Rigid.h"
#include "Linker.h"
#include <pthread.h>


#include "EdgeAggregator.h"
#include "EdgeAnnotation.h"
#include "IdFactory.h"
#include "Utilities.h"
#include "Constants.h"
#include "Options.h"
#include "FragmentGraph.h"


// Static allocation of the thread pool.
//Thread_Pool<OpenBabel::OBMol*, bool> Molecule::pool(THREAD_POOL_SIZE, OBGen::obgen);

// a static class member needs to be defined in some .cpp
pthread_mutex_t Molecule::openbabel_lock;

unsigned int Molecule::RIGID_INDEX_START = -1;
unsigned int Molecule::RIGID_INDEX_END = -1;
unsigned int Molecule::LINKER_INDEX_START = -1;
unsigned int Molecule::LINKER_INDEX_END = -1;
unsigned int Molecule::FRAGMENT_END_INDEX = -1;
unsigned int Molecule::NUM_UNIQUE_FRAGMENTS = -1;

std::vector<Molecule*> Molecule::baseMolecules;
IdFactory Molecule::connectionIdMaker(100);
static const unsigned int NO_CONNECTION = -1;



Molecule::Molecule() : lipinskiPredicted(false),
                       lipinskiEstimated(false)
{
    init_openbabel_lock();

    numLinkers = 0;
    numRigids = 0;
}

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

Molecule::Molecule(OpenBabel::OBMol* mol, const std::string& n, MoleculeT t) :
    uniqueIndexID(-1),
    numLinkers(-1),
    numRigids(-1),
    numUniqueLinkers(-1),
    numUniqueRigids(-1),
    obmol(mol),
    name(n),
    type(t),
    fragmentCounter(0),
    lipinskiPredicted(false),
    lipinskiEstimated(false) 
{

    init_openbabel_lock();

    // Locking open babel since it is not thread-safe (at all)
    pthread_mutex_lock(&Molecule::openbabel_lock);

    // Create the initial atom / bond data based on obmol.
    localizeOBMol();

    // Unlocking open babel
    pthread_mutex_unlock(&openbabel_lock);
}

void Molecule::init_openbabel_lock()
{
    //
    // initializing global openbabel lock (once)
    //
    static bool openbabel_lock_init = false;
    if (openbabel_lock_init)
    {
        pthread_mutex_init(&openbabel_lock, NULL);
        openbabel_lock_init = true;
    }
}

void Molecule::initFragmentDevices()
{
    initFragmentInfo();
    calcFragmentInfo();


    // Indicate we are using this fragment
    fragmentCounter[uniqueIndexID] = 1;

    //
    // Create the connection identifiers for this linker / rigid
    //
    // Find the connections for this molecule and create ids for them.
    // These ids are unique to the linker and rigid.
    for (int a = 0; a < atoms.size(); a++)
    {
        unsigned int id = atoms[a].getMaxConnect() > 0 ? connectionIdMaker.getNextId() 
                                                       : NO_CONNECTION;
        connectionIDs.push_back(id);
        atoms[a].setConnectionID(id);
    }
}

//
// Calculate the number of linkers / rigids (copies and unique)
//
void Molecule::initFragmentInfo()
{
    if (this->fragmentCounter == 0)
    {
        int sz = Molecule::NUM_UNIQUE_FRAGMENTS;

        // Create the reference count array
        fragmentCounter = new unsigned int[sz];

        // init the counters to zero
        memset(fragmentCounter, 0, sz * sizeof(unsigned int));
    }
}

void Molecule::calcFragmentInfo()
{
    this->numLinkers = 0;
    this->numUniqueLinkers = 0;
    this->numRigids = 0;
    this->numUniqueRigids = 0;

    for (int r = Molecule::RIGID_INDEX_START; r <= Molecule::RIGID_INDEX_END; r++)
    {
        if (fragmentCounter[r] != 0)
        {
            this->numUniqueRigids++;
            this->numRigids += fragmentCounter[r];
            this->rigids.push_back(static_cast<Rigid*>(baseMolecules[r]));
        }
    }

    for (int ell = Molecule::LINKER_INDEX_START; ell <= Molecule::LINKER_INDEX_END; ell++)
    {
        if (fragmentCounter[ell] != 0)
        {
            this->numUniqueLinkers++;
            this->numLinkers += fragmentCounter[ell];
            this->linkers.push_back(static_cast<Linker*>(baseMolecules[ell]));
        }
    }
}


void Molecule::initGraphRepresentation()
{
    // Create the graph (using the number of unique fragments)
    fingerprint = new FragmentGraph();

    // Add the new node (with sub-nodes) to the graph.
    unsigned int nodeIndex = fingerprint->AddInitialNode(this);

    // Add the node to all the connection points
    for (int a = 0; a < atoms.size(); a++)
    {
        if (connectionIDs[a] != NO_CONNECTION)
        {
            atoms[a].setGraphNodeIndex(std::make_pair(this->uniqueIndexID, nodeIndex));
        }
        // else atoms[a].setGraphNodeIndex(std::make_pair(uniqueIndexID, -1));
    }
}

void Molecule::SetBaseMoleculeInfo(const std::vector<Molecule*> baseMols,
                                  unsigned int numRigids, unsigned int numLinkers)
{
    baseMolecules = baseMols;

    // Set the base molecule indices.
    Molecule::RIGID_INDEX_START = 0;
    Molecule::RIGID_INDEX_END = numRigids - 1;
    Molecule::LINKER_INDEX_START = numRigids;
    Molecule::LINKER_INDEX_END = numRigids + numLinkers - 1;
    Molecule::FRAGMENT_END_INDEX = Molecule::LINKER_INDEX_END;
    Molecule::NUM_UNIQUE_FRAGMENTS = numRigids + numLinkers;
}

void Molecule::openBabelPredictLipinski()
{
    pthread_mutex_lock(&Molecule::openbabel_lock);

    // calculate the molecular weight, H donors and acceptors and the plogp
    OpenBabel::OBDescriptor* pDesc1 = OpenBabel::OBDescriptor::FindType("HBD");
    OpenBabel::OBDescriptor* pDesc2 = OpenBabel::OBDescriptor::FindType("HBA1");
    OpenBabel::OBDescriptor* pDesc4 = OpenBabel::OBDescriptor::FindType("logP");

    if (!pDesc1) cerr << "HBD not found" << endl;
    if (!pDesc2) cerr << "HBA1 not found" << endl;
    if (!pDesc4) cerr << "logP not found" << endl;
    if (!pDesc1 || !pDesc2 || !pDesc4) return;

    MolWt = (this->obmol)->GetMolWt(); // the standard molar mass given by IUPAC atomic masses (amu)
    HBD = pDesc1->Predict(this->obmol);
    HBA1 = pDesc2->Predict(this->obmol);
    logP = pDesc4->Predict(this->obmol);

    pthread_mutex_unlock(&Molecule::openbabel_lock);

    lipinskiPredicted = true;
    lipinskiEstimated = false;
}


//
// Near the end of the synthesis process, there is little benefit 
// to composing molecules if the two molecules will exceed the additive molecular weight.  
//
bool Molecule::willExceedMolecularWeight(const Molecule &mol1, const Molecule &mol2)
{
    return 6.6746 + 0.95965 * (mol1.getMolWt() + mol2.getMolWt()) > MOLWT_UPPERBOUND;
}

void Molecule::estimateLipinski(const Molecule &mol1, const Molecule &mol2)
{
    if (!mol1.lipinskiPredicted && !mol1.lipinskiEstimated)
    {
        cerr << "estimateLipinski: mol1 has no lipinski coefficients available" << endl;
        return;
    }
    if (!mol2.lipinskiPredicted && !mol2.lipinskiEstimated)
    {
        cerr << "estimateLipinski: mol2 has no lipinski coefficients available" << endl;
        return;
    }

    double calc_MolWt = mol1.getMolWt() + mol2.getMolWt();
    double calc_HBD = mol1.getHBD() + mol2.getHBD();
    double calc_HBA1 = mol1.getHBA1() + mol2.getHBA1();
    double calc_logP = mol1.getlogP() + mol2.getlogP();

    MolWt = 6.6746 + 0.95965 * calc_MolWt; // the molar mass given by IUPAC atomic masses (amu)
    HBD = 0.41189 + 0.4898 * calc_HBD;
    HBA1 = 0.278 + 0.93778 * calc_HBA1;
    logP = 0.84121 + 0.59105 * calc_logP;

    lipinskiPredicted = false;
    lipinskiEstimated = true;
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

        this->addBond((int)oneObBond->GetBeginAtom()->GetId(), 
                      (int)oneObBond->GetEndAtom()->GetId());
    }
}

//
// Molecular comparison is through the use of a local fingerprinting scheme.
// We construct a fingerprint by noting the connection anchors for fragments
// and constructing a graph based on those anchor points. A fingerprint equality
// check performs graph isomorphism.
//
bool Molecule::operator==(const Molecule& that) const
{
// std::cout << "Comparing: " << *this << " and " << that << std::endl;

    //
    // The fragment counter maintains the number of instances of each specific fragment;
    // if any of those counts differ, we have non-isomorphism.
    //
    for (int f = 0; f < Molecule::FRAGMENT_END_INDEX; f++)
    {
        if (this->fragmentCounter[f] != that.fragmentCounter[f])
        {
            return false;
        }
    }

    //
    // If we reach this point in the code, we can expect the two molecules to have the
    // same number of (1) linkers, (2) rigids, (3) unique rigids, (4) unique linkers,
    // (5) bonds, and (6) # atoms
    //
    //
    // Fingerprint verification.
    //
    // Fingerprint checking is last since it is slow; check other characteristics first.
    //
    return this->fingerprint->IsIsomorphicTo(that.getFingerprint());
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
bool Molecule::exceedsMaxEstimatedThresholds()
{
    if (!lipinskiPredicted && !lipinskiEstimated) this->openBabelPredictLipinski();

    if (!lipinskiPredicted && !lipinskiEstimated)
    {
        cerr << "didnt predict and failed to estimate lipinski coefficients" << endl;
        return false;
    }

    if (MolWt > MOLWT_UPPERBOUND)
    {
        std::cout << "# Fragments (" << this->size() << ")"
                  << MolWt << " > " << MOLWT_UPPERBOUND << std::endl;
    }

    // (b) Hydrogen Bond donors
    if (HBD > HBD_UPPERBOUND)
    {
        //std::cerr << "Failed To Add due to HBD" << std::endl;
        return false;
    }

    // (c) Hydrogen Bond Acceptors
    if (HBA1 > HBA1_UPPERBOUND)
    {
        //std::cerr << "Failed To Add due to HBA1" << std::endl;
        return false;
    }

    return MolWt > MOLWT_UPPERBOUND;
}

bool Molecule::isOpenBabelLipinskiCompliant()
{
    if (!lipinskiPredicted) this->openBabelPredictLipinski();

    return this->isLipinskiCompliant();
}

//
// Static function to check whether an OpenBabel Mol is Lipinski compliant.
//
// No need for locks since locks should go AROUND the function call.
bool Molecule::isOpenBabelLipinskiCompliant(OpenBabel::OBMol& mol)
{
    // calculate the molecular weight, H donors and acceptors and the plogp
    OpenBabel::OBDescriptor* pDesc1 = OpenBabel::OBDescriptor::FindType("HBD");
    OpenBabel::OBDescriptor* pDesc2 = OpenBabel::OBDescriptor::FindType("HBA1");
    OpenBabel::OBDescriptor* pDesc4 = OpenBabel::OBDescriptor::FindType("logP");

    if (!pDesc1) throw "HBD not found";
    if (!pDesc2) throw "HBA1 not found";
    if (!pDesc4) throw "logP not found";

    // (b) Hydrogen Bond donors
    if (pDesc1->Predict(&mol) > HBD_UPPERBOUND)
    {
        // std::cerr << "Failed HBD" << std::endl;
        return false;
    }

    // (c) Hydrogen Bond Acceptors
    if (pDesc2->Predict(&mol) > HBA1_UPPERBOUND)
    {
        // std::cerr << "Failed HBA1" << std::endl;
        return false;
    }

    // Octanol-water partition coefficient log P not greater than 5
    if (pDesc4->Predict(&mol) > LOGP_UPPERBOUND)
    {
        // std::cerr << "Failed PLOGP" << std::endl;
        return false;
    }

    return true;
}

bool Molecule::isLipinskiCompliant()
{
    if (!lipinskiPredicted && !lipinskiEstimated) this->openBabelPredictLipinski();

    if (!lipinskiPredicted && !lipinskiEstimated)
    {
        cerr << "didnt predict and failed to estimate lipinski coefficients" << endl;
        return false;
    }

    // (b) Hydrogen Bond donors
    if (HBD > HBD_UPPERBOUND) return false;

    // (c) Hydrogen Bond Acceptors
    if (HBA1 > HBA1_UPPERBOUND) return false;

    // Octanol-water partition coefficient log P not greater than 5
    if (logP > LOGP_UPPERBOUND) return false;

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
    if (exceedsMaxEstimatedThresholds()) return false;

    // Do the linkers / rigids create a loop in the molecule?
    if (ContainsLoops()) return false;

    return true;
}

std::vector<EdgeAggregator*>* Molecule::Compose(const Molecule& that) const
{
    std::vector<EdgeAggregator*>* newMolecules = new std::vector<EdgeAggregator*>();


    //
    // Pre-emptively check the molecular weight to see if there is a benefit
    // to composing these molecules.
    //
    if (Molecule::willExceedMolecularWeight(*this, that)) return newMolecules;

    //
    // For each atom in this molecule, does it connect to an atom in that molecule?
    //
    for (unsigned int thisA = 0; thisA < atoms.size(); thisA++)
    {
        for (unsigned int thatA = 0; thatA < that.atoms.size(); thatA++)
        {
            //
            // We've established the fact that these two particular atoms are connectable
            // Can we actually connect these two molecules at these two atoms?
            //
            if (atoms[thisA].CanConnectTo(that.atoms[thatA]))
            {

                if (g_debug_output)
                {
                    std::cerr << "Connection Possible: " << std::endl;
                    std::cerr << "\t" << atoms[thisA].toString() << std::endl;
                    std::cerr << "\t" << that.atoms[thatA].toString() << std::endl;
                }

                // Create a new molecule;
                // The indices are the new indices when the atoms and bonds are combined together.
                Molecule* newMol = ComposeToNewMolecule(that, thisA + 1,
                                                        thatA + this->atoms.size() + 1);

                // Check if the new Molecule satisfies all criteria.
                if (newMol->satisfiesMoleculeSynthesisCriteria())
                {
                    // std::cerr << "Created: " << *newMol << std::endl;

                    // Antecedent
                    std::vector<Molecule> ante;
                    ante.push_back(*this);
                    ante.push_back(that);

                    // Add the new molecule / edge to the list of new molecules
                    newMolecules->push_back(new EdgeAggregator(ante, newMol, new EdgeAnnotationT()));
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
    pthread_mutex_lock(&openbabel_lock);
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

    // Unlocking open babel; the next constructor call performs a relock
    pthread_mutex_unlock(&openbabel_lock);

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

    int firstThatIndex = newAtomCount;
    for (int a = 0; a < that.atoms.size(); a++, newAtomCount++)
    {
        newLocal->atoms[newAtomCount].SetBasedOn(that.atoms[a]);
    }

    // Init the fragment counter container.
    newLocal->initFragmentInfo();

    // Combine all the linkers and rigids into this molecule.
    for (int f = 0; f <= FRAGMENT_END_INDEX; f++)
    {
        newLocal->fragmentCounter[f] = this->fragmentCounter[f] + that.fragmentCounter[f];

    }
    // Calculate all the fragment values: summary data.
    newLocal->calcFragmentInfo();

    //
    // Add local information to the new molecule.
    // Bonds in open babel start indexing at 1.
    //
    newLocal->atoms[thisAtomIndex-1].addExternalConnection(thatAtomIndex-1);
    newLocal->atoms[thatAtomIndex-1].addExternalConnection(thisAtomIndex-1);

//std::cout << "Adding to fingerprint" << *this << that << std::endl;

    // Create the fingerprint graph for the new molecule by:
    //   (1) copying this fignerprint graph
    newLocal->fingerprint = this->fingerprint->copy();



/*
std::cerr << "Copy: ";

std::cout << *this->fingerprint << std::endl << "+++++++++++" << std::endl;
std::cout << *that.fingerprint << std::endl << "===========" << std::endl;
std::cout << *newLocal->fingerprint << std::endl;

std::cout << "Graph node index: ("
          << newLocal->atoms[thisAtomIndex - 1].getGraphNodeIndex().first
          << ", " << newLocal->atoms[thisAtomIndex - 1].getGraphNodeIndex().second
          << ")" << std::endl;
*/

    // Add the new linker / rigid connection to the graph
    std::pair<unsigned int, unsigned int> toIndex;
    toIndex = newLocal->fingerprint->AddEdgeAndNode(
                          newLocal->atoms[thisAtomIndex - 1].getConnectionID(),
                          newLocal->atoms[thisAtomIndex - 1].getGraphNodeIndex(),
                          newLocal->atoms[thatAtomIndex - 1],
                          that);

    // Update the atoms of the new 'to' node to reflect the proper indices in the graph.
    for (int a = 0; a < that.atoms.size(); a++)
    {
        newLocal->atoms[firstThatIndex++].UpdateIndices(toIndex);
    }

/*
std::cout << *this->fingerprint << std::endl << "+++++++++++" << std::endl;
std::cout << *that.fingerprint << std::endl << "===========" << std::endl;
std::cout << *newLocal->fingerprint << std::endl;
*/
    // Estimate the Lipinski parameters.
    newLocal->estimateLipinski(*this, that);

// exit(0);

    return newLocal;
}

// *****************************************************************************

//
// On-Demand acquisition of the fingerprint one time.
//
FragmentGraph* Molecule::getFingerprint() const
{
    return this->fingerprint;
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

/*
    atoms[xIndex].addConnection(yIndex);
    atoms[yIndex].addConnection(xIndex);
*/

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
}

// *****************************************************************************

std::string Molecule::toString() const
{
    std::ostringstream oss;

    oss << "Molecule: " << uniqueIndexID << " ";

    if (IsLinker()) oss << " is a linker.";
    else if (IsRigid())  oss << " is a rigid.";
    else if (IsComplex())
    {
        oss << "There are ";
        oss << NumRigids() << " rigids, and ";
        oss << NumLinkers();
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

std::ostream& operator<< (std::ostream& os, const Molecule& mol)
{
    os << mol.toString() << std::endl;

    return os;
}
