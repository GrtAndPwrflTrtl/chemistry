#ifndef _MOLECULE_GUARD
#define _MOLECULE_GUARD 1


#include <string>
#include <cstring> // for memset
#include <vector>
#include <memory>
#include <map>
#include <pthread.h>


#include <openbabel/mol.h>


#include "Bond.h"
#include "Atom.h"
#include "IdFactory.h"
#include "obgen.h"
#include "Constants.h"


class EdgeAggregator;
class Rigid;
class Linker;
class FragmentGraph;

class Molecule
{
  public:
    Molecule();
    Molecule(OpenBabel::OBMol* mol, const std::string& name, MoleculeT t);

    ~Molecule();

    void setUniqueIndexID(unsigned int id) { uniqueIndexID = id; }
    unsigned int getUniqueIndexID() const { return uniqueIndexID; }

    //
    // Get Functions
    //
    //unsigned int getMoleculeID() const { return this->moleculeID; }
    bool IsLinker() const { return type == LINKER; }
    bool IsComplex() const { return type == COMPLEX; }
    bool IsRigid() const { return type == RIGID; }

    bool islipinskiPredicted() { return lipinskiPredicted; };
    bool islipinskiEstimated() { lipinskiEstimated; };
    bool isLipinskiCompliant();
    bool isOpenBabelLipinskiCompliant();
    double getMolWt() const { return MolWt; }
    double getHBD() const { return HBD; }
    double getHBA1() const { return HBA1; }
    double getlogP() const { return logP; }

    int getNumberOfAtoms() const { return this->atoms.size(); }
    int getNumberOfBonds() const { return this->bonds.size(); }

    Atom getAtom(int id) const;
    Bond getBond(int id) const;
    Bond getBond(int xID, int yID) const;

    OpenBabel::OBMol* getOpenBabelMol() const { return obmol; }
    FragmentGraph* getFingerprint() const;

    bool addBond(int xID, int yID); //, eTypeOfBondT bt, eStatusBitT s);
    void addAtom(const Atom& a);

    std::string toString() const;
    friend std::ostream& operator<< (std::ostream& os, const Molecule& mol);
    virtual bool operator==(const Molecule& that) const;
    std::vector<EdgeAggregator*>* Compose(const Molecule&) const;

    void getRigids(std::vector<Rigid*>& rigids) const { rigids = this->rigids; }
    void getLinkers(std::vector<Linker*>& linkers) const { linkers = this->linkers; } 

    unsigned int NumRigids() const { return numRigids; }
    unsigned int NumLinkers() const { return numLinkers; }
    unsigned int NumUniqueRigids() const { return numUniqueRigids; }
    unsigned int NumUniqueLinkers() const { return numUniqueLinkers; }

    std::string getName() const { return name; }

    void openBabelPredictLipinski();
    static bool isOpenBabelLipinskiCompliant(OpenBabel::OBMol& mol);
    void estimateLipinski(const Molecule &mol1, const Molecule &mol2);
    static bool willExceedMolecularWeight(const Molecule &mol1, const Molecule &mol2);

    // The 'size' of a molecule is based on the number of total fragments.
    unsigned int size() const { return numLinkers + numRigids; }

    // Initialize any containers to track fragments (linkers / rigids)
    void initFragmentDevices();

    // Calculate the number of linkers / rigids (copies and unique)
    void calcFragmentInfo();

    // Initialize the fragment container
    void initFragmentInfo();

    // Initialize the graph-based representation of the fragment
    void initGraphRepresentation();

    // Collection of linkers and rigids for this synthesis.
    static std::vector<Molecule*> baseMolecules;
    static void SetBaseMoleculeInfo(const std::vector<Molecule*> baseMols,
                                    unsigned int numRigids, unsigned int numLinkers); 

    virtual unsigned int getFragmentId() const { return -1; }

    void getConnectionIDs(std::vector<unsigned int>& conns) const { conns = connectionIDs; }

    static unsigned int NUM_UNIQUE_FRAGMENTS;

    // Lock openbabel
    void init_openbabel_lock();
    static pthread_mutex_t openbabel_lock;

  private:
    void localizeOBMol();

    bool exceedsMaxEstimatedThresholds();
    bool ContainsLoops() const;
    bool satisfiesMoleculeSynthesisCriteria();
    Molecule* ComposeToNewMolecule(const Molecule& that,
                                   int thisAtomIndex,
                                   int thatAtomIndex) const;

    static unsigned int RIGID_INDEX_START;
    static unsigned int RIGID_INDEX_END;
    static unsigned int LINKER_INDEX_START;
    static unsigned int LINKER_INDEX_END;
    static unsigned int FRAGMENT_END_INDEX;

  protected:
    unsigned int uniqueIndexID;
    std::string name;
    MoleculeT type;
    
    IdFactory atomIdMaker;
    IdFactory bondIdMaker;

    // Unique index for an instance of a linker / rigid; not applicable to a complex molecule. 
    unsigned int localMoleculeID;

    // Each linker / rigid has connection points; we create unique ids for those connections. 
    static IdFactory connectionIdMaker;

    // Array of connection ids for the linker / rigid; parallels the atom vector
    std::vector<unsigned int> connectionIDs;

    // An array used to count the number of each specific linker / rigid in this molecule
    unsigned int* fragmentCounter;
    // Descriptors of this molecule (used as heuristics in the fingerprint)
    unsigned int numLinkers;
    unsigned int numUniqueLinkers;
    unsigned int numRigids;
    unsigned int numUniqueRigids;

    // The rigids and linkers in this molecule
    std::vector<Rigid*> rigids;
    std::vector<Linker*> linkers;

    // Open Babel representation of this molecule.
    OpenBabel::OBMol* obmol;

    // Used for molecular comparison
    FragmentGraph* fingerprint;

    // Local atoms and bonds
    std::vector<Atom> atoms;
    std::vector<Bond> bonds;
   
    int getAtomIndex(int id) const;
    int getBondIndex(int id) const;
    int getBondIndex(int xID, int yID) const;

    // Lipinski Descriptors
    // whether Lipinski values are calculated with Predict(), or simple estimation
    bool lipinskiPredicted, lipinskiEstimated;
    double MolWt;
    double HBD;
    double HBA1;
    double logP;

    virtual void parseAppendix(std::string& comment)
    {
        std::cerr << "Called Wrong parseAppendix::MOLECULE" << std::endl;
    }
};

#endif
