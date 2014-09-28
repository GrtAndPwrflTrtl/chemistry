#ifndef _MOLECULE_GUARD
#define _MOLECULE_GUARD 1


#include <string>
#include <vector>
#include <memory>


#include <openbabel/mol.h>


#include "Bond.h"
#include "Atom.h"
#include "IdFactory.h"
#include "obgen.h"
// #include "Thread_Pool.h"
#include "Constants.h"


class EdgeAggregator;
class Rigid;
class Linker;


class Molecule
{
  public:
    Molecule();
    Molecule(OpenBabel::OBMol* mol, const std::string& name, MoleculeT t);

    ~Molecule();

    void setGraphID(unsigned int id) { graphID = id; }
    unsigned int getGraphID() const { return graphID; }

    //
    // Get Functions
    //
    unsigned int getMoleculeID() const { return this->moleculeID; }
    bool IsLinker() const { return type == LINKER; }
    bool IsComplex() const { return type == COMPLEX; }
    bool IsRigid() const { return type == RIGID; }
    
    bool islipinskiPredicted() { return lipinskiPredicted; };
    bool islipinskiEstimated() { lipinskiEstimated; };
    bool isLipinskiCompliant();
    bool isOpenbabelLipinskiCompliant();
    double getMolWt() const { return MolWt;};
    double getHBD() const { return HBD;};
    double getHBA1() const { return HBA1;};
    double getlogP() const { return logP;};

    int getNumberOfAtoms() const { return this->atoms.size(); }
    int getNumberOfBonds() const { return this->bonds.size(); }
    int getNumberOfRigids() const { return this->rigids.size(); }
    int getNumberOfLinkers() const { return this->linkers.size(); }

    Atom getAtom(int id) const;
    Bond getBond(int id) const;
    Bond getBond(int xID, int yID) const;

    OpenBabel::OBMol* getOpenBabelMol() const { return obmol; }
    void GetFingerprint(std::vector<unsigned int>& fp) const;


    //
    // Set Functions
    //
    void setMoleculeID(int x) { this->moleculeID = x; }
    void makeLinker() { type = LINKER; }
    void makeRigid() { type = RIGID; }
    void makeComplex() { type = COMPLEX; }
    
    bool addBond(int xID, int yID); //, eTypeOfBondT bt, eStatusBitT s);
    void addAtom(const Atom& a);

    std::string toString() const;
    friend std::ostream& operator<< (std::ostream& os, Molecule& mol);
    virtual bool operator==(const Molecule& that) const;
    std::vector<EdgeAggregator*>* Compose(const Molecule&) const;

    const std::vector<Rigid*>& getRigids() const { return rigids; }
    const std::vector<Linker*>& getLinkers() const { return linkers; }
    std::string getName() const { return name; }

    void predictLipinski();
    void estimateLipinski(const Molecule &mol1, const Molecule &mol2);

  private:
    void localizeOBMol();

    static const unsigned int RIGID_UPPER_BOUND = 10;
    static const unsigned int LINKER_UPPER_BOUND = 6;
    bool exceedsMaxMolecularMass();
    bool ContainsLoops() const;
    bool satisfiesMoleculeSynthesisCriteria();
    Molecule* ComposeToNewMolecule(const Molecule& that,
                                   int thisAtomIndex,
                                   int thatAtomIndex) const;

//    static Thread_Pool<OpenBabel::OBMol*, bool> pool; // (THREAD_POOL_SIZE, OBGen::obgen);

  protected:
    unsigned int graphID;
    unsigned int moleculeID;
    std::string name;
    MoleculeT type;
    

    IdFactory atomIdMaker;
    IdFactory bondIdMaker;

    OpenBabel::OBMol* obmol;
    // Used for molecular comparison of this OBMol object.
    std::vector<unsigned int> fingerprint; 


    std::vector<Atom> atoms;
    std::vector<Bond> bonds;
   
    std::vector<Rigid*> rigids;
    std::vector<Linker*> linkers;

    int getAtomIndex(int id) const;
    int getBondIndex(int id) const;
    int getBondIndex(int xID, int yID) const;

    // Lipinski Descriptors
    bool lipinskiPredicted, lipinskiEstimated; // whether Lipinski values are calculated with Predict(), or simple estimation
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
