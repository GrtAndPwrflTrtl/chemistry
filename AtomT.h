#ifndef _ATOM_TYPE_GUARD
#define _ATOM_TYPE_GUARD 1

#include <string>
#include <cctype>
#include <cstdlib>

typedef enum
{
  CARBON,
  CHLORINE,
  HYDROGEN,
  NITROGEN,
  OXYGEN,
  PHOSPHORUS,
  SULFUR,
  UNKNOWN
} AtomEnumT;

typedef enum
{
  AM,
  AROMATIC,
  NONE,
  PL,
} SpecialEnumT;

//
// Light-weight aggregator for Atoms and their types
//
class AtomT
{
  public:
    AtomEnumT atomType;
    int specificNum;
    SpecialEnumT specialT;

    AtomT(AtomEnumT type = UNKNOWN, int val = -1, SpecialEnumT spec = NONE);
    AtomT(const std::string&);
 
    std::string toString() const;
    friend std::ostream& operator<< (std::ostream& os, const AtomT& atomType);
    bool operator==(const AtomT& that) const;
    bool operator!=(const AtomT& that) const { return !(*this == that); }

  private:
    static AtomEnumT convertToAtomEnum(std::string&);
    static SpecialEnumT convertToSpecialEnum(std::string& s);
};

#endif
