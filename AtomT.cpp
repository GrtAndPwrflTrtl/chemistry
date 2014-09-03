#include <string>
#include <sstream>
#include <cctype>
#include <cstdlib>
#include <iostream>

#include "AtomT.h"

AtomT::AtomT(AtomEnumT type, int val, SpecialEnumT spec) : atomType(type),
                                                           specificNum(val),
                                                           specialT(spec)
{
}

AtomEnumT AtomT::convertToAtomEnum(std::string& s)
{
    if (s == "C") return CARBON;
    else if (s == "Cl") return CHLORINE;
    else if (s == "H") return HYDROGEN;
    else if (s == "N") return NITROGEN;
    else if (s == "O") return OXYGEN;
    else if (s == "P") return PHOSPHORUS;
    else if (s == "S") return SULFUR;

    return UNKNOWN;
}

SpecialEnumT AtomT::convertToSpecialEnum(std::string& s)
{
    if (s == "am") return AM;
    else if (s == "ar") return AROMATIC;
    else if (s == "pl") return PL;

    return NONE;
}

//
// Split the string into the constituent elements: <ELEMENT>.(<Specification> | <VAL>)
//
void parseString(const std::string& s, std::string& prefix, std::string& suffix, int& val)
{
    suffix = "";
    val = -1;
    
    // Look for the period
    int pos = s.find(".");

    // If no period, then we have just an element
    if (pos == std::string::npos)
    {
        prefix = s;
        return;
    }

    prefix = s.substr(0, pos);

    suffix = s.substr(pos+1);

// std::cerr << "Type Suffix: " << suffix << std::endl;

    // Check if the last position is a numeric value
    std::string lastChar = s.substr(s.size() - 1);

    if (isdigit(lastChar[0]))
    {
        val = atoi(lastChar.c_str());
        // Cut down the suffix string by one character.
        suffix = suffix.substr(0, suffix.size() - 1);
    }
}

//
// An AtomT is of the form: <Atom-Type>.<int-type> OR
//                          <Atom-Type>.<string-type>
AtomT::AtomT(const std::string& in)
{
    std::string prefix = "";
    std::string suffix = "";

    parseString(in, prefix, suffix, this->specificNum);

    this->atomType = convertToAtomEnum(prefix);
    this->specialT = convertToSpecialEnum(suffix);
}

// **************************************************************************************

bool AtomT::operator==(const AtomT& that) const
{
//std::cerr << *this << " == ";
//std::cerr << that << ": ";

    if (this->atomType != that.atomType) return false;

    if (this->specificNum != that.specificNum) return false;

    if (this->specialT != that.specialT) return false;

//std::cerr << "true" << std::endl;

    return true;
}


// **************************************************************************************

std::string AtomT::toString() const
{
    std::ostringstream oss;

    switch(atomType)
    {
      case CARBON:      oss << "C"; break;
      case CHLORINE:    return "Cl";
      case HYDROGEN:    oss << "H"; break;
      case NITROGEN:    oss << "N"; break;
      case OXYGEN:      oss << "O"; break;
      case PHOSPHORUS:  oss << "P"; break;
      case SULFUR:      oss << "S"; break;
    }

    oss << ".";

    switch(specialT)
    {
      case AM:          oss << "am"; break;
      case AROMATIC:    oss << "ar"; break;
      case PL:          oss << "pl"; break;
    }

    if (specificNum != -1) oss << specificNum;

    return oss.str();
}


std::ostream& operator<< (std::ostream& os, const AtomT& atomType)
{
    os << atomType.toString();

    return os;
}
