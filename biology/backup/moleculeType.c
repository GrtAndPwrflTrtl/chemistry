#ifndef _MOLECULE_TYPE_H_
#define _MOLECULE_TYPE_H_ 1
#include "moleculeType.h"
#endif

#ifndef _MOLECULE_TYPE_
#define _MOLECULE_TYPE_ 1

void Molecule_set_values(moleculeType *m, char name_p[], int numAtoms_p, int numBonds_p, int numSubStructs_p, int numFeatures_p, int numSets_p, char moleculeType_p[], char chargeType_p[], char statusBit_p[], char comments_p[], enum eFileSources fSource)
{
	/* set literal values */
	strcpy(m->name, name_p);
	m->numAtoms = numAtoms_p;
	m->numBonds = numBonds_p;
	m->numSubStructs = numSubStructs_p;
	m->numFeatures = numFeatures_p;
	m->numSets = numSets_p;
	strcpy(m->comments, comments_p);
	m->fSource = fSource;

	/* Because the moleculeType is a string with only 5 possible values, I'll store it as an integer */
	m->typeOfMolecule = -1; /* default value, if the moleculeType is unassigned */
	if(strcmp(moleculeType_p, "SMALL") == 0)
		m->typeOfMolecule = 1;
	if(strcmp(moleculeType_p, "BIOPOLYMER") == 0)
		m->typeOfMolecule = 2;
	if(strcmp(moleculeType_p, "PROTEIN") == 0)
		m->typeOfMolecule = 3;
	if(strcmp(moleculeType_p, "NUCLEIC_ACID") == 0)
		m->typeOfMolecule = 4;
	if(strcmp(moleculeType_p, "SACCHARIDE") == 0)
		m->typeOfMolecule = 5;

	/* Because the chargeType is a string with only 12 possible values, I'll store it as an integer */
	m->chargeType = -1; /* default value, if the chargeType is unassigned */
	if(strcmp(chargeType_p, "NO_CHARGES") == 0)
		m->chargeType = 1;
	if(strcmp(chargeType_p, "DEL_RE") == 0)
		m->chargeType = 2;
	if(strcmp(chargeType_p, "GASTEIGER") == 0)
		m->chargeType = 3;
	if(strcmp(chargeType_p, "GAST_HUCK") == 0)
		m->chargeType = 4;
	if(strcmp(chargeType_p, "HUCKEL") == 0)
		m->chargeType = 5;
	if(strcmp(chargeType_p, "PULLMAN") == 0)
		m->chargeType = 6;
	if(strcmp(chargeType_p, "GAUSS80_CHARGES") == 0)
		m->chargeType = 7;
	if(strcmp(chargeType_p, "AMPAC_CHARGES") == 0)
		m->chargeType = 8;
	if(strcmp(chargeType_p, "MULLIKEN_CHARGES") == 0)
		m->chargeType = 9;
	if(strcmp(chargeType_p, "DICT_CHARGES") == 0)
		m->chargeType = 10;
	if(strcmp(chargeType_p, "MMFF94_CHARGES") == 0)
		m->chargeType = 11;
	if(strcmp(chargeType_p, "USER_CHARGES") == 0)
		m->chargeType = 12;

	/* Because the statusBits is a string with only 6 possible values, I'll store it as an integer */
	m->statusBits = -1; /* default value, if the statusBits is unassigned */
	if(strcmp(statusBit_p, "system") == 0 || strcmp(statusBit_p, "System") == 0 || strcmp(statusBit_p, "SYSTEM") == 0)
		m->statusBits = 1;
	if(strcmp(statusBit_p, "invalid_charges") == 0 || strcmp(statusBit_p, "Invalid_charges") == 0 || strcmp(statusBit_p, "INVALID_CHARGES") == 0)
		m->statusBits = 2;
	if(strcmp(statusBit_p, "analyzed") == 0 || strcmp(statusBit_p, "Analyzed") == 0 || strcmp(statusBit_p, "ANALYZED") == 0)
		m->statusBits = 3;
	if(strcmp(statusBit_p, "substituted") == 0 || strcmp(statusBit_p, "Substituted") == 0 || strcmp(statusBit_p, "SUBSTITUTED") == 0)
		m->statusBits = 4;
	if(strcmp(statusBit_p, "altered") == 0 || strcmp(statusBit_p, "Altered") == 0 || strcmp(statusBit_p, "ALTERED") == 0)
		m->statusBits = 5;
	if(strcmp(statusBit_p, "ref_angle") == 0 || strcmp(statusBit_p, "Ref_angle") == 0 || strcmp(statusBit_p, "REF_ANGLE") == 0)
		m->statusBits = 6;

} /* Molecule_set_values */

/****************************************************************************************/

int parseMolecule(moleculeType *m, char moleculeName[], char moleculeStats[], char typeOfMolecule[], char moleculeCharge[], char moleculeStatusBits[], char comments[], enum eFileSources fSource)
{
	/*
		Input: line from file, and open file that is currently in a molecule section
		Output: molecule, and the file in the input has advanced to the next section
	*/

	int numAtoms_p; 
	int numBonds_p; 
	int numSubStructs_p; 
	int numFeatures_p; 
	int numSets_p;
	int index = 0;

	/* get number of atoms from stats */
	while(moleculeStats[index] == ' ' || moleculeStats[index] == '\t')	//get past white space
		index++;
	numAtoms_p = moleculeStats[index] - '0';
	index++;
	while(moleculeStats[index] != ' ' && moleculeStats[index] != '\t' && moleculeStats[index] != '\n' && moleculeStats[index] != '\0')
	{ /* parse number */
		numAtoms_p *= 10;
		numAtoms_p += moleculeStats[index] - '0';
		index++;
	} /* parse number */

	while(moleculeStats[index] == ' ' || moleculeStats[index] == '\t')	//get past white space but not newline
		index++;
	
	if(moleculeStats[index] != '\n' && moleculeStats[index] != '\0')
	{ /* if not end of line */

		/* get number of bonds */
		numBonds_p = moleculeStats[index] - '0';
		index++;
		while(moleculeStats[index] != ' ' && moleculeStats[index] != '\t' && moleculeStats[index] != '\n' && moleculeStats[index] != '\0')
		{ /* parse number */
			numBonds_p *= 10;
			numBonds_p += moleculeStats[index] - '0';
			index++;
		} /* parse number */

	} /* if not end of line */

	while(moleculeStats[index] == ' ' || moleculeStats[index] == '\t')	//get past white space but not newline
		index++;

	if(moleculeStats[index] != '\n' && moleculeStats[index] != '\0')
	{ /* if not end of line */

		numSubStructs_p = moleculeStats[index] - '0';
		index++;
		while(moleculeStats[index] != ' ' && moleculeStats[index] != '\t' && moleculeStats[index] != '\n' && moleculeStats[index] != '\0')
		{ /* get num_subst */
			numSubStructs_p *= 10;
			numSubStructs_p += moleculeStats[index] - '0';
			index++;	
		} /* parse number */
	} /* if not end of line */

	while(moleculeStats[index] == ' ' || moleculeStats[index] == '\t')	//get past white space but not newline
		index++;

	if(moleculeStats[index] != '\n' && moleculeStats[index] != '\0')
	{ /* if not end of line */
		
		numFeatures_p = moleculeStats[index] - '0';
		index++;
		while(moleculeStats[index] != ' ' && moleculeStats[index] != '\t' && moleculeStats[index] != '\n' && moleculeStats[index] != '\0')
		{ /* get num_feat */
			numFeatures_p *= 10;
			numFeatures_p += moleculeStats[index] - '0';
			index++;
		} /* parse number */
	} /* if not end of line */

	while(moleculeStats[index] == ' ' || moleculeStats[index] == '\t')	//get past white space but not newline
		index++;

	if(moleculeStats[index] != '\n' && moleculeStats[index] != '\0')
	{ /* if not end of line */
		
		numSets_p = moleculeStats[index] - '0';
		index++;
		while(moleculeStats[index] != ' ' && moleculeStats[index] != '\t' && moleculeStats[index] != '\n' && moleculeStats[index] != '\0')
		{ /* get num_set */
			numSets_p *= 10;
			numSets_p += moleculeStats[index] - '0';
			index++;
		} /* parse number */
	} /* if not end of line */

	/* set molecule values */
	Molecule_set_values(m, moleculeName, numAtoms_p, numBonds_p, numSubStructs_p, numFeatures_p, numSets_p, typeOfMolecule, moleculeCharge, moleculeStatusBits, comments, fSource);

} /* parseMolecule */

/****************************************************************************************/

#endif