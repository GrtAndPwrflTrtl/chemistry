#ifndef _BOND_TYPE_H_
#define _BOND_TYPE_H_ 1
#include "bondType.h"
#endif

#ifndef _BOND_TYPE_
#define _BOND_TYPE_ 1

void Bond_set_values(bondType *b, int id_p, int originAtom_p, int targetAtom_p, char bondType_p[],  char statusBit_p[], enum eFileSources fSource)
{
	/* set literal values */
	b->id = id_p;
	b->originAtomID = originAtom_p;
	b->targetAtomID = targetAtom_p;
	b->fSource = fSource;
	
	/* Because the bondType is a string with only 8 possible values, I'll store it as an integer */
	if(strcmp(bondType_p, "1") == 0)
		b->typeOfBond = 1;
	if(strcmp(bondType_p, "2") == 0)
		b->typeOfBond = 2;
	if(strcmp(bondType_p, "3") == 0)
		b->typeOfBond = 3;
	if(strcmp(bondType_p, "am") == 0)
		b->typeOfBond = 4;
	if(strcmp(bondType_p, "ar") == 0)
		b->typeOfBond = 5;
	if(strcmp(bondType_p, "du") == 0)
		b->typeOfBond = 6;
	if(strcmp(bondType_p, "un") == 0)
		b->typeOfBond = 7;
	if(strcmp(bondType_p, "nc") == 0)
		b->typeOfBond = 8;
	
	/* Because the statusBit is a string with only 6 possible values, I'll store it as an integer */
	b->statusBit = -1;	/* the default value, if the status bit is unassigned */
	if(strcmp(statusBit_p, "TYPECOL") == 0)
		b->statusBit = 1;
	if(strcmp(statusBit_p, "GROUP") == 0)
		b->statusBit = 2;
	if(strcmp(statusBit_p, "CAP") == 0)
		b->statusBit = 3;
	if(strcmp(statusBit_p, "BACKBONE") == 0)
		b->statusBit = 4;
	if(strcmp(statusBit_p, "DICT") == 0)
		b->statusBit = 5;
	if(strcmp(statusBit_p, "INTERRES") == 0)
		b->statusBit = 6;

} /* Bond_set_values */


/****************************************************************************************/

int parseBond(bondType *b, char buffer[], enum eFileSources fSource)
{
	/* 
		Input: line from file
		Output: bond
	*/
	
	int id;
	int originAtom;
	int targetAtom;
	char bondType[500];
	char statusBit[500];
	int index = 0;

	/* initialize */
	strcpy(statusBit, "");

	/* get id */
	while(buffer[index] == ' ' || buffer[index] == '\t')		//get past white space
		index++;

	id = buffer[index] - '0';
	index++;

	while(buffer[index] != ' ' && buffer[index] != '\t' && buffer[index] != '\n' && buffer[index] != '\0')
	{ /* parse number */

		id *= 10;
		id += buffer[index] - '0';
		index++;
	} /* parse number */

	while(buffer[index] == ' ' || buffer[index] == '\t')		//get past white space
		index++;

	/* get origin atom */
	originAtom = buffer[index] - '0';
	index++;

	while(buffer[index] != ' ' && buffer[index] != '\t' && buffer[index] != '\n' && buffer[index] != '\0')
	{ /* parse number */

		originAtom *= 10;
		originAtom += buffer[index] - '0';
		index++;

	} /* parse number */

	while(buffer[index] == ' ' || buffer[index] == '\t')		//get past white space
		index++;

	/* get target atom */
	targetAtom = buffer[index] - '0';
	index++;

	while(buffer[index] != ' ' && buffer[index] != '\t' && buffer[index] != '\n' && buffer[index] != '\0')
	{ /* parse number */

		targetAtom *= 10;
		targetAtom += buffer[index] - '0';
		index++;

	} /* parse number */

	while(buffer[index] == ' ' || buffer[index] == '\t')		//get past white space
		index++;

	/* get bond type */
	sscanf(buffer + index, " %s", bondType);
	index += strlen(bondType);

	if(buffer[index] != '\n' && buffer[index] != '\0')
	{ /* if not end of line */

		/* get status bit */
		sscanf(buffer + index, " %s", statusBit);

	} /* if not end of line */
	
	/* set bond values */
	Bond_set_values(b, id, originAtom, targetAtom, bondType, statusBit, fSource);

} /* parseBond */

/****************************************************************************************/

#endif