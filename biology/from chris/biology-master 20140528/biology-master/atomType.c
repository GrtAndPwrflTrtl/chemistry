#ifndef _ATOM_TYPE_H_
#define _ATOM_TYPE_H_ 1
#include "atomType.h"
#endif

#ifndef _ATOM_TYPE_
#define _ATOM_TYPE_ 1

void Atom_set_values(atomType *a, int id_p, char name_p[], float x_p, float y_p, float z_p, char atomType_p[], int subst_id_p, char subst_name_p[], float charge_p, char statusBit_p[], enum eFileSources fSource)
{
	/* set literal values */
	a->id = id_p;
	strcpy(a->name, name_p);
	a->x = x_p;
	a->y = y_p;
	a->z = z_p;
	strcpy(a->typeOfAtom, atomType_p);
	a->subst_id = subst_id_p;
	strcpy(a->subst_name, subst_name_p);
	a->charge = charge_p;
	a->fSource = fSource;
	
	/* Because the statusBit is a string with only 8 possible values, I'll store it as an integer */
	a->statusBit = -1;	/* the default value, if the status bit is unassigned */
	if(strcmp(statusBit_p, "DSPMOD") == 0)
		a->statusBit = 1;
	if(strcmp(statusBit_p, "TYPECOL") == 0)
		a->statusBit = 2;
	if(strcmp(statusBit_p, "CAP") == 0)
		a->statusBit = 3;
	if(strcmp(statusBit_p, "BACKBONE") == 0)
		a->statusBit = 4;
	if(strcmp(statusBit_p, "DICT") == 0)
		a->statusBit = 5;
	if(strcmp(statusBit_p, "ESSENTIAL") == 0)
		a->statusBit = 6;
	if(strcmp(statusBit_p, "WATER") == 0)
		a->statusBit = 7;
	if(strcmp(statusBit_p, "DIRECT") == 0)
		a->statusBit = 8;
} /* Atom_set_values */

/****************************************************************************************/

int parseAtom(atomType *a, char buffer[], enum eFileSources fSource)
{
	/*
		Input: line from file 
		Output: atom
	*/

	int index = 0;
	int endOfString;
	char temp[500];
	int tempLength;
	int id, subst_id;
	char name[100], typeOfAtom[10], subst_name[100], statusBit[100];
	float x, y, z, charge;

	/* initialization*/
	subst_id = -1;
	strcpy(subst_name, "");
	charge = 0.0;
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

	/* get name */
	sscanf(buffer + index, " %s", name);
	index += strlen(name);

	while(buffer[index] == ' ' || buffer[index] == '\t')		//get past white space
		index++;

	/* get x */
	sscanf(buffer + index, " %f", &x);
	while(buffer[index] != ' ' && buffer[index] != '\t')		//get past non-white space
		index++;

	while(buffer[index] == ' ' || buffer[index] == '\t')		//get past white space
		index++;

	/* get y */
	sscanf(buffer + index, " %f", &y);
	while(buffer[index] != ' ' && buffer[index] != '\t')		//get past non-white space
		index++;

	while(buffer[index] == ' ' || buffer[index] == '\t')		//get past white space
		index++;

	/* get z */
	sscanf(buffer + index, " %f", &z);
	while(buffer[index] != ' ' && buffer[index] != '\t')		//get past non-white space
		index++;

	while(buffer[index] == ' ' || buffer[index] == '\t')		//get past white space
		index++;

	/* get type of item */
	sscanf(buffer + index, " %s", typeOfAtom);
	index += strlen(typeOfAtom);

	while(buffer[index] == ' ' || buffer[index] == '\t')		//get past white space
		index++;

	if(buffer[index] != '\n' && buffer[index] != '\0')
	{ /* if not end of line */

		/* get subst_id */
		subst_id = buffer[index] - '0';
		index++;

		while(buffer[index] != ' ' && buffer[index] != '\t' && buffer[index] != '\n' && buffer[index] != '\0')
		{ /* parse number */

			subst_id *= 10;
			subst_id += buffer[index] - '0';
			index++;

		} /* parse number */
	} /* if not end of line */

	while(buffer[index] == ' ' || buffer[index] == '\t')		//get past white space
		index++;

	if(buffer[index] != '\n' && buffer[index] != '\0')
	{ /* if not end of line */

		/* get subst_name */
		sscanf(buffer + index, " %s", subst_name);
		index += strlen(subst_name);

	} /* if not end of line */

	while(buffer[index] == ' ' || buffer[index] == '\t')		//get past white space
		index++;

	if(buffer[index] != '\n' && buffer[index] != '\0')
	{ /* if not end of line */

		/* get charge */
		sscanf(buffer + index, " %f", &charge);
		while(buffer[index] != ' ' && buffer[index] != '\t' && buffer[index] != '\n' && buffer[index] != '\0')		//get past non-white space
			index++;

	} /* if not end of line */

	while(buffer[index] == ' ' || buffer[index] == '\t')		//get past white space
		index++;

	if(buffer[index] != '\n' && buffer[index] != '\0')
	{ /* if not end of line */

		/* get status_bit */
		sscanf(buffer + index, " %s", statusBit);

	} /* if not end of line */

	/* set atom values */
	Atom_set_values(a, id, name, x, y, z, typeOfAtom, subst_id, subst_name, charge, statusBit, fSource);

} /* parseAtom */

/****************************************************************************************/

#endif