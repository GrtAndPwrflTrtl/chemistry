/* a small change in the biology folder */
#ifndef _STDIO_H_
#define _STDIO_H_ 1
#include <stdio.h>
#endif

#ifndef _STRING_H_
#define _STRING_H_ 1
#include <string.h>
#endif

#ifndef _MOLECULE_TYPE_H_
#define _MOLECULE_TYPE_H_ 1
#include "moleculeType.h"
#endif

#ifndef _BOND_TYPE_H_
#define _BOND_TYPE_H_ 1
#include "bondType.h"
#endif

#ifndef _ATOM_TYPE_H_
#define _ATOM_TYPE_H_ 1
#include "atomType.h"
#endif

#ifndef _AVL_TYPE_H_
#define _AVL_TYPE_H_1
#include "avl.h"
#endif

/****************************************************************************************/

#ifndef _ENUM_FILE_SOURCE
#define _ENUM_FILE_SOURCE 1
enum eFileSources
{
	active=1, decoy=2, crystal=3
};
#endif

/****************************************************************************************/

int compareAtoms(atomType *a, atomType *b);
int compareBonds(bondType *a, bondType *b);
int compareMolecules(moleculeType *a, moleculeType *b);
int getAtoms(FILE *f, int *numberOfAtoms, int *numberOfBonds, int *numberOfMolecules, enum eFileSources fSource);
int getBonds(FILE *f, int *numberOfAtoms, int *numberOfBonds, int *numberOfMolecules, enum eFileSources fSource);
int getMolecules(FILE *f, int *numberOfAtoms, int *numberOfBonds, int *numberOfMolecules, enum eFileSources fSource);
void readFile(FILE *f, int *numberOfAtoms, int *numberOfBonds, int *numberOfMolecules, enum eFileSources fSource);
int showFile(FILE *f, int x);
void visitAtom(atomType *a);

/****************************************************************************************/

int main(void)
{
	FILE *fileNames;
	FILE *currentFile;
	FILE *f;
	//AVL_TREE *decoyAtoms;
	AVL_TREE *crystalAtoms;
	//AVL_TREE *activeAtoms;
	AVL_TREE *crystalBonds;
	AVL_TREE *crystalMolecules;

	atomType *currentAtom;
	bondType *currentBond;
	moleculeType *currentMolecule;
	enum eFileSources fSource;
	char buffer[500], nameOfFile[20];
	char moleculeName[500], moleculeStats[500], typeOfMolecule[500], moleculeCharge[500], moleculeStatusBits[500], comments[500];
	int sectionFlag, numberOfAtoms = 0, numberOfBonds = 0, numberOfMolecules = 0;

	if(1)
	{ /* comment this out */

		/* open crystal file */
		//f = fopen("crystal_ligand.mol2", "r");
		//f = fopen("actives_final.mol2", "r");
		f = fopen("decoys_final.mol2", "r");

		/* create crystals */
		crystalAtoms = AVL_Create(compareAtoms);
		crystalBonds = AVL_Create(compareBonds);
		crystalMolecules = AVL_Create(compareMolecules);

		/* read first line in file */
		fscanf(f, " %[^\n]", buffer);
		if(buffer[9] == 'A')
			sectionFlag = 1;
		if(buffer[9] == 'B')
			sectionFlag = 2;
		if(buffer[9] == 'M')
			sectionFlag = 3;

		while(!feof(f) && sectionFlag != -1)
		{ /* while not eof */
		
			/* if file is in bonds section, then get bonds */
			if(sectionFlag == 2)
			{ /* bonds */
				//sectionFlag = getBonds(f, &numberOfAtoms, &numberOfBonds, &numberOfMolecules, crystal);

				/* read line from file */
				fscanf(f, " %[^\n]", buffer);

				while(!feof(f) && buffer[0] != '@')
				{ /* while not eof and line doesn't start with @ */

					/* parse bond */
					currentBond = (bondType *)malloc(sizeof(bondType));
					parseBond(currentBond, buffer, crystal);

					currentBond->id = 1 + AVL_Count(crystalBonds);
					/* add bond to crystal */
					if(!AVL_Insert(crystalBonds, currentBond))
					{
						printf("ran out of memory\n");
						exit(0);
					}

					/* read line from file */
					fscanf(f, " %[^\n]", buffer);
					numberOfBonds = AVL_Count(crystalBonds);
					if(currentBond->id % 1000 == 0)
						printf("currentBond = %d\tnumber of bonds = %d\n", currentBond->id, numberOfBonds);

				} /* while not eof and line doesn't start with @ */

				if(feof(f))
					sectionFlag = -1;
				if(buffer[9] == 'A')
					sectionFlag = 1;
				if(buffer[9] == 'B')
					sectionFlag = 2;
				if(buffer[9] == 'M')
					sectionFlag = 3;

			} /* bonds */

			/* if file is in molecules section, then get molecules */
			else if(sectionFlag == 3)
			{ /* molecules */

				//sectionFlag = getMolecules(f, &numberOfAtoms, &numberOfBonds, &numberOfMolecules, crystal);

				/* initialize strings */
				strcpy(moleculeName, "");
				strcpy(moleculeStats, "");
				strcpy(typeOfMolecule, "");
				strcpy(moleculeCharge, "");
				strcpy(moleculeStatusBits, "");
				strcpy(comments, "");

				/* read name of molecule */
				fscanf(f, " %[^\n]", moleculeName);

				/* read molecule stats line */
				fscanf(f, " %[^\n]", moleculeStats);

				/* read molecule type */
				fscanf(f, " %[^\n]", typeOfMolecule);

				/* read charge type */
				fscanf(f, " %[^\n]", moleculeCharge);

				/* read line into buffer */
				fscanf(f, " %[^\n]", buffer);

				if(buffer[0] != '@')
				{ /* if buffer doesn't start with @ */

					/* copy buffer into status bits */
					strcpy(moleculeStatusBits, buffer);

					/* read line into buffer */
					fscanf(f, " %[^\n]", buffer);

					if(buffer[0] != '@')
					{ /* if buffer doesn't start with @ */

						/* copy buffer into comments */
						strcpy(comments, buffer);

					} /* if buffer doesn't start with @ */

				} /* if buffer doesn't start with @ */

				/* parse molecule */
				currentMolecule = (moleculeType *)malloc(sizeof(moleculeType));
				parseMolecule(currentMolecule, moleculeName, moleculeStats, typeOfMolecule, moleculeCharge, moleculeStatusBits, comments, fSource);
				currentMolecule->id = 1 + AVL_Count(crystalMolecules);

				/* add molecule to crystal */
				if(!AVL_Insert(crystalMolecules, currentMolecule))
				{
					printf("ran out of memory\n");
					exit(0);
				}

				numberOfMolecules = AVL_Count(crystalMolecules);
				if(numberOfMolecules % 1000 == 0)
					printf("currentMolecule = %s\tnumber of molecules = %d\n", currentMolecule->name, numberOfMolecules);

				/* if next section can't be determined, read line from file to find next section */
				while(!feof(f) && buffer[0] != '@')
					fscanf(f, " %[^\n]", buffer);

				if(feof(f))
					sectionFlag = -1;
				if(buffer[9] == 'A')
					sectionFlag = 1;
				if(buffer[9] == 'B')
					sectionFlag = 2;
				if(buffer[9] == 'M')
					sectionFlag = 3;

			} /* molecules */

			else if(sectionFlag == 1)
			{ /* atoms */

				/* read line from file */
				fscanf(f, " %[^\n]", buffer);

				while(!feof(f) && buffer[0] != '@')
				{ /* while not eof and line doesn't start with @ */
				
					/* parse atom */
					currentAtom = (atomType *)malloc(sizeof(atomType));
					parseAtom(currentAtom, buffer, crystal);

					currentAtom->id = 1 + AVL_Count(crystalAtoms);

					/* add atom to crystal */
					if(!AVL_Insert(crystalAtoms, currentAtom))
					{
						printf("ran out of memory\n");
						exit(0);
					}

					/* read line from file */
					fscanf(f, " %[^\n]", buffer);
					numberOfAtoms = AVL_Count(crystalAtoms);
					if(numberOfAtoms % 1000 == 0)
						printf("currentAtom = %d\tnumber of atoms = %d\n", currentAtom->id, numberOfAtoms);
			
				} /* while not eof and line doesn't start with @ */

				if(feof(f))
					sectionFlag = -1;
				if(buffer[9] == 'A')
					sectionFlag = 1;
				if(buffer[9] == 'B')
					sectionFlag = 2;
				if(buffer[9] == 'M')
					sectionFlag = 3;
						
			} /* atoms */

		} /* while not eof */

		/* close file */
		fclose(f);

		/* traverse list */
		//AVL_traverse(crystalAtoms, visitAtom);

		/* count crystals */
		numberOfAtoms = AVL_Count(crystalAtoms);
		printf("Number of atoms = %d\n", numberOfAtoms);
		numberOfBonds = AVL_Count(crystalBonds);
		printf("Number of bonds = %d\n", numberOfBonds);
		numberOfMolecules = AVL_Count(crystalMolecules);
		printf("Number of molecules = %d\n", numberOfMolecules);

		/* destroy lists */
		crystalAtoms = AVL_Destroy(crystalAtoms);
		crystalBonds = AVL_Destroy(crystalBonds);
		crystalMolecules = AVL_Destroy(crystalMolecules);
	} /* comment this part out */

	if(0) /* comment this part out */
	{
		/* open list of files */
		fileNames = fopen("listOfFiles.txt", "r");
	
		/* for each filename in the list of files */
		fscanf(fileNames, " %[^\n]", nameOfFile);
		while(!feof(fileNames))
		{
			printf("%s\n", nameOfFile);

			/* document filesource name */
			switch(nameOfFile[0])
			{
				case 'c' :
					fSource = crystal;
					break;
				case 'a' :
					fSource = active;
					break;
				case 'd' :
					fSource = decoy;
					break;
			}
		
			/* open file */
			currentFile = fopen(nameOfFile, "r");

			/* read file */
			readFile(currentFile, &numberOfAtoms, &numberOfBonds, &numberOfMolecules, fSource);

			fscanf(fileNames, " %[^\n]", nameOfFile);
		}

		/* close list of files */
		fclose(fileNames);

		/* give output */
		printf("%d atoms\n", numberOfAtoms);
		printf("%d bonds\n", numberOfBonds);
		printf("%d molecules\n", numberOfMolecules);
	} /* end of commented out part */

	printf("Done.\n");
	return 0;
} /* main */

/****************************************************************************************/

int showFile(FILE *f, int x)
{
  char buffer[1000];
  int lines = x;

  fscanf(f, " %[^\n]", buffer);
  while(!feof(f))
  {
    printf("%s\n", buffer);
    lines++;
    fscanf(f, " %[^\n]", buffer);
  }
  printf("lines = %d\n", lines);
  return lines;

} /* show file */

/****************************************************************************************/

int getAtoms(FILE *f, int *numberOfAtoms, int *numberOfBonds, int *numberOfMolecules, enum eFileSources fSource)
{
	/*
		Input: An open file that is currently in an atom section
		Output: All atoms are read from that section, and the next section type or eof flag is returned
	*/

	char buffer[500];
	atomType a;

	/* read line from file */
	fscanf(f, " %[^\n]", buffer);

	while(!feof(f) && buffer[0] != '@')
	{ /* while not eof and line doesn't start with @ */

		(*numberOfAtoms)++;

		/* parse atom */
		parseAtom(&a, buffer, fSource);

		/* spot checking atoms */
		if((*numberOfAtoms) % 10000 == 0)
		{ /* every 10,000 atoms will be displayed */

			printf("Atom:\tid=%d name=%s\tcoord=(%.3f,%.3f,%.3f)\ttype=%s\tsubstID=%d substName=%s \tcharge=%.3f \tstatusBit=%d\tSource=%d\n",
				a.id, a.name, a.x, a.y, a.z, a.typeOfAtom, a.subst_id, a.subst_name, a.charge, a.statusBit, a.fSource);

		} /* every 10,000 atoms will be displayed */


		/* read line from file */
		fscanf(f, " %[^\n]", buffer);

	} /* while not eof and line doesn't start with @ */

	/* return value should be:
		1	if the file enters an atom section
		2	if the file enters a bond section
		3	if the file enters a molecule section
		-1	if eof
	*/

	/* return */
	if(feof(f))
		return -1;
	if(buffer[9] == 'A')
		return 1;
	if(buffer[9] == 'B')
		return 2;
	if(buffer[9] == 'M')
		return 3;

	/* error */
	return -2;

} /* getAtoms */

/****************************************************************************************/

int getBonds(FILE *f, int *numberOfAtoms, int *numberOfBonds, int *numberOfMolecules, enum eFileSources fSource)
{
	/*
		Input: An open file that is currently in a bond section
		Output: All bonds are read from that section, and the next section type or eof flag is returned
	*/

	char buffer[500];
	bondType b;

	/* read line from file */
	fscanf(f, " %[^\n]", buffer);

	while(!feof(f) && buffer[0] != '@')
	{ /* while not eof and line doesn't start with @ */

		(*numberOfBonds)++;

		/* parse bond */
		parseBond(&b, buffer, fSource);


		/* spot checking bonds */
		if(0 && (*numberOfBonds) % 10000 == 0)
		{ /* every 10,000 bonds will be displayed */

			printf("Bond:\tid=%d originID=%d targetID=%d bondTypeID=%d statusBitID=%d\tsource=%d\n",
				b.id, b.originAtomID, b.targetAtomID, b.typeOfBond, b.statusBit, b.fSource);

		} /* every 10,000 bonds will be displayed */

		/* read line from file */
		fscanf(f, " %[^\n]", buffer);

	} /* while not eof and line doesn't start with @ */

	/* return value should be:
		1	if the file enters an atom section
		2	if the file enters a bond section
		3	if the file enters a molecule section
		-1	if eof
	*/

	/* return */
	if(feof(f))
		return -1;
	if(buffer[9] == 'A')
		return 1;
	if(buffer[9] == 'B')
		return 2;
	if(buffer[9] == 'M')
		return 3;

	/* error */
	return -2;

} /* getBonds */

/****************************************************************************************/

int getMolecules(FILE *f, int *numberOfAtoms, int *numberOfBonds, int *numberOfMolecules, enum eFileSources fSource)
{
	/*
		Input: An open file that is currently in a molecule section
		Output: All molecules are read from that section, and the next section type or eof flag is returned
	*/

	moleculeType m;
	char moleculeName[500], moleculeStats[500], typeOfMolecule[500], moleculeCharge[500], moleculeStatusBits[500], comments[500], buffer[500];

	/* initialize strings */
	strcpy(moleculeName, "");
	strcpy(moleculeStats, "");
	strcpy(typeOfMolecule, "");
	strcpy(moleculeCharge, "");
	strcpy(moleculeStatusBits, "");
	strcpy(comments, "");

	/* read name of molecule */
	fscanf(f, " %[^\n]", moleculeName);

	/* read molecule stats line */
	fscanf(f, " %[^\n]", moleculeStats);

	/* read molecule type */
	fscanf(f, " %[^\n]", typeOfMolecule);

	/* read charge type */
	fscanf(f, " %[^\n]", moleculeCharge);

	/* read line into buffer */
	fscanf(f, " %[^\n]", buffer);

	if(buffer[0] != '@')
	{ /* if buffer doesn't start with @ */

		/* copy buffer into status bits */
		strcpy(moleculeStatusBits, buffer);

		/* read line into buffer */
		fscanf(f, " %[^\n]", buffer);

		if(buffer[0] != '@')
		{ /* if buffer doesn't start with @ */

			/* copy buffer into comments */
			strcpy(comments, buffer);

		} /* if buffer doesn't start with @ */

	} /* if buffer doesn't start with @ */

	(*numberOfMolecules)++;

	/* parse molecule */
	parseMolecule(&m, moleculeName, moleculeStats, typeOfMolecule, moleculeCharge, moleculeStatusBits, comments, fSource);


	/* spot checking molecules */
	if(0 && (*numberOfBonds) % 1000 == 0)
	{ /* every 1,000 bonds will be displayed */

		printf("Molecule:\tname=%s numOfAtoms=%d numOfBonds=%d numOfSubst=%d numOfFeat=%d numOfSets=%d molTypeId=%d chargeTypeId=%d statusBitId=%d comment=%s\tSourse=%d\n",
			m.name, m.numAtoms, m.numBonds, m.numSubStructs, m.numFeatures, m.numSets, m.typeOfMolecule, m.chargeType, m.statusBits, m.comments, m.fSource);

	} /* every 10,000 bonds will be displayed */


	/* if next section can't be determined, read line from file to find next section */
	while(!feof(f) && buffer[0] != '@')
		fscanf(f, " %[^\n]", buffer);

	/* return value should be:
		1	if the file enters an atom section
		2	if the file enters a bond section
		3	if the file enters a molecule section
		-1	if eof
	*/

	/* return */
	if(feof(f))
		return -1;
	if(buffer[9] == 'A')
		return 1;
	if(buffer[9] == 'B')
		return 2;
	if(buffer[9] == 'M')
		return 3;

	/* error */
	return -2;

} /* getMolecules */

/****************************************************************************************/

void readFile(FILE *f, int *numberOfAtoms, int *numberOfBonds, int *numberOfMolecules, enum eFileSources fSource)
{
	/*
		Input: An open file with the pointer at the beginning
		Output: The file is parsed, then closed 
	*/

	char buffer[500];
	int sectionFlag;

	/* read the first line of the file, and determine the section */
	fscanf(f, " %[^\n]", buffer);
	if(buffer[9] == 'A')
		sectionFlag = 1;
	if(buffer[9] == 'B')
		sectionFlag = 2;
	if(buffer[9] == 'M')
		sectionFlag = 3;

	while(!feof(f) && sectionFlag != -1)
	{ /* while not eof */

		/* if file is in atom section, then get atoms */
		if(sectionFlag == 1)
			sectionFlag = getAtoms(f, numberOfAtoms, numberOfBonds, numberOfMolecules, fSource);

		/* if file is in bonds section, then get bonds */
		else if(sectionFlag == 2)
			sectionFlag = getBonds(f, numberOfAtoms, numberOfBonds, numberOfMolecules, fSource);

		/* if file is in molecules section, then get molecules */
		else if(sectionFlag == 3)
			sectionFlag = getMolecules(f, numberOfAtoms, numberOfBonds, numberOfMolecules, fSource);

	} /* while not eof */

	/* close file */
	fclose(f);

} /* readFile */

/****************************************************************************************/

int compareAtoms(atomType *a, atomType *b)
{
	if(a->id == b->id)
	{
		if(strcmp(a->name, b->name) > 0)
			return 1;
		if(strcmp(a->name, b->name) < 0)
			return -1;
		return 0;
	}
	if(a->id > b->id)
		return 1;
	return -1;
} /* compare atoms */

/****************************************************************************************/

int compareBonds(bondType *a, bondType *b)
{
	if(a->id == b->id)
		return 0;
	if(a->id > b->id)
		return 1;
	return -1;
} /* compare bonds */

/****************************************************************************************/

int compareMolecules(moleculeType *a, moleculeType *b)
{
	/*
	if(strcmp(a->name, b->name) > 0)
		return 1;
	if(strcmp(a->name, b->name) < 0)
		return -1;
	return 0;
	*/

	if(a->id == b->id)
		return 0;
	if(a->id > b->id)
		return 1;
	return -1;
} /* compare molecules */

/****************************************************************************************/

void visitAtom(atomType *a)
{
	
	//if(a->id % 10000 == 0)
		printf("Atom:\tid=%d name=%s\tcoord=(%.3f,%.3f,%.3f)\ttype=%s\tsubstID=%d substName=%s \tcharge=%.3f \tstatusBit=%d\tSource=%d\n",
					a->id, a->name, a->x, a->y, a->z, a->typeOfAtom, a->subst_id, a->subst_name, a->charge, a->statusBit, a->fSource);
				
	/*
	printf(".%d.\t", a->id);
	if(a->id % 10 == 0)
		printf("\n");
		*/
}

/****************************************************************************************/