#include <stdio.h>
#include "avl.h"

/****************************************************************************************/

int compareInt(int *argu1, int *argu2);
int menu(void);
void printInt(int *dataPtr);

/****************************************************************************************/

int main(void)
{
	AVL_TREE *tree;
	int menuChoice, x, *dataPtr;

	tree = AVL_Create(compareInt);

	menuChoice = menu();
	while(menuChoice != -1)
	{
		switch(menuChoice)
		{
			case 1 :
				printf("Enter number: ");
				scanf("%d", &x);
				dataPtr = AVL_Retrieve(tree, &x);
				if(dataPtr == NULL)
				{
					dataPtr = malloc(sizeof (int));
					(*dataPtr) = x;
					AVL_Insert(tree, dataPtr);
				}
				else
					printf("%d already exists\n", x);
				break;
			case 2 :
				printf("Enter number: ");
				scanf("%d", &x);
				dataPtr = AVL_Retrieve(tree, &x);
				if(dataPtr == NULL)
					printf("%d does not exitst\n", x);
				else
				{
					dataPtr = malloc(sizeof (int));
					(*dataPtr) = x;
					AVL_Delete(tree, dataPtr);
					free(dataPtr);
				}
				break;
			case 3 :
				printf("There are %d nodes in the tree\n", AVL_Count(tree));
				break;
			case 4 :
				printf("Enter number: ");
				scanf("%d", &x);
				dataPtr = AVL_Retrieve(tree, &x);
				if(dataPtr == NULL)
					printf("%d was NOT found in the tree\n", x);
				else
					printf("%d was found in the tree\n", x);
				break;
			case 5 :
				printf("\n{ ");
				AVL_traverse(tree, printInt);
				printf(" }\n");
				break;
			default :
				printf("error\n");
				break;
		}
		menuChoice = menu();
	}

	tree = AVL_Destroy(tree);

	printf("Done.\n");
	return 0;
} /* main */

/****************************************************************************************/

int menu(void)
{
	int retVal;
	printf("1:	insert\n");
	printf("2:	delete\n");
	printf("3:	count\n");
	printf("4:	search\n");
	printf("5:	traverse\n");
	printf("-1:	exit\n");

	printf("Enter choice: ");
	scanf("%d", &retVal);
	while(retVal != -1 && (retVal < 1 || retVal > 5))
	{
		printf("Not valid choice\n");
		printf("1:	insert\n");
		printf("2:	delete\n");
		printf("3:	count\n");
		printf("4:	search\n");
		printf("5:	traverse\n");
		printf("-1:	exit\n");
		printf("Enter choice: ");
		scanf("%d", &retVal);
	}

	return retVal;

} /* menu */

/****************************************************************************************/

int compareInt(int *argu1, int *argu2)
{
	if(*argu1 == *argu2)
		return 0;
	if(*argu1 < *argu2)
		return -1;
	return 1;
} /* compare int */

/****************************************************************************************/

void printInt(int *dataPtr)
{
	printf(" %d ", (*dataPtr));
} /* print int */

/****************************************************************************************/
/****************************************************************************************/
