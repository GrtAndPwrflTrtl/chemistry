#include <stdio.h>
#include "avl.h"

/****************************************************************************************/

int compareInt(int *argu1, int *argu2);
void destroyPtr(void *dataPtr);
int menu(void);
void printInt(int *dataPtr);

/****************************************************************************************/

int main(void)
{
	AVL_TREE *tree;
	int menuChoice, x, *dataPtr;

	tree = AVL_Create(compareInt, destroyPtr);

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
			case 6 :
				printf("Enter number: ");
				scanf("%d", &x);
				dataPtr = malloc(sizeof (int));
				(*dataPtr) = x;
				AVL_AddFirst(tree, dataPtr);
				break;
			case 7 :
				printf("Enter number: ");
				scanf("%d", &x);
				dataPtr = malloc(sizeof (int));
				(*dataPtr) = x;
				AVL_AddLast(tree, dataPtr);
				break;
			case 8 :
				dataPtr = AVL_getFirst(*tree);
				printf("Deleting %d\n", *dataPtr);
				AVL_Delete(tree, dataPtr);
				free(dataPtr);
				break;
			case 9 :
				dataPtr = AVL_getLast(*tree);
				printf("Deleting %d\n", *dataPtr);
				AVL_Delete(tree, dataPtr);
				free(dataPtr);
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
	printf("6:	add left\n");
	printf("7:	add right\n");
	printf("8:	remove left\n");
	printf("9:	remove right\n");
	printf("-1:	exit\n");

	printf("Enter choice: ");
	scanf("%d", &retVal);
	while(retVal != -1 && (retVal < 1 || retVal > 9))
	{
		printf("Not valid choice\n");
		printf("1:	insert\n");
		printf("2:	delete\n");
		printf("3:	count\n");
		printf("4:	search\n");
		printf("5:	traverse\n");
		printf("6:	add left\n");
		printf("7:	add right\n");
		printf("8:	remove left\n");
		printf("9:	remove right\n");
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

void destroyPtr(void *dataPtr)
{
	return;
} /* destroy ptr */
/****************************************************************************************/
