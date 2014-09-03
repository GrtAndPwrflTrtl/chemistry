#ifndef _HYPER_EDGE_H_
#define _HYPER_EDGE_H_ 1
#include "HyperEdge.h"
#endif

/*************************************************************************/

int HyperEdge_create(int src[], int target, EdgeAnnotationType annot, int sourceCount, HyperEdgeType *hEdge)
{
	/*
		creates hyper edge
	*/
	int x;

	hEdge->sourceNodes = (int *)malloc(sourceCount * sizeof(int));
	if(hEdge->sourceNodes == NULL)
		/* failure */
		return 0;

	/* hEdge->sourceNodes = src;*/
	for(x = 0; x < sourceCount; x++)
		hEdge->sourceNodes[x] = src[x];

	hEdge->targetNodes = target;
	hEdge->sourceCount = sourceCount;
	EdgeAnnotation_createEdgeAnnotation(annot.justification, annot.active, &(hEdge->annotation));

	/* success */
	return 1;

} /* HyperEdge_create */

/*************************************************************************/

int HyperEdge_definesEdge(int antecedent[], int consequent, int antecedentLength, HyperEdgeType hEdge)
{
	/*
		compares antecedent and consequent to the hyper edge
	*/
	int a, x, foundIt;

	for(a = 0; a < antecedentLength; a++)
	{ /* for each a in antecedent */

		/* if a is not in hEdge.sourceNodes, then return false */
		x = 0;
		foundIt = 0;
		while(!foundIt && x < hEdge.sourceCount)
		{
			if(antecedent[a] == hEdge.sourceNodes[x])
				foundIt = 1;
			x++;
		}
		if(!foundIt)
			return 0;

	} /* for each a in antecedent */

	/* At this point in the function, defines edge depends on if target node matches consequent */
	return consequent == hEdge.targetNodes;

} /* HyperEdge_definesEdge */

/*************************************************************************/

void HyperEdge_destroy(HyperEdgeType *hEdge)
{
	/*
		destroys hyper edge and frees dynamically allocated memory
	*/

	free(hEdge->sourceNodes);
	hEdge->sourceNodes = NULL;
	hEdge->sourceCount = -1;

} /* HyperEdge_destroy */

/*************************************************************************/

int HyperEdge_equals(HyperEdgeType a, HyperEdgeType b)
{
	/*
		determine if a = b based on the nodes
	*/

	int aN, bN, foundIt;

	/* check to see if there are the same number of source nodes */
	if(a.sourceCount != b.sourceCount)
		return 0;

	for(aN = 0; aN < a.sourceCount; aN++)
	{ /* for each node aN in a.sourceNodes */


		/* if aN is in b.sourceNodes, then return false */
		bN = 0;
		foundIt = 0;
		while(!foundIt && bN < b.sourceCount)
		{
			if(a.sourceNodes[aN] == b.sourceNodes[bN])
				foundIt = 1;
			bN++;
		}
		if(!foundIt)
			return 0;

	} /* for each node aN in a.sourceNodes */

	/* At this point in the function, equals depends on if the targetNodes match */
	return a.targetNodes == b.targetNodes;

} /* HyperEdge_equals */

/*************************************************************************/

void HyperEdge_toString(HyperEdgeType hEdge, char retString[])
{
	int x;
	char buffer[10];

	/* set retString to " { " */
	strcpy(retString, " { ");
	
	for(x = 0; x < hEdge.sourceCount; x++)
	{ /* for each node in sourceNodes */

		/* append integer to retString, then append a comma */
		//itoa(hEdge.sourceNodes[x], buffer, 10);
		sprintf(buffer, "%d", hEdge.sourceNodes[x]);
		strcat(retString, buffer);
		strcat(retString, ",");

	} /* for each node in sourceNodes */

	/* remove last comma, if there are commas */
	if(hEdge.sourceCount > 0)
		retString[strlen(retString) - 1] = '\0';

	/* append " } " and target node to retString */
	strcat(retString, " } -> ");
	sprintf(buffer, "%d", hEdge.targetNodes);
	strcat(retString, buffer);

} /* HyperEdge_toString */

/*************************************************************************/

int EdgeAnnotation_isActive(EdgeAnnotationType *ea)
{
	return ea->active;
} /* EdgeAnnotation_isActive */

/*************************************************************************/

void EdgeAnnotation_createEdgeAnnotation(char just[], int active, EdgeAnnotationType *ret)
{
	strcpy(ret->justification, just);
	ret->active = active;
} /* EdgeAnnotation_createEdgeAnnotation */

/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/