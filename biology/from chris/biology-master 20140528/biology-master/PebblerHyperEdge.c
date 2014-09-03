#ifndef PEBBLER_HYPER_NODE_H_ 
#define PEBBLER_HYPER_NODE_H_ 1
#include "PebblerHyperEdge.h"
#endif

/*************************************************************************/

int PebblerHyperEdgeCreate(PebblerHyperEdgeType *e, int *src, int numOfSrc, int target, EdgeAnnotationType annot)
{
	/*
		src must be dynamically allocated
	*/
	e->sourceNodes = src;
	e->sourcePebbles = NULL;
	e->numOfSourceNodes = numOfSrc;
	e->targetNode = target;
	e->annot.active = annot.active;
	strcpy(e->annot.justification, annot.justification);
	e->numOfSourcePebbles = 0;
	e->pebbleColor = NO_PEBBLE;
	e->pebbles = 0;
	return 1;
} /* create */

/*************************************************************************/

void PebblerHyperEdgeDestroy(PebblerHyperEdgeType *e)
{
	free(e->sourcePebbles);
	free(e->sourceNodes);
} /* destroy */

/*************************************************************************/

int PebblerHyperEdgeIsFullyPebbled(PebblerHyperEdgeType e)
{
	int srcNodeIndex, srcNode, srcPebbleIndex, foundIt;

	for(srcNodeIndex = 0; srcNodeIndex < e.numOfSourceNodes; srcNodeIndex++)
	{ /* for each srcNode in e.sourceNodes */

		srcNode = e.sourceNodes[srcNodeIndex];

		/* if srcNode is not in e.sourcePebbles, then return false */
		foundIt = 0;
		for(srcPebbleIndex = 0; !foundIt && srcPebbleIndex < e.numOfSourcePebbles; srcPebbleIndex++)
			if(srcNode == e.sourcePebbles[srcPebbleIndex])
				foundIt = 1;

		if(!foundIt)
			return 0;

	} /* for each srcNode in e.sourceNodes */

	/* return if source nodes is equal to source pebbles */
	return e.numOfSourceNodes == e.numOfSourcePebbles;

} /* isFullyPebbled */

/*************************************************************************/

void PebblerHyperEdgeSetColor(PebblerHyperEdgeType *e, PebblerColorType color)
{
	e->pebbleColor = color;
} /* set color */

/*************************************************************************/

int PebblerHyperEdgeHasNotBeenVisited(PebblerHyperEdgeType e)
{
	return e.pebbleColor == NO_PEBBLE || !PebblerHyperEdgeIsFullyPebbled(e);
} /* has not been visited */

/*************************************************************************/

int PebblerHyperEdgeEquals(PebblerHyperEdgeType a, PebblerHyperEdgeType b)
{
	int src, srcIndex, foundIt, bSrc;

	for(srcIndex = 0; srcIndex < a.sourceNodes; srcIndex++)
	{ /* for each src in a.sourceNodes */

		foundIt = 0;
		src = a.sourceNodes[srcIndex];

		/* if src is not in b.sourceNodes, then return false */
		for(bSrc = 0; !foundIt && bSrc < b.sourceNodes; bSrc++)
			if(src == b.sourceNodes[bSrc])
				foundIt = 1;

		if(!foundIt)
			return 0;

	} /* for each src in a.sourceNodes */

	return a.targetNode == b.targetNode;
} /* equals */

/*************************************************************************/

void PebblerHyperEdgeToString(PebblerHyperEdgeType e, char retString[])
{
	int x;
	char buffer[10];

	/* set retString to " { " */
	strcpy(retString, " { ");

	for(x = 0; x < e.numOfSourceNodes; x++)
	{ /* for each node in sourceNodes */

		sprintf(buffer, "%d", e.sourceNodes[x]);
		strcat(retString, buffer);
		strcat(retString, ",");
	} /* for each node in sourceNodes */

	/* remove last comma, if there are commas */
	if(e.numOfSourceNodes > 0)
		retString[strlen(retString) - 1] = '\0';

	/* append " } " and target node to retString */
	strcat(retString, " } -> ");
	sprintf(buffer, "%d", e.targetNode);
	strcat(retString, buffer);

} /* to string */

/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
