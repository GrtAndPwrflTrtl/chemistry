#ifndef PEBBLER_HYPER_NODE_H_ 
#define PEBBLER_HYPER_NODE_H_ 1
#include "PebblerHyperNode.h"
#endif

/*************************************************************************/

int PebblerHyperNodeCreate(PebblerHyperNodeType *n, void *dataPtr, int id)
{
	/* internal prototype declaration */
	int compareEdges(PebblerHyperEdgeType *a, PebblerHyperEdgeType *b);
	void PebblerHyperEdgeDestroy(PebblerHyperEdgeType *e);

	n->edges = AVL_Create(compareEdges, PebblerHyperEdgeDestroy);

	if(n->edges == NULL)
		/* failure */
		return 0;
	
	n->nodes = NULL;
	n->numOfNodes = 0;
	n->pebble = NO_PEBBLE;
	n->pebbled = 0;
	n->id = id;
	n->dataPtr = dataPtr;
	return 1;

} /* create */

/*************************************************************************/

void PebblerHyperNodeDestroy(PebblerHyperNodeType *n)
{
	n->edges = AVL_Destroy(n->edges);
	free(n->dataPtr);
	free(n->nodes);

} /* destroy */

/*************************************************************************/

int PebblerHyperNodeAddEdge(PebblerHyperNodeType *n, PebblerHyperEdgeType *e)
{
	return AVL_Insert(n, e);
} /* add edge */

/*************************************************************************/

int PebblerHyperNodeAddEdgeSourceTarget(PebblerHyperNodeType *n, EdgeAnnotationType annot, int *sourceNodes, int numOfSourceNodes, int targetNode)
{
	PebblerHyperEdgeType *e;

	int PebblerHyperEdgeCreate(PebblerHyperEdgeType *e, int *src, int numOfSrc, int target, EdgeAnnotationType annot);
	void PebblerHyperEdgeDestroy(PebblerHyperEdgeType *e);

	e = (PebblerHyperEdgeType *)malloc(sizeof(PebblerHyperEdgeType));
	if(e == NULL)
		return 0;

	if(!PebblerHyperEdgeCreate(e, sourceNodes, numOfSourceNodes, targetNode, annot))
	{
		free(e);
		return 0;
	}

	if(!PebblerHyperNodeAddEdge(n, e))
	{
		PebblerHyperEdgeDestroy(e);
		free(e);
		return 0;
	}

	return 1;

} /* add edge from sources and target */

/*************************************************************************/

void PebblerHyperNodeToString(PebblerHyperNodeType hEdge, char retString[])
{
	char buffer[500];
	int n;

	/* internal prototype declaration */
	void processEdgeForNodeTraverse(PebblerHyperEdgeType *e, char stringInVisit[]);

	strcpy(retString, " { ");
	sprintf(buffer, "%d", hEdge.id);
	strcat(retString, buffer);
	strcat(retString, ", Pebbled(");
	if(hEdge.pebble == NO_PEBBLE)
		strcat(retString, "NONE), SuccN={");
	if(hEdge.pebble == RED_FORWARD)
		strcat(retString, "RED), SuccN={");
	if(hEdge.pebble == BLUE_BACKWARD)
		strcat(retString, "BLUE), SuccN={");
	if(hEdge.pebble == PURPLE_BOTH)
		strcat(retString, "PURPLE), SuccN={");

	for(n = 0; n < hEdge.numOfNodes; n++)
	{
		sprintf(buffer, "%d", hEdge.nodes[n]);
		strcat(retString, buffer);
		strcat(retString, ",");
	}

	if(hEdge.numOfNodes > 0)
		retString[strlen(retString) - 1] = '\0';

	strcat(retString, "}, SuccE = { ");

	AVL_traverse_string(hEdge.edges, processEdgeForNodeTraverse, buffer);
	strcat(retString, buffer);
	if(AVL_Count(hEdge.edges) > 0)
		retString[strlen(retString) - 2] = '\0';

	strcat(retString, " } }");

} /* to string */

/*************************************************************************/

int compareEdges(PebblerHyperEdgeType *a, PebblerHyperEdgeType *b)
{
	return strcmp(a->annot.justification, b->annot.justification);
} /* compare edges */

/*************************************************************************/

void processEdgeForNodeTraverse(PebblerHyperEdgeType *e, char stringInVisit[])
{
	char buffer[500];

	PebblerHyperEdgeToString(*e, buffer);
	strcat(stringInVisit, buffer);
	strcat(stringInVisit, ", ");

} /* processEdgeForNodeTraverse */

/*************************************************************************/