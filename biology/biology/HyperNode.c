#ifndef _HYPER_NODE_H_
#define _HYPER_NODE_H_ 1
#include "HyperNode.h"
#endif

#ifndef _HYPER_EDGE_H_
#define _HYPER_EDGE_H_ 1
#include "HyperEdge.h"
#endif

#ifndef _HYPER_NODE_
#define _HYPER_NODE_ 1

/*************************************************************************/

int AddBackwardEdge(HyperNodeType *hNode, HyperEdgeType edge)
{
	/*
		Adds a backward edge to hNode
	*/

	HyperEdgeType *edgePtr;

	edgePtr = (HyperEdgeType *)malloc(sizeof(HyperEdgeType));
	if(edgePtr == NULL)
		return 0;

	if(!AVL_Insert(hNode->backwardEdges, edgePtr))
	{
		free(edgePtr);
		return 0;
	}

	return 1;
} /* AddBackwardEdge */

/*************************************************************************/

int AddForwardEdge(HyperNodeType *hNode, HyperEdgeType edge)
{
	/*
		Adds a forward edge to hNode
	*/

	HyperEdgeType *edgePtr;

	edgePtr = (HyperEdgeType *)malloc(sizeof(HyperEdgeType));
	if(edgePtr == NULL)
		return 0;

	if(!AVL_Insert(hNode->forwardEdges, edgePtr))
	{
		free(edgePtr);
		return 0;
	}

	return 1;
} /* AddForwardEdge */

/*************************************************************************/

int createHyperNode(void *data, int id, HyperNodeType *hNode)
{
	/* 
		Allocates dynamic memory for a hyper node, and creates lists of forward edges and backward edges
		memory for data must be already allocated
	*/

	/* internal prototype declaration */
	int compareHyperEdges(void *argu1, void *argu2);
	void destroyEdges(void *v);

	/* create forward edges */
	hNode->forwardEdges = AVL_Create(compareHyperEdges, destroyEdges);
	if(hNode->forwardEdges == NULL)
		return 0;

	/* create backward edges */
	hNode->backwardEdges = AVL_Create(compareHyperEdges, destroyEdges);
	if(hNode->backwardEdges  == NULL)
	{
		hNode->forwardEdges = AVL_Destroy(hNode->forwardEdges);
		return 0;
	}

	/* assign id */
	hNode->id = id;

	/* assign data pointer */
	hNode->dataPtr = data;

} /* createHyperNode */

/*************************************************************************/

void HyperNodeTypeToString(HyperNodeType hNode, char retString[])
{
	/*
		Puts the list of source nodes into the string retString
	*/

	int x;

	/* internal prototype declaration */
	void edgeStringAppend(HyperEdgeType *dataPtr, char stringInVisit[]);

	/* set retString = "\t\t\t\t= { " */
	strcpy(retString, "\t\t\t\t= { ");

	/* append "SuccE = { " to retString */
	strcat(retString, "SuccE = { ");

	/* for each edge in forward edges */
		/* append edge and comma to retString */
	AVL_traverse_string(hNode.forwardEdges, edgeStringAppend, retString);

	/* if there are any forward edges, then remove last comma */
	if(AVL_Count(hNode.forwardEdges) > 0)
		retString[strlen(retString) - 1] = '\0';

	/* append " } }, BackE={" to retString */
	strcat(retString, " } }, BackE={");

	/* for each edge in forward edges */
		/* append edge and comma to retString */
	AVL_traverse_string(hNode.backwardEdges, edgeStringAppend, retString);

	/* if there are any backward edges, then remove last comma */
	if(AVL_Count(hNode.backwardEdges) > 0)
		retString[strlen(retString) - 1] = '\0';

	/* set retString = " } }" */
	strcat(retString, " } }");

} /* HyperNodeTypeToString */

/*************************************************************************/

void destroyHyperNode(HyperNodeType *hNode)
{
	/* 
		Destroys lists of forward edges and backward edges, and frees data in the hypernode.
	*/

	/* destroy forward edges */
	hNode->forwardEdges = AVL_Destroy(hNode->forwardEdges);

	/* destroy backward edges */
	hNode->backwardEdges = AVL_Destroy(hNode->backwardEdges);

	/* destroy data pointer */
	free(hNode->dataPtr);

} /* destroyHyperNode */

/*************************************************************************/

int compareHyperEdges(void *argu1, void *argu2)
{
	/*
		This function is crated for the purpose of storing hyper edges in an AVL tree.
		Since the edges are not in any particular order, the comparison is arbitrary.
		The second argument always falls to the right of the first argument.
	*/
	return 1;

} /* compare hyper edges */

/*************************************************************************/

void edgeStringAppend(HyperEdgeType *dataPtr, char stringInVisit[])
{
	/* 
		for a particular edge in either forward edges or backward edges
		 append edge and comma to retString
	*/

	int x;
	char intBuf[50];

	/* put list of source nodes, followed by colon, followed by target nodes within parentheses */

	strcat(stringInVisit, "(");
	for(x = 0; x < dataPtr->sourceCount - 1; x++)
	{ /* for each source node except the last one, add the node and a comma */

		sprintf(intBuf, "%d", dataPtr->sourceNodes[x]);
//		itoa(dataPtr->sourceNodes[x], intBuf, 10);
		strcat(stringInVisit, intBuf);
		strcat(stringInVisit, ",");

	} /* for each source node except the last one, add the node and a comma */

	/* add last source node followed by colon */
//	itoa(dataPtr->sourceNodes[x], intBuf, 10);
	sprintf(intBuf, "%d", dataPtr->sourceNodes[x]);
	strcat(stringInVisit, intBuf);
	strcat(stringInVisit, " : ");

	/* add target node, followed by close parentheses and a comma */
//	itoa(dataPtr->targetNodes, intBuf, 10);
	sprintf(intBuf, "%d", dataPtr->targetNodes);
	strcat(stringInVisit, intBuf);
	strcat(stringInVisit, "),");

} /* edge string append */

/*************************************************************************/

void destroyEdges(void *v)
{
	return;
}

/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/

#endif