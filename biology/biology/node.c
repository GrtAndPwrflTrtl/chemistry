#ifndef _NODE_H_
#define _NODE_H_ 1
#include "node.h"
#endif

#ifndef _EDGE_H_
#define _EDGE_H_ 1
#include "edge.h"
#endif

/*************************************************************************/

int NodeAddBackwardEdge(NodeType *hNode, EdgeType *edge)
{
	return AVL_Insert(hNode->backwardEdges, edge);
} /* add backward edge */

/*************************************************************************/

int NodeAddForwardEdge(NodeType *hNode, EdgeType *edge)
{
	return AVL_Insert(hNode->forwardEdges, edge);
} /* add forward edge */

/*************************************************************************/

int NodeCreate(void *data, int id, NodeType *n, void (*destroyPtr)(void *dataPtr))
{
	AVL_TREE *forwardEdges;
	AVL_TREE *backwardEdges;

	/* internal prototype declaration */
	int NodeCompareForCreate(void *argu1, void *argu2);

	forwardEdges = AVL_Create(NodeCompareForCreate, EdgeDestroy);
	if(forwardEdges == NULL)
		return 0;

	backwardEdges = AVL_Create(NodeCompareForCreate, EdgeDestroy);
	if(backwardEdges == NULL)
	{
		forwardEdges = AVL_Destroy(forwardEdges);
		return 0;
	}

	n->id = id;
	n->dataPtr = data;
	n->forwardEdges = forwardEdges;
	n->backwardEdges = backwardEdges;
	n->destroyPtr = destroyPtr;

	return 1;
} /* create */

/*************************************************************************/

void NodeDestroy(NodeType *n)
{
	n->destroyPtr(n->dataPtr);
	free(n->dataPtr);
	n->forwardEdges = AVL_Destroy(n->forwardEdges);
	n->backwardEdges = AVL_Destroy(n->backwardEdges);
	free(n->forwardEdges);
	free(n->backwardEdges);
} /* destroy */

/*************************************************************************/

int NodeCompareForCreate(void *argu1, void *argu2)
{
	EdgeType *a, *b;

	a = argu1;
	b = argu2;
	return EdgeCompare(*a, *b);
} /* NodeCompareForCreate */

/*************************************************************************/