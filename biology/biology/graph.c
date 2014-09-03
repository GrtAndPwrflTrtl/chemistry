#ifndef _GRAPH_H_
#define _GRAPH_H_ 1
#include "graph.h"
#endif

/*************************************************************************/

int GraphCreate(GraphType g, void (*nodeDestroyer)(void *dataPtr), int (*compare) (void *argu1, void *argu2))
{
	AVL_TREE *nodes;
	AVL_TREE *edges;
	
	/* internal prototype declaration */
	int GraphCompareForCreate(void *argu1, void *argu2);

	nodes = AVL_Create(compare, nodeDestroyer);
	if(nodes == NULL)
		return 0;
	edges = AVL_Create(GraphCompareForCreate, EdgeDestroy);
	if(edges == NULL)
	{
		edges = AVL_Destroy(edges);
		return 0;
	}

	g.nodes = nodes;
	g.edges = edges;
	g.nodeDestroyer = nodeDestroyer;
	g.compare = compare;
	return 1;
} /* create */

/*************************************************************************/

void GraphDestroy(GraphType g)
{
	/* destroy and free edges */
	g.nodes = AVL_Destroy(g.nodes);
	free(g.nodes);

	/* destroy and free nodes */
	g.edges = AVL_Destroy(g.edges);
	free(g.edges);

} /* destroy */

/*************************************************************************/

int GraphAddNode(GraphType h, NodeType *x)
{
	return AVL_Insert(h.nodes, x);
} /* add node */

/*************************************************************************/

int GraphAddEdge(GraphType h, int sourceNodeID, int targetNodeID)
{
	EdgeType *e;
	NodeType *nSourcePtr, *nTargetPtr, nSource, nTarget;

	/* find source node */
	nSource.id = sourceNodeID;
	nSourcePtr = &nSource;
	nSourcePtr = AVL_Retrieve(h.nodes, nSourcePtr);
	if(nSourcePtr == NULL)
		return 0;

	/* find targetNode */
	nTarget.id = targetNodeID;
	nTargetPtr = &nTarget;
	nTargetPtr = AVL_Retrieve(h.nodes, nTargetPtr);
	if(nTargetPtr == NULL)
		return 0;

	/* make new edge */
	e = (EdgeType *)malloc(sizeof(EdgeType));
	if(e == NULL)
		return 0;
	e->targetNode = targetNodeID;
	e->sourceNode = sourceNodeID;

	/* add forward edge to source node */
	if(!AVL_Insert(nSourcePtr->forwardEdges, e))
	{
		free(e);
		return 0;
	}

	/* add back edge to target node */
	if(!AVL_Insert(nTargetPtr->backwardEdges, e))
	{
		AVL_Delete(nSourcePtr->forwardEdges, e);
		free(e);
		return 0;
	}

	/* add edge to h.edges */
	if(!AVL_Insert(h.edges, e))
	{
		AVL_Delete(nTargetPtr->backwardEdges, e);
		AVL_Delete(nSourcePtr->forwardEdges, e);
		free(e);
		return 0;
	}

	return 1;
} /* add edge */

/*************************************************************************/

void *GraphGetNode(GraphType h, NodeType *x)
{
	return AVL_Retrieve(h.nodes, x);
} /* get node */

/*************************************************************************/

int GraphGetNodeFromId(GraphType h, int id, NodeType *n)
{
	NodeType *x, y;

	/* make fake node */
	x = &y;
	y.id = id;

	/* run avl retrieve on nodes using fake node */
	x = AVL_Retrieve(h.nodes, x);

	/* set n to fake node */
	if(x == NULL)
		return 0;
	n = x;

	return 1;

} /* get node from id */

/*************************************************************************/

int GraphHasEdge(GraphType h, int sourceNodeID, int targetNodeID)
{
	EdgeType e, *ePtr;

	/* make fake edge */
	e.targetNode = targetNodeID;
	e.sourceNode = sourceNodeID;
	ePtr = &e;

	/* run avl retrieve on edges using fake edge */
	ePtr = AVL_Retrieve(h.edges, ePtr);

	/* return boolean of found */
	return ePtr != NULL;

} /* has edge */

/*************************************************************************/

int GraphGetNumOfEdges(GraphType h)
{
	return AVL_Count(h.edges);
} /* get num of edges */

/*************************************************************************/

int GraphGetNumOfNodes(GraphType h)
{
	return AVL_Count(h.nodes);
} /* get num of nodes */

/*************************************************************************/

int GraphCompareForCreate(void *argu1, void *argu2)
{
	EdgeType *a, *b;

	a = argu1;
	b = argu2;
	return EdgeCompare(*a, *b);
} /* compare for create */

/*************************************************************************/
/*
void GraphBreadthFirstTraverse(GraphType g, void(*process)(void *dataPtr))
{
} 
*/
/*************************************************************************/
/*
void GraphDepthFirstTraverse(GraphType g, void(*process)(void *dataPtr))
{
}
*/
/*************************************************************************/