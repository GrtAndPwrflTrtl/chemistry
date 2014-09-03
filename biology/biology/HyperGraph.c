#ifndef _HYPER_GRAPH_H_
#define _HYPER_GRAPH_H_ 1
#include "HyperGraph.h"
#endif

#ifndef _HYPER_EDGE_H_
#define _HYPER_EDGE_H_ 1
#include "HyperEdge.h"
#endif

#ifndef _AVL_TYPE_H_
#define _AVL_TYPE_H_ 1
#include "avl.h"
#endif

#ifndef _HYPER_NODE_
#define _HYPER_NODE_ 1

/*************************************************************************/

int hyperGraphSize(HyperGraphType g)
{ 
	/* returns number of vertices */
	return AVL_Count(g.vertices);
}  /* hyper graph size */

/*************************************************************************/

int createHyperGraph(HyperGraphType g)
{
	/* internal prototype declaration */
	int hyperNodeCompareForCreate(void *a, void *b);

	/* creates a new empty hypergraph */
	g.vertices = AVL_Create(hyperNodeCompareForCreate, destroyHyperNode);
	if(g.vertices == NULL)
		return 0;
	return 1;

} /* create hyper graph */

/*************************************************************************/

int hyperNodeCompareForCreate(void *a, void *b)
{
	/*
		function needed for the create function
	*/
	return 1;
} /* hyperEdgeCompareForCreate */

/*************************************************************************/

void *hyperGraphGetNode(HyperGraphType h, HyperNodeType *x)
{
	/*
		if there exists a node in h that matches x,
			then return pointer to that node
		else
			return NULL

		can also be used as as hasNode.  
	*/

	return AVL_Retrieve(h.vertices, x);
} /* hyper graph get node */

/*************************************************************************/

int hyperGraphAddNode(HyperGraphType h, HyperNodeType *x)
{
	/*
		adds x to h
	*/
	
	return AVL_Insert(h.vertices, x);
} /* hyper graph add node */

/*************************************************************************/

/* Has Local Forward Edge */
int hyperGraphHasLocalForwardEdge(HyperGraphType h, HyperEdgeType *e)
{
	/*
		Searches for e within the nodes in h
	*/

	int foundIt = 0;

	/* internal prototype declaration */
	void visitHasLocalForwardEdge(HyperNodeType *nodePtr, HyperEdgeType *vPtr, int *foundIt);
	
	AVL_traverse_int(h.vertices, visitHasLocalForwardEdge, e, &foundIt);
	return foundIt;
} /* hyper graph has local forward edge */

/*************************************************************************/

void visitHasLocalForwardEdge(HyperNodeType *nodePtr, HyperEdgeType *vPtr, int *foundIt)
{
	/*
		*x is a boolean value that is true, if the edge has been found, and false otherwise 
		*vPtr is pointer to the forward edge that is being sought
	*/

	if(!(*foundIt))
	{ /* if we haven't found it */

		/* search for vPtr in the node's forward edges */
		if(AVL_Retrieve(nodePtr->forwardEdges, vPtr) == NULL)
		{ /* if the node is found */

			/* update foundIt */
			(*foundIt) = 1;

		} /* if the node is found */

	} /* if we haven't found it */

} /* visitHasLocalForwardEdge */

/*************************************************************************/

int hyperGraphHasForwardEdge(HyperGraphType h, int sourceNodes[], int sourceNodeCount, int targetNode)
{
	/*
		searches for an edge defined by a many to one clause mapping within h
		returns true if found, false if not found
	*/

	int foundIt;
	HyperEdgeType e;

	/* internal prototype declaration */
	void checkNodesForward(HyperNodeType *currentNode, HyperEdgeType *edgeToCheck, int *foundIt);

	e.sourceNodes = sourceNodes;
	e.sourceCount = sourceNodeCount;
	e.targetNodes = targetNode;
	
	/* traverse */
	AVL_traverse_int(h.vertices, checkNodesForward, &e, &foundIt);
	
	return foundIt;

} /* hyperGraphHasForwardEdge */

/*************************************************************************/

int hyperGraphGetNodeFromId(HyperGraphType h, int id, HyperNodeType *n)
{
	/*	search for a node with a matching id in h
		if found, let n be that node
		else, let n be null
	*/

	HyperNodeType key;

	key.id = id;
	n = AVL_Retrieve(h.vertices, &key);
	return (n != NULL);

} /* hyper graph get node */

/*************************************************************************/

int hyperGraphAddForwardEdge(HyperGraphType h, int sourceNodes[], int sourceNodeCount, int targetNode, EdgeAnnotationType annot)
{
	int x;
	HyperEdgeType *e;
	HyperNodeType n;

	/* if h has forward edge, success */
	if(hyperGraphHasForwardEdge(h, sourceNodes, sourceNodeCount, targetNode))
		return 1;

	/* add new hyperedge to h */
	e = (HyperEdgeType *) malloc(sizeof (HyperEdgeType));
	if(e == NULL)
		return 0;
	
	if(!HyperEdge_create(sourceNodes, targetNode, annot, sourceNodeCount, e))
	{
		free(e);
		return 0;
	}

	for(x = 0; x < sourceNodeCount; x++)
	{ /* for each sorce node */

		/* find node in graph that matches source node */
		if(hyperGraphGetNodeFromId(h, sourceNodes[x], &n))
		{
			/* add forward edge */
			int AddForwardEdge(HyperNodeType *hNode, HyperEdgeType edge);
			if(!AddForwardEdge(&n, *e))
			{
				HyperEdge_destroy(e);
				free(e);
				return 0;
			}
		}

	} /* for each sorce node */

	/* success */
	HyperEdge_destroy(e);
	free(e);
	return 1;
} /* hyper graph add forward edge */

/*************************************************************************/

void checkNodesForward(HyperNodeType *currentNode, HyperEdgeType *edgeToCheck, int *foundIt)
{
	/* internal prototype declaration */
	void checkEdges(HyperEdgeType *currentEdge, HyperEdgeType *edgeToCheck, int *foundIt);

	/* if not foundIt, check edges within node */
	if(!(*foundIt))
	{
		AVL_traverse_int(currentNode->forwardEdges, checkEdges, edgeToCheck, foundIt);
	}
} /* check nodes */

/*************************************************************************/

void checkEdges(HyperEdgeType *currentEdge, HyperEdgeType *edgeToCheck, int *foundIt)
{
	if(!(*foundIt))
		if(HyperEdge_equals(*currentEdge, *edgeToCheck))
			(*foundIt) = 1;
}

/*************************************************************************/
/*************************************************************************/
/*************************************************************************/

#endif