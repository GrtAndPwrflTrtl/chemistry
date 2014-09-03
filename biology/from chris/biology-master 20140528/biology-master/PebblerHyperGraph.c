#ifndef PEBBLER_HYPER_GRAPH_H_ 
#define PEBBLER_HYPER_GRAPH_H_ 1
#include "PebblerHyperGraph.h"
#endif

/*************************************************************************/

int PebblerHyperGraphCreate(PebblerHyperGraphType *x, AVL_TREE *vertices)
{
	/* internal prototype declaration */
	int PebblerHyperGraphLocalCompare(void *argu1, void *argu2);
	void hyperEdgeDestroyerForCreate(void *x);

	x->forwardPebbledEdges = AVL_Create(PebblerHyperGraphLocalCompare, hyperEdgeDestroyerForCreate);
	if(x->forwardPebbledEdges == NULL)
		return 0;

	x->backwardPebbledEdges = AVL_Create(PebblerHyperGraphLocalCompare, hyperEdgeDestroyerForCreate);
	if(x->backwardPebbledEdges == NULL)
	{
		x->forwardPebbledEdges = AVL_Destroy(x->forwardPebbledEdges);
		free(x->forwardPebbledEdges);
		return 0;
	}

	x->backwardEdgeReachableNodes = AVL_Create(PebblerHyperGraphLocalCompare, hyperEdgeDestroyerForCreate);
	if(x->backwardEdgeReachableNodes == NULL)
	{
		x->forwardPebbledEdges = AVL_Destroy(x->forwardPebbledEdges);
		free(x->forwardPebbledEdges);
		x->backwardPebbledEdges = AVL_Destroy(x->backwardPebbledEdges);
		free(x->backwardPebbledEdges);
		return 0;
	}

	x->vertices = vertices;
	return 1;

} /* create */

/*************************************************************************/

void hyperEdgeDestroyerForCreate(void *x)
{
	PebblerHyperEdgeType *e = x;
	PebblerHyperEdgeDestroy(e);
} /* hyperEdgeDestroyerForCreate */

/*************************************************************************/

int PebblerHyperGraphLocalCompare(void *argu1, void *argu2)
{
	return 1;
} /* PebblerHyperGraphLocalCompareEdges */

/*************************************************************************/

void PebblerHyperGraphDestroy(PebblerHyperGraphType *x)
{
	x->forwardPebbledEdges = AVL_Destroy(x->forwardPebbledEdges);
	free(x->forwardPebbledEdges);
	x->backwardPebbledEdges = AVL_Destroy(x->backwardPebbledEdges);
	free(x->backwardPebbledEdges);
	x->backwardEdgeReachableNodes = AVL_Destroy(x->backwardEdgeReachableNodes);
	free(x->backwardEdgeReachableNodes);
	x->vertices = AVL_Destroy(x->vertices);
	free(x->vertices);
} /* destroy */

/*************************************************************************/

int PebblerHyperGraphClearPebbles(PebblerHyperGraphType *g)
{
	/* internal prototype declaration */
	void clearPebblesInANode(PebblerHyperNodeType *n);

	/* for each node in g->vertices, change node to no pebble */
	AVL_traverse(g->vertices, clearPebblesInANode);
		
} /* ClearPebbles */

/*************************************************************************/

void clearPebblesInANode(PebblerHyperNodeType *n)
{
	/* internal prototype declaration */
	void clearPebblesInEdge(PebblerHyperEdgeType *e);

	/* clear node pebble */
	n->pebbled = 0;
	n->pebble = NO_PEBBLE;

	/* for each pebbled edge, clear pebbles */
	AVL_traverse(n->edges, clearPebblesInEdge);

} /* clearPebblesInANode */

/*************************************************************************/

void clearPebblesInEdge(PebblerHyperEdgeType *e)
{
	e->pebbleColor = NO_PEBBLE;
} /* clearPebblesInEdge */

/*************************************************************************/

int PebblerHyperGraphNumVertices(PebblerHyperGraphType g)
{
	return AVL_Count(g.vertices);
} /* PebblerHyperGraphNumVertices */

/*************************************************************************/

void PebblerHyperGraphGetNode(PebblerHyperGraphType g, int id, PebblerHyperNodeType *retNode)
{
	PebblerHyperNodeType keyNode;

	keyNode.id = id;
	retNode = AVL_Retrieve(g.vertices, &keyNode);

} /* get node */

/*************************************************************************/

int PebblerHyperGraphIsNodePebbledForward(PebblerHyperGraphType g, int id)
{
	PebblerHyperNodeType n;

	PebblerHyperGraphGetNode(g, id, &n);
	return n.pebble == RED_FORWARD;

} /* IsNodePebbledForward */

/*************************************************************************/

int PebblerHyperGraphIsNodePebbledBackward(PebblerHyperGraphType g, int id)
{
	PebblerHyperNodeType n;

	PebblerHyperGraphGetNode(g, id, &n);
	return n.pebble == BLUE_BACKWARD;

} /* IsNodePebbledBackward */

/*************************************************************************/

int PebblerHyperGraphIsNodeNotPebbled(PebblerHyperGraphType g, int id)
{
	PebblerHyperNodeType n;

	PebblerHyperGraphGetNode(g, id, &n);
	return n.pebble == NO_PEBBLE;

} /* IsNodeNotPebbled */

/*************************************************************************/

/* INCOMPLETE, LACKS ESSENTIAL FUNCTIONS */
int PebblerHyperGraphPebble(PebblerHyperGraphType *g, int figure[], int figureLength, int givens[], int givenLength)
{
	HyperNodeType *axiomaticNodes, *reflexiveNodes, *obviousDefinitionNodes;

	/* create lists of hyper nodes */

	/* for each hyper node n in g->hyperGraph */
		/* add hyper node n to one of the three avl trees, depending on what kind of node it is */
		/*!!! this part of the algorithm requires IsAxiomatic, IsReflexive, and IsClearDefinition */

	/* pebble forward */

	/* [ebble backward */

} /* Pebble */  /* INCOMPLETE, LACKS ESSENTIAL FUNCTIONS */

/*************************************************************************/
/*
   AVL_TREE *vertices;
	AVL_TREE *forwardPebbledEdges;
	AVL_TREE *backwardPebbledEdges;
	AVL_TREE *backwardEdgeReachableNodes;

	typedef struct
{
	void *dataPtr;
	int id;
	int pebbled; 
	PebblerColorType pebble;
	int *nodes;
	int numOfNodes;
	AVL_TREE *edges;
}PebblerHyperNodeType;

 
*/

/* PebbleForward */ /* need forward traversal before I can write it */

/*************************************************************************/

/* PebbleBackward */

/*************************************************************************/

/* ForwardTraversal */	/* requires hyper edge multimap */


/*************************************************************************/

/* BackwardPebbleFigure */

/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
