#ifndef _EDGE_H_
#define _EDGE_H_ 1
#include "edge.h"
#endif

/*************************************************************************/

int EdgeCompare(EdgeType a, EdgeType b)
{
	if(a.sourceNode > b.sourceNode)
		return 1;
	if(a.sourceNode < b.sourceNode)
		return -1;
	if(a.targetNode > b.targetNode)
		return 1;
	if(a.targetNode < b.targetNode)
		return -1;
	return 0;
} /* EdgeCompare */

/*************************************************************************/

EdgeType EdgeCreateWithWeight(int sourceNodeId, int targetNodeId, int weight)
{
	EdgeType ret;
	ret.sourceNode = sourceNodeId;
	ret.targetNode = targetNodeId;
	ret.weight = weight;
	return ret;
} /* EdgeCreate with weight */

/*************************************************************************/

EdgeType EdgeCreate(int sourceNodeId, int targetNodeId)
{
	EdgeType ret;
	ret.sourceNode = sourceNodeId;
	ret.targetNode = targetNodeId;
	return ret;
} /* create edge */

/*************************************************************************/

void EdgeDestroy(EdgeType *hEdge)
{
	return;
} /* EdgeDestroy */

/*************************************************************************/

int Edge_equals(EdgeType a, EdgeType b)
{
	return a.targetNode == b.targetNode && a.sourceNode == b.sourceNode;
} /* Edge_equals */

/*************************************************************************/