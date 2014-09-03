#ifndef _AVL_TYPE_H_
#define _AVL_TYPE_H_ 1
#include "avl.h"
#endif

/*************************************************************************/

AVL_TREE *AVL_Create(int (*compare) (void *argu1, void *argu2), void (*destroyPtr)(void *dataPtr))
{
	/*
		Allocates dynamic memory for an AVL tree head node.
	*/

    AVL_TREE *tree;

    tree = (AVL_TREE *) malloc(sizeof (AVL_TREE));
    if(tree)
    {
        tree->root = NULL;
        tree->count = 0;
        tree->compare = compare;
		tree->destroyPtr = destroyPtr;
    }

    return tree;
} /* avl create */

/*************************************************************************/

int AVL_Insert(AVL_TREE *tree, void *dataInPtr)
{
	/*
		inserts new data into AVL tree.
	*/

	/* internal prototype declaration */
    AVL_NODE *AVL__insert(AVL_TREE *tree, AVL_NODE *root, AVL_NODE *newPtr, int *taller);

    AVL_NODE *newPtr;
    int forTaller;

	/* allocates memory for a new pointer */
    newPtr = (AVL_NODE *)malloc(sizeof(AVL_NODE));
    if(!newPtr)
        return FALSE;

	/* a new node is always evenly balanced, and has nulls for left and right child */
    newPtr->bal = EH;
    newPtr->right = NULL;
    newPtr->left = NULL;

	/* put data into new node */
    newPtr->dataPtr = dataInPtr;

	/* run the recursive insert function on the root node of the tree */
    tree->root = AVL__insert(tree, tree->root, newPtr, &forTaller);

	/* increment number of nodes in the tree */
    (tree->count)++;

    return TRUE;
} /* avl insert */

/*************************************************************************/

AVL_NODE *AVL__insert(AVL_TREE *tree, AVL_NODE *root, AVL_NODE *newPtr, int *taller)
{
	/* 
		recursively inserts new node into a tree based on the node 
	*/

	/* internal prototype declaration */
    AVL_NODE *AVL_leftBalance (AVL_NODE *root, int *taller);
    AVL_NODE *AVL_rightBalance (AVL_NODE *root, int *taller);

    if(!root)
    { /* if the current node is null, then insert new node here */

        root = newPtr;
        *taller = TRUE;
        return root;

    } /* if the current node is null, then insert new node here */

    if(tree->compare(newPtr->dataPtr, root->dataPtr) < 0)
    { /* if the new node belongs to the left child of the current node */

		/* recursively add new node to the left child of the current node */
        root->left = AVL__insert(tree, root->left, newPtr, taller);

		/* if the left child is now taller than the right child */
        if(*taller)
            switch(root->bal)
            {
                case LH:	//if the current node was already left heavy, then balance the node
                    root = AVL_leftBalance(root, taller);
                    break;
                case EH:	//if the current node was evenly balanced before, then it is now left heavy
                    root->bal = LH;
                    break;
                case RH:	//if the current node was right heavy before, then it is now evenly balanced, and the tree is not taller
                    root->bal = EH;
                    *taller = FALSE;
                    break;
            } /* switch root balance */

        return root;

    } /* if the new node belongs to the left child of the current node */

    else
    { /* if the new node belongs to the right child of the current node */

		/* recursively add new node to the right child of the current node */
        root->right =  AVL__insert(tree, root->right, newPtr, taller);

		/* if the right child is now taller than the left child */
        if(*taller)
            switch(root->bal)
            {
                case LH:	//if the current node was left heavy before, then it is now evenly balanced, and the tree is not taller
					root->bal = EH;
                    *taller = FALSE;
                    break;
                case EH:	//if the current node was evenly balanced before, then it is now right heavy
                    root->bal = RH;
                    break;
                case RH:	//if the current node was already right heavy, then balance the node
                    root = AVL_rightBalance(root, taller);
                    break;
            } /* switch root balance */

        return root;

    } /* if the new node belongs to the left child of the current node */

    return root;

} /* AVL__insert */

/*************************************************************************/

AVL_NODE *AVL_leftBalance(AVL_NODE *root, int *taller)
{
	/*
		Balance a node that is unbalanced to the left
	*/

	/* internal prototype declaration */
    AVL_NODE *AVL_rotateLeft (AVL_NODE *root);
    AVL_NODE *AVL_rotateRight (AVL_NODE *root);

    AVL_NODE *rightTree;
    AVL_NODE *leftTree;

    leftTree = root->left;
    switch(leftTree->bal)
    {
        case LH:	// if the unbalance is left-left, then do a single rotation
            root->bal = EH;
            leftTree->bal = EH;
            root = AVL_rotateRight(root);
            *taller = FALSE;
            break;
        case EH:
            printf("error in left balance\n");
            exit(100);
        case RH:	// if the unbalance is left-right, then do double rotation

            rightTree = leftTree->right;

			/* relabel the balances */
            switch (rightTree->bal)
            {
                case LH:
                    root->bal = RH;
                    leftTree->bal = EH;
                    break;
                case EH:
                    root->bal = EH;
                    leftTree->bal = EH;
                    break;
                case RH:
                    root->bal = EH;
                    leftTree->bal = LH;
                    break;
            }
            rightTree->bal = EH;

			/* change it from left-right inbalance to a left-left inbalance */
            root->left = AVL_rotateLeft(leftTree);

			/* balance it as you would a left-left */
            root = AVL_rotateRight(root);
            *taller = FALSE;
    }
    return root;

} /* left balance */

/*************************************************************************/

AVL_NODE *AVL_rightBalance(AVL_NODE *root, int *taller)
{
	/*
		Balance a node that is unbalanced to the right
	*/

	/* internal prototype declaration */
    AVL_NODE *AVL_rotateLeft (AVL_NODE *root);
    AVL_NODE *AVL_rotateRight (AVL_NODE *root);

    AVL_NODE *rightTree;
    AVL_NODE *leftTree;

    rightTree = root->right;
    switch(rightTree->bal)
    {
        case RH:	// if the unbalance is right-right, then do a single rotation
            root->bal = EH;
            rightTree->bal = EH;
            root = AVL_rotateLeft(root);
            *taller = FALSE;
            break;
        case EH:
            printf("error in right balance\n");
            exit(100);
        case LH:	// if the unbalance is right-left, then do double rotation

            leftTree = rightTree->left;

			/* relabel the balances */
            switch (leftTree->bal)
            {
                case RH:
                    root->bal = RH;
                    rightTree->bal = EH;
                    break;
                case EH:
                    root->bal = EH;
                    rightTree->bal = EH;
                    break;
                case LH:
                    root->bal = EH;
                    rightTree->bal = LH;
                    break;
            }
            leftTree->bal = EH;

			/* change it from right-left inbalance to a right-right inbalance */
            root->right = AVL_rotateRight(rightTree);

			/* balance it as you would a right-right */
            root = AVL_rotateLeft(root);
            *taller = FALSE;
    }
    return root;

} /* right balance */

/*************************************************************************/

AVL_NODE *AVL_rotateLeft(AVL_NODE *root)
{
	/* exchange pointsers to rotate a node to the left */
    AVL_NODE *tempPtr;

    tempPtr = root->right;
    root->right = tempPtr->left;
    tempPtr->left = root;

    return tempPtr;
} /* rotate left */

/*************************************************************************/

AVL_NODE *AVL_rotateRight(AVL_NODE *root)
{
	/* exchange pointsers to rotate a node to the right */

    AVL_NODE *tempPtr;

    tempPtr = root->left;
    root->left = tempPtr->right;
    tempPtr->right = root;

    return tempPtr;
} /* rotate right */

/*************************************************************************/

int AVL_Delete (AVL_TREE *tree, void *dltKey)
{
	/* delete a node from a tree, and rebalance if necessary */

	/* internal prototype declaration */ 
    AVL_NODE *AVL__delete(AVL_TREE *tree, AVL_NODE *root, void *dltKey, int *shorter, int *success);

    int shorter;
    int success;

    AVL_NODE *newRoot;

	/* delete the target node from the root node using the recursive delete function */
    newRoot = AVL__delete(tree, tree->root, dltKey, &shorter, &success);

	/* if successful, then decriment the count */
    if(success)
    {
        tree->root = newRoot;
        (tree->count)--;
        return TRUE;
    }
    else
        return FALSE;
} /* avl delete */

/*************************************************************************/

AVL_NODE *AVL__delete(AVL_TREE *tree, AVL_NODE *root, void *dltKey, int *shorter, int *success)
{
	/*
		recursively delete a targeted node from the tree if it exists 
	*/

	/* internal prototype declaration */ 
    AVL_NODE *AVL_dltRightBalance(AVL_NODE *root, int *shorter);
    AVL_NODE *AVL_dltLeftBalance(AVL_NODE *root, int *shorter);

    AVL_NODE *dltPtr;
    AVL_NODE *exchPtr;
    AVL_NODE *newRoot;

	/* if the root is empty, then there is nothing to delete */
    if(!root)
    {
        *shorter = FALSE;
        *success = FALSE;
        return NULL;
    }

	/* if the target belongs to the left sub-tree, then recursively search for and delete it from the left subtree */
    if(tree->compare(dltKey, root->dataPtr) < 0)
    {
        root->left = AVL__delete(tree, root->left, dltKey, shorter, success);
        if(*shorter)
            root = AVL_dltRightBalance(root, shorter);
    }

	/* if the target belongs to the right sub-tree, then recursively search for and delete it from the right subtree */
    else if (tree->compare(dltKey, root->dataPtr) > 0)
    {
        root->right = AVL__delete(tree, root->right, dltKey, shorter, success);
        if(*shorter)
            root = AVL_dltLeftBalance(root, shorter);
    }

    else
    { /* WE HAVE A WINNER!!! Delete the current node */

        dltPtr = root;
        if(!root->right)
        { /* if there is no right child, then replace the doomed node with its left child */

            newRoot = root->left;
            *success = TRUE;
            *shorter = TRUE;
			tree->destroyPtr(dltPtr);
            free(dltPtr);
            return newRoot;

        } /* if there is no right child, then replace the doomed node with its left child */
        else
            if(!root->left)
            { /* if there is no left child, then replace the doomed node with its right child */

                newRoot = root->right;
                *success = TRUE;
                *shorter = TRUE;
				tree->destroyPtr(dltPtr);
                free(dltPtr);
                return newRoot;

            } /* if there is no left child, then replace the doomed node with its right child */
            else
            { /* the current node has 2 children */

				/* find the right-most descendant of the left child */
                exchPtr = root->left;
                while(exchPtr->right)
                    exchPtr = exchPtr->right;

				/* swap current node with the node found in the previous step */
                root->dataPtr = exchPtr->dataPtr;

				/* recursively search for and delete the node from the left subtree */
                root->left = AVL__delete(tree, root->left, exchPtr->dataPtr, shorter, success);
                if(*shorter)
                    root = AVL_dltRightBalance(root, shorter);

            } /* the current node has 2 children */

    } /* WE HAVE A WINNER!!! Delete the current node */

    return root;

} /* AVL__delete */

/*************************************************************************/

AVL_NODE *AVL_dltRightBalance(AVL_NODE *root, int *shorter)
{
	/*
		A deletion has occured on the left subtree, and now must be rebalanced.
	*/

	/* internal prototype declaration */ 
    AVL_NODE *AVL_rotateLeft(AVL_NODE *root);
    AVL_NODE *AVL_rotateRight(AVL_NODE *root);

    AVL_NODE *rightTree;
    AVL_NODE *leftTree;

    switch(root->bal)
    { /* switch root balance */

        case LH:	// if it was left-heavy before the deletion, then the deletion caused it to be evenly balanced
            root->bal = EH;
            break;
        case EH:	// if it was evenly balanced before the deletion, then the deletion caused it to be right-heavy
            root->bal = RH;
            *shorter = FALSE;
            break;
        case RH:	// if it was right-heavy before the deletion, then it is unbalanced to the right

            rightTree = root->right;

            if(rightTree->bal == LH)
            { /* if the right tree is left-heavy */

				/* relabel the balances */
                leftTree = rightTree->left;
                switch(leftTree->bal)
                {
                    case LH:
                        rightTree->bal = RH;
                        root->bal = EH;
                        break;
                    case EH:
                        root->bal = EH;
                        rightTree->bal = EH;
                        break;
                    case RH:
                        root->bal = LH;
                        rightTree->bal = EH;
                        break;
                }
                leftTree->bal = EH;

				/* do a double-rotation */
                root->right = AVL_rotateRight(rightTree);
                root = AVL_rotateLeft(root);

            } /* if the right tree is left-heavy */
            else
            { /* if the right tree is not left-heavy */

				/* relabel the balances */
                switch(rightTree->bal)
                {
                    case LH:
                    case RH:
                        root->bal = EH;
                        rightTree->bal = EH;
                        break;
                    case EH:
                        root->bal = RH;
                        rightTree->bal = LH;
                        *shorter = FALSE;
                        break;
                }

				/* do a single-rotation */
                root = AVL_rotateLeft(root);

            } /* if the right tree is not left-heavy */
    } /* switch root balance */

    return root;

} /* AVL_dltRightBalance */

/*************************************************************************/

AVL_NODE *AVL_dltLeftBalance(AVL_NODE *root, int *shorter)
{
	/*
		A deletion has occured on the right subtree, and now must be rebalanced.
	*/

	/* internal prototype declaration */ 
    AVL_NODE *AVL_rotateLeft(AVL_NODE *root);
    AVL_NODE *AVL_rotateRight(AVL_NODE *root);

    AVL_NODE *rightTree;
    AVL_NODE *leftTree;

    switch(root->bal)
    { /* switch root balance */

        case RH:	// if it was right-heavy before the deletion, then the deletion caused it to be evenly balanced
            root->bal = EH;
            break;
        case EH:	// if it was evenly balanced before the deletion, then the deletion caused it to be left-heavy
            root->bal = LH;
            *shorter = FALSE;
            break;
        case LH:	// if it was left-heavy before the deletion, then it is unbalanced to the left

            leftTree = root->left;

            if(leftTree->bal == LH)
            { /* if the left tree is right-heavy */

				/* relabel the balances */
                rightTree = leftTree->right;
                switch(rightTree->bal)
                {
                    case RH:
                        leftTree->bal = LH;
                        root->bal = EH;
                        break;
                    case EH:
                        root->bal = EH;
                        leftTree->bal = EH;
                        break;
                    case LH:
                        root->bal = RH;
                        leftTree->bal = EH;
                        break;
                }
                rightTree->bal = EH;

				/* do a double-rotation */	
                root->left = AVL_rotateLeft(leftTree);
                root = AVL_rotateRight(root);

            } /* if the left tree is right-heavy */
            else
            { /* if the left tree is not right-heavy */

				/* relabel the balances */
                switch(leftTree->bal)
                {
                    case LH:
                    case RH:
                        root->bal = EH;
                        leftTree->bal = EH;
                        break;
                    case EH:
                        root->bal = LH;
                        leftTree->bal = RH;
                        *shorter = FALSE;
                        break;
                }

				/* do a single-rotation */
                root = AVL_rotateRight(root);

            } /* if the left tree is not right-heavy */
    } /* switch root balance */

    return root;

} /* AVL_dltRightBalance */

/*************************************************************************/

void *AVL_Retrieve(AVL_TREE *tree, void *keyPtr)
{
	/* find a node in the tree */

	/* internal prototype declaration */
    void *AVL__retrieve(AVL_TREE *tree, void *keyPtr, AVL_NODE *root);

    if(tree->root)
        return AVL__retrieve(tree, keyPtr, tree->root);
    else
        return NULL;
} /* avl retrieve */

/*************************************************************************/

void *AVL__retrieve(AVL_TREE *tree, void *keyPtr, AVL_NODE *root)
{
	/* recursively find a node in the tree */

    if(root)
    {
        if(tree->compare(keyPtr, root->dataPtr) < 0)
            return AVL__retrieve(tree, keyPtr, root->left);
        else if(tree->compare(keyPtr, root->dataPtr) > 0)
            return AVL__retrieve(tree, keyPtr, root->right);
        else
            return(root->dataPtr);
    }
    else
        return NULL;
} /* AVL__retrieve */

/*************************************************************************/

void AVL_traverse(AVL_TREE *tree, void(*process)(void *dataPtr))
{
	/* visit each node in the tree with a function that the user passes in */

	/* internal prototype declaration */
    void AVL__traversal(AVL_NODE *root, void(*process) (void *dataPtr));

    AVL__traversal(tree->root, process);
    return;
} /* avl traverse */

/*************************************************************************/

void AVL__traversal(AVL_NODE *root, void(*process) (void *dataPtr))
{
	/* recursively visit each node in the tree with a function that the user passes in */

    if(root)
    {
        AVL__traversal(root->left, process);
        process(root->dataPtr);
        AVL__traversal(root->right, process);
    }
} /* AVL__traversal */

/*************************************************************************/

int AVL_Empty(AVL_TREE *tree)
{
    return(tree->count == 0);
} /* avl empty */

/*************************************************************************/

int AVL_Full(AVL_TREE *node)
{
    AVL_NODE *newPtr;

    newPtr = (AVL_NODE *)malloc(sizeof (AVL_NODE));
    if(newPtr)
    {
        free(newPtr);
        return FALSE;
    }
    else
        return TRUE;
} /* avl full */

/*************************************************************************/

int AVL_Count(AVL_TREE *tree)
{
    return tree->count;
} /* avl count */

/*************************************************************************/

AVL_TREE *AVL_Destroy(AVL_TREE *tree)
{
	/* frees memory taken up by tree */

	/* internal prototype declaration */
    void AVL__destroyAVL(AVL_TREE *tree, AVL_NODE *root);

    if(tree)
        AVL__destroyAVL(tree, tree->root);
    free(tree);
    return NULL;
} /* avl destroy */

/*************************************************************************/

void AVL__destroyAVL(AVL_TREE *tree, AVL_NODE *root)
{
	/* recursively frees memory taken up by tree */

    if(root)
    {
		/* free the left tree */
        AVL__destroyAVL(tree, root->left);

		/* free the data in current node */
		tree->destroyPtr(root->dataPtr);
        free(root->dataPtr);

		/* free the right tree */
        AVL__destroyAVL(tree, root->right);

		/* destroy the root */
        free(root);
    }
} /* _destroy avl */

/*************************************************************************/

void AVL_traverse_string(AVL_TREE *tree, void(*process)(void *dataPtr, char stringInVisit[]), char stringFromOut[])
{
	/* visit each node in the tree with a function that the user passes in */

	/* internal prototype declaration */
    void AVL__traversal_string(AVL_NODE *root, void(*process) (void *dataPtr, char stringInVisit[]), char stringFromOut[]);

    AVL__traversal_string(tree->root, process, stringFromOut);
    return;
} /* avl traverse with string */

/*************************************************************************/

void AVL__traversal_string(AVL_NODE *root, void(*process) (void *dataPtr, char stringInVisit[]), char stringFromOut[])
{
	/* recursively visit each node in the tree with a function that the user passes in which can manipulate a string */

    if(root)
    {
        AVL__traversal_string(root->left, process, stringFromOut);
        process(root->dataPtr, stringFromOut);
        AVL__traversal_string(root->right, process, stringFromOut);
    }
} /* AVL__traversal_string */

/*************************************************************************/

void *AVL_getAny(AVL_TREE tree)
{
	/* if tree is empty, return NULL */
	if(tree.root == NULL)
		return NULL;

	/* return data in root */
	return (tree.root->dataPtr);

} /* get any */

/*************************************************************************/

void *AVL_getFirst(AVL_TREE tree)
{
	 AVL_NODE *root;

	/* if tree is empty, return NULl */
	if(tree.root == NULL)
		return NULL;

	 /* find data in left-most node */
	 root = tree.root;
	 while(root->left != NULL)
		 root = root->left;

	 return (root->dataPtr);

} /* get first */

/*************************************************************************/

void *AVL_getLast(AVL_TREE tree)
{
	 AVL_NODE *root;

	/* if tree is empty, return NULl */
	if(tree.root == NULL)
		return NULL;

	 /* find data in left-most node */
	 root = tree.root;
	 while(root->right != NULL)
		 root = root->right;

	 return (root->dataPtr);

} /* get last */

/*************************************************************************/

void AVL_traverse_int(AVL_TREE *tree, void(*process)(void *dataPtr, void *auxPtr, int *x), void *auxPtr, int *x)
{
	/* visit each node in the tree with a function that the user passes in */

	/* internal prototype declaration */
    void AVL__traversal_int(AVL_NODE *root, void(*process)(void *dataPtr, void *auxPtr, int x), void *auxPtr, int *x);

    AVL__traversal_int(tree->root, process, auxPtr, x);
    return;
} /* avl traverse with string */

/*************************************************************************/

void AVL__traversal_int(AVL_NODE *root, void(*process)(void *dataPtr, void *auxPtr, int x), void *auxPtr, int *x)
{
	/* recursively visit each node in the tree with a function that the user passes in which can manipulate a string */

    if(root)
    {
        AVL__traversal_int(root->left, process, auxPtr, x);
        process(root->dataPtr, auxPtr, x);
        AVL__traversal_int(root->right, process, auxPtr, x);
    }
} /* AVL__traversal_string */

/*************************************************************************/

int AVL_AddFirst(AVL_TREE *tree, void *dataPtr)
{
	AVL_NODE *newPtr, *curPtr, *prevPtr;
    int taller, hadToBalance = 0;

	/* internal prototype declaration */
    AVL_NODE *AVL_leftBalance (AVL_NODE *root, int *taller);

	/* allocates memory for a new pointer */
	newPtr = (AVL_NODE *)malloc(sizeof(AVL_NODE));
    if(!newPtr)
        return FALSE;

	/* a new node is always evenly balanced, and has nulls for left and right child */
    newPtr->bal = EH;
    newPtr->right = NULL;
    newPtr->left = NULL;

	/* put data into new node */
    newPtr->dataPtr = dataPtr;

	if(!(tree->root))
	{ /* if the current node is null, then insert new node here */

		tree->root = newPtr;
        taller = TRUE;
		return TRUE;

	} /* if the current node is null, then insert new node here */

	/* find place for new node */
	curPtr = tree->root;
	while(curPtr->left != NULL)
		curPtr = curPtr->left;

	/* insert node */
	curPtr->left = newPtr;

	while(!hadToBalance && curPtr != tree->root)
	{ /* work your way back up to the root or balance done */

		if(taller)
			switch(curPtr->bal)
			{
				case LH:	//if the current node was already left heavy, then balance the node
					hadToBalance = 1;
					curPtr = AVL_leftBalance(curPtr, &taller);
					break;
				case EH:	//if the current node was evenly balanced before, then it is now left heavy
					curPtr->bal = LH;
					break;
				case RH: //if the current node was right heavy before, then it is now evenly balanced, and the tree is not taller
					curPtr->bal = EH;
					taller = FALSE;
					break;
			} /* switch current balance */

		/* find parent of curPtr */
		prevPtr = tree->root;
		while(prevPtr->left != curPtr)
			prevPtr = prevPtr->left;
		curPtr = prevPtr;

	} /* work your way back up to the root or balance done */

	/* balance root, if necessary */
	if(taller && !hadToBalance)
		switch(curPtr->bal)
		{
			case LH:	//if the current node was already left heavy, then balance the node
				curPtr = AVL_leftBalance(curPtr, &taller);
				break;
			case EH:	//if the current node was evenly balanced before, then it is now left heavy
				curPtr->bal = LH;
				break;
			case RH: //if the current node was right heavy before, then it is now evenly balanced, and the tree is not taller
				curPtr->bal = EH;
				taller = FALSE;
				break;
		} /* switch current balance */

	/* increment number of nodes in the tree */
    (tree->count)++;

    return TRUE;

} /* AVL_AddLeft */

/*************************************************************************/

int AVL_AddLast(AVL_TREE *tree, void *dataPtr)
{
	AVL_NODE *newPtr, *curPtr, *prevPtr;
    int taller, hadToBalance = 0;

	/* internal prototype declaration */
	AVL_NODE *AVL_rightBalance (AVL_NODE *root, int *taller);

	/* allocates memory for a new pointer */
	newPtr = (AVL_NODE *)malloc(sizeof(AVL_NODE));
    if(!newPtr)
        return FALSE;

	/* a new node is always evenly balanced, and has nulls for left and right child */
    newPtr->bal = EH;
    newPtr->right = NULL;
    newPtr->left = NULL;

	/* put data into new node */
    newPtr->dataPtr = dataPtr;

	if(!(tree->root))
	{ /* if the current node is null, then insert new node here */

		tree->root = newPtr;
        taller = TRUE;
		return TRUE;

	} /* if the current node is null, then insert new node here */

	/* find place for new node */
	curPtr = tree->root;
	while(curPtr->right != NULL)
		curPtr = curPtr->right;

	/* insert node */
	curPtr->right = newPtr;

	while(!hadToBalance && curPtr != tree->root)
	{ /* work your way back up to the root or balance done */

		if(taller)
			switch(curPtr->bal)
			{
				case LH:	//if the current node was left heavy before, then it is now evenly balanced, and the tree is not taller
					curPtr->bal = EH;
                    taller = FALSE;
                    break;
                case EH:	//if the current node was evenly balanced before, then it is now right heavy
                    curPtr->bal = RH;
                    break;
                case RH:	//if the current node was already right heavy, then balance the node
					hadToBalance = 1;
                    curPtr = AVL_rightBalance(curPtr, &taller);
                    break;
			} /* switch current balance */

		/* find parent of curPtr */
		prevPtr = tree->root;
		while(prevPtr->right != curPtr)
			prevPtr = prevPtr->right;
		curPtr = prevPtr;

	} /* work your way back up to the root or balance done */

	/* balance root, if necessary */
	if(taller && !hadToBalance)
		switch(curPtr->bal)
		{
			case LH:	//if the current node was already left heavy, then balance the node
				curPtr = AVL_leftBalance(curPtr, &taller);
				break;
			case EH:	//if the current node was evenly balanced before, then it is now left heavy
				curPtr->bal = LH;
				break;
			case RH: //if the current node was right heavy before, then it is now evenly balanced, and the tree is not taller
				curPtr->bal = EH;
				taller = FALSE;
				break;
		} /* switch current balance */

	/* increment number of nodes in the tree */
    (tree->count)++;

    return TRUE;

} /* AVL_AddRight */

/*************************************************************************/

int AVL_RemoveFirst(AVL_TREE *tree)
{
	AVL_NODE *toDelete, *curPtr, *prevPtr, *rightTree, *leftTree;
    int shorter, hadToBalance = 0;

	/* internal prototype declaration */
	AVL_NODE *AVL_dltRightBalance(AVL_NODE *root, int *shorter);
	AVL_NODE *AVL_rotateLeft(AVL_NODE *root);
    AVL_NODE *AVL_rotateRight(AVL_NODE *root);
    
	/* find left most node, set to toDelete node */
	toDelete = tree->root;
	while(toDelete->left != NULL)
		toDelete = toDelete->left;

	/* destroy and free data in node */
	tree->destroyPtr(toDelete->dataPtr);
	free(toDelete->dataPtr);

	if(tree->root != toDelete)
	{ /* if not at root */

		/* find prev node */
		prevPtr = tree->root;
		while(prevPtr->left != toDelete)
			prevPtr = prevPtr->left;

		/* put toDelete node's right child into prev node's left child */
		prevPtr->left = toDelete->right;

		/* balance tree */
		shorter = TRUE;
		curPtr = prevPtr;
		while(!hadToBalance && curPtr != tree->root)
		{ /* work your way back up to the root or until balance done */

			if(shorter)
			{ /* if shorter */
				switch(curPtr->bal)
				{ /* switch root balance */

					case LH:	// if it was left-heavy before the deletion, then the deletion caused it to be evenly balanced
						curPtr->bal = EH;
						break;
					case EH:	// if it was evenly balanced before the deletion, then the deletion caused it to be right-heavy
						curPtr->bal = RH;
						shorter = FALSE;
						break;
					 case RH:	// if it was right-heavy before the deletion, then it is unbalanced to the right

						 rightTree = curPtr->right;

						 if(rightTree->bal == LH)
						 { /* if the right tree is left heavy */

							 /* relabel the balances */
							 leftTree = rightTree->left;
							 switch(leftTree->bal)
							 {
								 case LH:
									rightTree->bal = RH;
									root->bal = EH;
									break;
								case EH:
									root->bal = EH;
									rightTree->bal = EH;
									break;
								case RH:
									root->bal = LH;
									rightTree->bal = EH;
									break;
							 }
							 leftTree->bal = EH;

							 /* do a double-rotation */
							root->right = AVL_rotateRight(rightTree);
							root = AVL_rotateLeft(curPtr);

						 } /* if the right tree is left heavy */
						 else
						 { /* if the right tree is not left heavy */

							/* relabel the balances */
							switch(rightTree->bal)
							{
								case LH:
								case RH:
									curPtr->bal = EH;
									rightTree->bal = EH;
									break;
								case EH:
									curPtr->bal = RH;
									rightTree->bal = LH;
									shorter = FALSE;
									break;
							}

							/* do a single-rotation */
							curPtr = AVL_rotateLeft(curPtr);

						} /* if the right tree is not left heavy */
				} /* switch root balance */
			} /* if shorter */

			/* find current node's parent */
			/***/

		} /* work your way back up to the root or until balance done */

		/* balance root, if necessary */

	} /* if not at root */
	else
	{ /* else, if at root */
		/* put root's right child into root */
	} /* else, if at root */
	/* destroy current node */
	/* decrement number of nodes in tree */

} /* remove first */

/*************************************************************************/

/* remove last */
int AVL_RemoveLast(AVL_TREE *tree)
{
	AVL_NODE *curPtr, *prevPtr;
    int taller, hadToBalance = 0;

	/* internal prototype declaration */
	AVL_NODE *AVL_dltLeftBalance(AVL_NODE *root, int *shorter);

	/* find right most node, set to current node */
	/* destroy and free data in node */
	/* if not at root */
		/* find prev node */
		/* put current node's left child into prev node's right child */
		/* balance tree */
		/* balance root, if necessary */
	/* else, if at root */
		/* put root's left child into root */
		
	/* destroy current node */
	/* decrement number of nodes in tree */
} /* remove last */

/*************************************************************************/