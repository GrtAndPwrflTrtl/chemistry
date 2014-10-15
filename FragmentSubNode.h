#ifndef _FRAGMENT_SUB_NODE_GUARD
#define _FRAGMENT_SUB_NODE_GUARD 1


class FragmentGraphNode;


class FragmentSubNode
{
  public:
    FragmentSubNode();
    FragmentSubNode(unsigned int id, FragmentGraphNode* parent);

    virtual unsigned int degree() = 0;
    virtual FragmentSubNode* copy() const = 0;
    virtual void addConnection(FragmentSubNode* connector) = 0;
    virtual bool IsIsomorphicTo(FragmentSubNode* that) = 0;

    inline unsigned int getSubNodeID() const { return uniqueSubnodeID; }

    std::string toString() const;
    friend std::ostream& operator<< (std::ostream& os, const FragmentSubNode& node);

    FragmentGraphNode* getParentNode() const { return parentFragmentNode; }

  protected:
    unsigned int uniqueSubnodeID;
    FragmentGraphNode* parentFragmentNode;
};

#endif
