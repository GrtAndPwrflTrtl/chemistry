#ifndef _LINKER_FRAGMENT_SUB_NODE_GUARD
#define _LINKER_FRAGMENT_SUB_NODE_GUARD 1


#include <vector>


#include "FragmentSubNode.h"


class LinkerFragmentSubNode : public FragmentSubNode
{
  public:
    LinkerFragmentSubNode();
    LinkerFragmentSubNode(unsigned int id, FragmentGraphNode* parent);

    FragmentSubNode* copy() const;

    void addConnection(FragmentSubNode* connector);

    unsigned int degree() { return connections.size(); }

    bool IsIsomorphicTo(FragmentSubNode* that);

    std::string toString() const;

  private:
    std::vector<FragmentSubNode*> connections;
};

#endif
