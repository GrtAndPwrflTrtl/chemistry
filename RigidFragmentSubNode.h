#ifndef _RIGID_FRAGMENT_SUB_NODE_GUARD
#define _RIGID_FRAGMENT_SUB_NODE_GUARD 1


#include <vector>


#include "FragmentSubNode.h"


class RigidFragmentSubNode : public FragmentSubNode
{
  public:
    RigidFragmentSubNode();
    RigidFragmentSubNode(unsigned int id, FragmentGraphNode* parent);

    FragmentSubNode* copy() const;

    void addConnection(FragmentSubNode* connector);

    unsigned int degree() { return connection == 0 ? 0 : 1; }

    bool IsIsomorphicTo(FragmentSubNode* that);

    std::string toString() const;

  private:
    FragmentSubNode* connection;
};

#endif
