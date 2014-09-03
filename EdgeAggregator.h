#ifndef _EDGE_AGGREGATOR_GUARD
#define _EDGE_AGGREGATOR_GUARD 1

#include <vector>

#include "Molecule.h"
#include "EdgeAnnotation.h"

//
// An aggregation class of information to pass back from instantiation.
//
class EdgeAggregator
{
  public:
    std::vector<Molecule> antecedent;
    Molecule* consequent;
    EdgeAnnotationT* annotation;
        
    EdgeAggregator(const std::vector<Molecule>& ante, Molecule* c, EdgeAnnotationT* ann)
    {
        for (int a = 0; a < ante.size(); a++)
        {
            antecedent.push_back(ante[a]);
        }

        consequent = c;
        annotation = ann;
    }

    ~EdgeAggregator()
    {
    }
};

#endif
