#ifndef DEPTHFIRST_H
#define DEPTHFIRST_H

#include<dune/bnb/master.h>
#include<dune/bnb/subproblem.h>
#include<dune/bnb/enumeration.h>
#include<list>

namespace BranchAndBound {

   class DepthFirst : public EnumerationStrategy {

      public:
	 Subproblem* chooseSubproblem(std::list<Subproblem*>& subs);

   };
}

#endif
