#ifndef BREADTHFIRST_H
#define BREADTHFIRST_H

#include<dune/bnb/master.h>
#include<dune/bnb/subproblem.h>
#include<dune/bnb/enumeration.h>
#include<list>

namespace BranchAndBound {

   class BreadthFirst : public EnumerationStrategy {

      public:
	 Subproblem* chooseSubproblem(std::list<Subproblem*>& subs);

   };
}

#endif
