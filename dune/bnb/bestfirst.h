#ifndef BESTFIRST_H
#define BESTFIRST_H

#include<dune/bnb/master.h>
#include<dune/bnb/subproblem.h>
#include<dune/bnb/enumeration.h>
#include<list>

namespace BranchAndBound {

   //BestFirst chooses the open subproblem with the best dual bound

   class BestFirst : public EnumerationStrategy {

      private:
	 OptSense optSense;

      public:
	 BestFirst(OptSense optSense_) : optSense(optSense_) {};
	 Subproblem* chooseSubproblem(std::list<Subproblem*>& subs);

   };
}

#endif
