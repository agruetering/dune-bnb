// Copyright © Alexandra Grütering, see file LICENSE.md in module root
// License: GLP-3.0 or later
#ifndef BRANCHING_H
#define BRANCHING_H

#include<dune/bnb/subproblem.h>
#include<list>

namespace BranchAndBound
{
   class BranchingStrategy
   {
      public:
	 virtual std::list<Subproblem*> branch(Subproblem* sub)=0;
   };
}

#endif
