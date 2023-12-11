// Copyright © Alexandra Grütering, see file LICENSE.md in module root
// License: GLP-3.0 or later
#ifndef ENUMERATION_H
#define ENUMERATION_H

#include<dune/bnb/subproblem.h>
#include<list>

namespace BranchAndBound
{
   class EnumerationStrategy
   {
      public:
	 virtual Subproblem* chooseSubproblem(std::list<Subproblem*>& subs)=0;
   };
}

#endif
