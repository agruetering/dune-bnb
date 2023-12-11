#include<dune/bnb/bestfirst.h>

using namespace BranchAndBound;
using namespace std;

Subproblem* BestFirst::chooseSubproblem(list<Subproblem*>& subs){
   //initialize the parameters with first subproblem
   Subproblem* bestKnown = *subs.begin();
   double bestKnownValue = bestKnown->getDualBound();
   list<Subproblem*>::iterator bestKnownIt = subs.begin();

   if(optSense == Minimize){
      for(list<Subproblem*>::iterator it = ++(subs.begin()); it != subs.end(); it++) {
	 if((*it)->getDualBound() < bestKnownValue){
	    bestKnownValue = (*it)->getDualBound();
	    bestKnown = *it;
	    bestKnownIt = it;
	 }
      }
   }else if (optSense == Maximize){
      for(list<Subproblem*>::iterator it = ++(subs.begin()); it != subs.end(); it++) {
	 if((*it)->getDualBound() > bestKnownValue){
	    bestKnownValue = (*it)->getDualBound();
	    bestKnown = *it;
	    bestKnownIt = it;
	 }
      }
   }

   subs.erase(bestKnownIt);
   return bestKnown;
}
