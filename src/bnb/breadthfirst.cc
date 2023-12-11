#include<dune/bnb/breadthfirst.h>

using namespace BranchAndBound;
using namespace std;

Subproblem* BreadthFirst::chooseSubproblem(list<Subproblem*>& subs){
   //return the last opened subproblem
   Subproblem* sub = subs.front();
   subs.erase(subs.begin());
   return sub;
}
