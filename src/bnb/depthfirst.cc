#include<dune/bnb/depthfirst.h>

using namespace BranchAndBound;
using namespace std;

Subproblem* DepthFirst::chooseSubproblem(list<Subproblem*>& subs){
   //return the last opened subproblem
   Subproblem* sub = subs.back();
   subs.pop_back();
   return sub;
}
