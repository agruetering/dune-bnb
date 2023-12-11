// Copyright © Alexandra Grütering, see file LICENSE.md in module root
// License: GLP-3.0 or later
#include<dune/bnb/master.h>
#include<dune/bnb/subproblem.h>
#include<dune/bnb/branching.h>
#include<dune/bnb/enumeration.h>
#include<dune/bnb/depthfirst.h>

using namespace BranchAndBound;
using namespace std;

Master::Master(OptSense optSense_, double optEps_, int outFreq_) : optSense(optSense_), optEps(optEps_), outFreq(outFreq_), primalBoundFound(false), branchingStrategy(NULL), enumerationStrategy(NULL), primalBound(0), dualBound(0), subCount(0) {
   firstSub = new Subproblem(this);
}

Master::~Master(){
   if(enumerationStrategy != NULL)
      delete enumerationStrategy;
   if(branchingStrategy != NULL)
      delete branchingStrategy;
}

void Master::optimize(){

   clock_t cstart = clock();
   //if no branching strategy and no enumeration strategy was set
   //a default strategy is used.
   if(branchingStrategy== NULL)
      exit(1); 
   if(enumerationStrategy == NULL)
      enumerationStrategy = new DepthFirst();
   
    // add first subproblem to list of open problems
   subproblems.push_back(firstSub);
   Subproblem* activeSub;

   while(!subproblems.empty()){

	updateDualBound();

      	//choose new subproblem
      	activeSub = enumerationStrategy->chooseSubproblem(subproblems);

     	//get bounds for active subproblem 
 	bool status=getBounds(activeSub,subCount); 


      	//Output
      	subCount++;
      	if(!((subCount-1)%outFreq)) {
		cout<<"open subs: "<<::std::setw(6)<<::std::left<<subproblems.size();
		cout<<"processed subs: "<<::std::setw(10)<<::std::left<<subCount;
		cout<<"active SUB DB: "<<::std::setw(15);
		if(activeSub->isFeasible())
	  		cout<<activeSub->getDualBound();
		else
	  		cout<<"infeas";
		cout<<"PRIMAL BOUND: "<<::std::setw(15);
		if(primalBoundFound)
	  		cout<<primalBound<<endl;
		else
	  		cout<<"---"<<endl;
      	
	}

	 // check if subproblem can be fathomed 
	 if(fathomSubproblem(activeSub) && status){
		delete activeSub; 
		continue; // skip branching
	}

      	//Branch active subproblem 
      	list<Subproblem*> newSubs = branchingStrategy->branch(activeSub);  
      	delete activeSub;

      	if(!newSubs.empty()){
	 	//add new subproblems to list
		 for(list<Subproblem*>::iterator it = newSubs.begin(); it != newSubs.end(); it++)
	       		subproblems.push_back(*it);
	}
   }
   double time=(double) (clock()-cstart)/CLOCKS_PER_SEC; 
   cout.precision(30);
   cout<<"Number of subproblems: "<<subCount<<endl;
   cout <<"Time: " << time << endl;
   //Problem is solved
   cout<<"Optimum: "<<primalBound<<endl;
   for(int i = 0; i< primalSolution.size(); i++) {
      cout << primalSolution[i]<< endl;
   }

}


bool Master::fathomSubproblem(Subproblem* sub){
   if (!sub->isFeasible()){
      return true;
   }
   if(optSense == Minimize && primalBoundFound)
      if((primalBound-sub->getDualBound())/primalBound <= optEps)
	 return true;
   if(optSense == Maximize && primalBoundFound)
      if((sub->getDualBound()-primalBound)/primalBound <= optEps)
	 return true;
   return false;
}


void Master::updateDualBound(){
	double bestKnownValue = subproblems.front()->getDualBound();

   	if(optSense == Minimize){
      		for(list<Subproblem*>::iterator it = ++(subproblems.begin()); it != subproblems.end(); it++) {
	 		if((*it)->getDualBound() < bestKnownValue)
	    			bestKnownValue = (*it)->getDualBound();
      		}
   	}
	else if (optSense == Maximize){
      		for(list<Subproblem*>::iterator it = ++(subproblems.begin()); it != subproblems.end(); it++) {
	 		if((*it)->getDualBound() > bestKnownValue)
	    			bestKnownValue = (*it)->getDualBound();
      		}
   	}
  	dualBound = bestKnownValue;
}

