// Copyright © Alexandra Grütering, see file LICENSE.md in module root
// License: GLP-3.0 or later
#ifndef MASTER_H
#define MASTER_H

#include<list>
#include<vector>
#include<iostream>
#include<iomanip>
#include<stdlib.h>
#include<time.h>

namespace BranchAndBound
{
   class Subproblem;
   class BranchingStrategy;
   class EnumerationStrategy;

   enum VarType{
      Continuous,
      Integer
   };

   enum OptSense{
      Minimize,
      Maximize
   };

   class Master
   {
      protected:

	 OptSense optSense;
	 Subproblem* firstSub;

	 bool primalBoundFound;
	 double primalBound;
	 std::vector<double> primalSolution;
	 double dualBound;
         double optEps;
	 int outFreq;

	 int subCount;

	 std::list<Subproblem*> subproblems;

	 BranchingStrategy* branchingStrategy;
	 EnumerationStrategy* enumerationStrategy;

	 //Has to set dual and primal bounds and solution of "sub" if possible. 
	 //Otherwise set feasibility-flag
	 virtual int getBounds(Subproblem* sub,int) {return 0;};//, double& dualBound, vector<double>& dualSolution, bool& feasible);
	 /* return value: 
	  * false - pruning of subproblem not allowed 
	  * true - pruning of subproblem possible
	  */
	
	 void setPrimalBound(double& bound, std::vector<double> solution){
	    if(!primalBoundFound
	       || (optSense == Maximize && bound > primalBound)
	       || (optSense == Minimize && bound < primalBound)){
		  primalBound = bound;
		  primalBoundFound = true;
		  primalSolution = solution;
	    }
	 }

	 //Returns true if subproblem can be fathomed
	 bool fathomSubproblem(Subproblem* sub);

	 //Updates global dualBound
	 void updateDualBound();

      public:

	 Master(OptSense optSense_, double optEps_ = 1e-6,int outFreq_ = 1);
	 virtual ~Master();

	 void setBranchingStrategy(BranchingStrategy* branchingStrategy_){
	    branchingStrategy = branchingStrategy_;
	 }
	 void setEnumerationStrategy(EnumerationStrategy* enumerationStrategy_){
	    enumerationStrategy = enumerationStrategy_;
	 }
	 OptSense getOptSense(){
	    return optSense;
	 }

	 int getTotalSubs(){
	    return subCount;
	 }

	 double &getPrimalSolution(int i){
	    return primalSolution[i];
	 }

	 void optimize();

   };

}

#endif
