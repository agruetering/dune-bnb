// Copyright © Alexandra Grütering, see file LICENSE.md in module root
// License: GLP-3.0 or later
#ifndef SUBPROBLEM_H
#define SUBPROBLEM_H

#include<dune/bnb/master.h>
#include<vector>
#include<utility>

namespace BranchAndBound
{

   /* Branch and Bound subproblem class */
   class Subproblem
   {
	 Master* pMaster; // master problem 

      private:
	 
	 std::vector<std::pair<int,double>> tau; // fixation points (chronologically ordered by first switches and second time) 
	 std::vector<double> d; // fixation values
	 int numPoints; // number of fixations
	
         double dualBound;
	 std::vector<double> dualSolution;
	 bool feasible;

	bool oaflag; // flag if outer approximation stopped with feasible control
	

      public:

        Subproblem(Master* pMaster_,std::vector<std::pair<int,double>> tau_=std::vector<std::pair<int,double>>(), std::vector<double> d_=std::vector<double>()) 
	  : pMaster(pMaster_), tau(tau_), d(d_), feasible(true), dualBound(0), oaflag(0),numPoints(tau.size()) 
	 {};

        Subproblem(Subproblem* pParent_,std::vector<std::pair<int,double>> tau_, std::vector<double> d_) 
	  : pMaster(pParent_->pMaster), tau(tau_), d(d_), feasible(true), dualBound(pParent_->dualBound), oaflag(0),numPoints(tau.size()) 
	 {};

	 virtual ~Subproblem() {};
	 
	 //Setter
	 inline void setDualBound(double dualBound_) { dualBound = dualBound_; }
	 inline void setDualSolution(std::vector<double> dualSolution_) {
	    dualSolution = dualSolution_;
	 }
	 inline void setFeasible(bool feasible_) { feasible = feasible_; }
	 // Set that outer approximation stopped with feasible control
	 inline void OAstopped(bool stop){ oaflag=stop; }

	 //Getter
	 inline Master* master() { return pMaster; }
	 inline double getDualBound() { return dualBound; }
	 inline bool isFeasible() { return feasible; }
	 inline bool doesOAstopped() { return oaflag; }
	 inline std::vector<double> getDualSolution() { return dualSolution; }
	 inline std::vector<std::pair<int,double>> &getTau() { return tau; }
	 inline std::vector<double> &getFix() { return d; }
	 inline int getNumFixPoints() { return numPoints; }
   };
}

#endif
