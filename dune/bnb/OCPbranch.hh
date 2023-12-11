#ifndef OCPBRANCH_HH
#define OCPBRANCH_HH

#include<dune/bnb/subproblem.h>
#include<dune/bnb/branching.h>
#include<vector>
#include<list>

namespace BranchAndBound {
  
   // chooses the variable with highest deviation to 0/1 (minimal distance to 0/1 multiplied by interval length) 
   template<class SubType>
   class OCPBranch : public BranchingStrategy {
      public:
	 std::list<Subproblem*> branch(Subproblem* sub)	
	 {
		const int n=1; // number of switches
		assert(SubType::n==1 && "branching strategy only for a single switch");

		using OCPMaster=typename SubType::Master;
		typedef typename SubType::GFS GFS;
		typedef typename SubType::ControlVector ControlVector; 
		typedef typename SubType::Matrix Matrix; 
		typedef typename SubType::Vector Vector; 

		vector<double> dualSolution = sub->getDualSolution();
   		vector<pair<int,double>> tau= sub->getTau();
   		vector<double> c=sub->getFix();
   		vector<double> dt=((SubType*)sub)->getdeltaT();
   
   		//initialize parameters
   		int branchCandidate = -1; // which switch 
		int branchPoint = -1; // index of time point
  		double topdist = 0;
 		double tmpdist;


   		// find variable with highest deviation from 0/1 (minimal distance to 0/1 multiplied by interval length) 
   		for(int j =0 ; j < n; j++){
			for(int i=0; i < dt.size(); i++){
				tmpdist =min(abs(dualSolution[i*n+j]),abs(1-dualSolution[i*n+j]));
				tmpdist*=dt[i];
	 			if(tmpdist > topdist){
					// prove that time points is not fixed yet
	    				double t=accumulate(dt.begin(),dt.begin()+i+1,0.0); // candidate for time point

					// current fixationsn for j_th control
					vector<double> h;
					for(auto& entry :tau){
						if(entry.first==j)
							h.push_back(entry.second);
					}
	    				for(double& element : h) 
						element-=t; 
	    				bool fixed=false; 
            				for(int i=0; i < h.size(); i++){
						if(abs(h[i])<=1e-10){	// time point already fixed
							fixed =true; 
							break;
						}
					}
					if(!fixed){
	    					branchCandidate = j;
						branchPoint = i;
	    					topdist = tmpdist;
					}
				}
   			}
		}

   		// create new subproblems
		typedef typename SubType::GFS GFS;
		typedef typename SubType::vtype vtype;
   		list<Subproblem*> newSubs;
   		if(branchCandidate != -1){
			double t=accumulate(dt.begin(),dt.begin()+branchPoint+1,0.0); // new fixation point
			int i=0;
			while(i< tau.size() && (tau[i].first<branchCandidate || tau[i].second < t))
					i++;
 	
	
			// insert new fixation point at right position
			vector<pair<int,double>> ntau=tau; 
			ntau.emplace(ntau.begin()+i,pair(branchCandidate,t));
			// insert fixation values at right position
			vector<double> c1=c; 
			vector<double> c2=c; 
			c1.emplace(c1.begin()+i,0); 
			c2.emplace(c2.begin()+i,1);
			GFS gfs=((SubType*)sub)->getGFS();
			Subproblem *s1 = new SubType(sub,gfs,ntau,c1);
   			newSubs.push_back(s1);
   			Subproblem *s2 =  new SubType(sub,gfs,ntau,c2);
   			newSubs.push_back(s2);		
   		}
   		else {
			/* all possible branching points fixed */
				// get subproblem data
				ControlVector u=((SubType*)sub)->getControl();
				GFS gfs=((SubType*)sub)->getGFS();
				ControlVector Gstar_q,Gstar_eta;
				vector<vtype> Sf;
				Matrix A; 
				Vector v, mu,w, nu, b;
				vector<pair<int,int>> index;
				((SubType*)sub)->getLambda(index,Sf,Gstar_eta,Gstar_q, A, b,v, mu, w, nu);
				
				// refine problem
				vector<double> dt_old=dt;
				Master* pmaster=sub->master();
				vector<vtype> g=((SubType*)sub)->getStaryd();
				vector<double> e(dt.size());
				BlockVector<FieldVector<double,1>> mu_=mu; 
				mu_*=((OCPMaster*)pmaster)->Rho();
				((OCPMaster*)pmaster)->Error(gfs,dt,Sf,g,index,A,mu_,u,e);
				int refined=((OCPMaster*)pmaster)->refine(sub,gfs,dt,e);
				if(refined==1){
					((OCPMaster*)pmaster)->updateData(index,gfs,dt,dt_old,Sf,g,Gstar_eta,Gstar_q,A,w,nu,u); 

					// update data of subproblem
					((SubType*)sub)->setData(index, gfs, dt, Sf,g, Gstar_eta,Gstar_q,A, w, nu, u);
				
					Subproblem *s1 = new SubType(sub,gfs,tau,c);
					((SubType*)s1)->GridChanged();
					newSubs.push_back(s1);
				}
   		//	}
		}
   		return newSubs;
	}
   };
}
#endif
