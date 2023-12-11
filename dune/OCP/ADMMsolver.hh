// Copyright © Alexandra Grütering, see file LICENSE.md in module root
// License: GLP-3.0 or later
#ifndef ADMM_SOLVER_HH
#define ADMM_SOLVER_HH

#include <dune/bnb/master.h>
#include <dune/OCP/ADMMsystem.hh>

using namespace BranchAndBound;

/* Solver for subproblems in branch-and-bound algorithm */ 
template<class MType>
class ADMMsolver
{
public: 
	typedef MType Master;
	typedef typename Master::GV GV;
	typedef typename Master::GFS GFS;
	typedef typename Master::RangeType RangeType; // range type discrete grid function
	typedef typename Master::vtype vtype; // type of Dof vector 

	typedef typename Master::ControlVector ControlVector; // type of control vector
	typedef typename ControlVector::field_type field_type; // field type of control vector
	typedef typename Master::Matrix Matrix; // matrix type

	typedef typename Master::AdjProblem AdjProblem;
	typedef typename Master::CGProblem CGProblem;

	using Vector = BlockVector<FieldVector<double,1>>; 

	typedef InverseOperator2Preconditioner<InverseP<Matrix,Vector,Vector>> Precond;

	enum{n=Master::n}; // number of switches
private: 

	const HeatProblemInterface<GV,n>& heat; //  problem class for Sf (with boundary data)
	AdjointProblemInterface<GV>& YdProblem; // problem class for yd (with boundary data)
	const vector<double>& dt; // temporal grid 
	const vector<double>& tgrid; // temporal grid of yd and ud
	const GFS&gfs; // spatial grid
	vector<pair<int,int>> I; // free control variables
	const vector<vtype>& Sf;
	const ControlVector& Gstar_eta; //  encodes  G*S*(Sf-yd) for target state
	const ControlVector& Gstar_q; // encodes G*S*SGu_fix for target state 
	const ControlVector& ud; // encodes ud for target state
	const ControlVector& ua; // lower bound for control variable
	const ControlVector& ub; // upper bound for control variable
	double alpha; // Thikonov term 

	double rho; // penalty term cutting plane
	double beta; // penalty term box constraints
	double eabs; // absolute tolerance for termination criteria 
	double erel; // relative tolerance for termination cirteria
	double tol; // error tolerance for termination criteria
	double reduction; // reduction factor for cg method
	double iter_lin; // maximum iterations of cg method
public: 
	ADMMsolver(const HeatProblemInterface<GV,n>& heat_,
 	AdjointProblemInterface<GV>& desired,
	const GFS& gfs_,const vector<double>& dt_, const vector<double>& tgrid_,
	double alpha_,
	double rho_, double beta_,
	double tabs, double trel, double tol_,double factor, double iter,
	const vector<pair<int,int>>& index, // fixed indices of sub problem
	const vector<vtype>& Sf_,
	const ControlVector& Gstar_eta_,
	const ControlVector& Gstar_q_,  
	const ControlVector& ud_, 
	const ControlVector& ua_,
	const ControlVector& ub_): heat(heat_), YdProblem(desired), gfs(gfs_), dt(dt_),tgrid(tgrid_), alpha(alpha_), rho(rho_), beta(beta_), eabs(tabs), erel(trel), tol(tol_), reduction(factor),iter_lin(iter), Sf(Sf_), Gstar_eta(Gstar_eta_),Gstar_q(Gstar_q_), ud(ud_), ua(ua_), ub(ub_)
	{	

		for(int i=0; i < dt.size(); i++){ 
			for(int j=0; j < n; j++){
				if (find(index.begin(), index.end(), make_pair(j,i)) == index.end()) // (j,i) not fixed
					I.push_back(make_pair(j,i));
			}
		}
	}

	// ADMM Algorithm root node relaxation (no switching constraints)
	void apply(ControlVector& u, Vector&w, Vector&nu, double&opt, int& iteration); 
	
	// ADMM Algorithm root node relaxation (if switching constraints existent)
	bool apply(ControlVector& u, Vector& v, Vector& p,Vector& x, Vector &nu, const Matrix &A, const Vector& b,double& opt, int& iterations, double time); 

	//  ADMM Algorithm for B&B (no switching constraints)
	bool apply(ControlVector& u,Vector& w, Vector& nu,double& opt,const OptSense& optSense,const bool& primalBoundFound, const double& primalBound,const double& optEps, int &iteration);
	
	// ADMM algorithm for B&B (if switching constraints existent)
	bool apply(ControlVector& u, Vector& v, Vector& p,Vector& x, Vector &nu, const Matrix &A, const Vector& b,double& opt,const OptSense& optSense,const bool& primalBoundFound, const double& primalBound,const double& optEps, int &iterations);  

	// primal objective value (if switching constraints existent) 
	void primalObj(const ControlVector& u, double& opt);  

};

#include <src/OCP/ADMMsolver.cc>

#endif //ifndef ADMM_SOLVER_HH
