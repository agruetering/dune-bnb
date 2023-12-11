// Copyright © Alexandra Grütering, see file LICENSE.md in module root
// License: GLP-3.0 or later
#ifndef OCPDMAX_HH
#define OCPDMAX_HH

#include <dune/bnb/OCPMaster.hh>
#include <dune/OCP/Dmax.hh>

using namespace BranchAndBound; 
using namespace std; 
using namespace Dune; 

/* master problem for
 *    min J(y,u)=1/2|y-y_d|_{L^2(Q)}^2 + \apha/2 |u-u_d|_{L^2(0,T)}^2
 *    s.t.  \partial_t y - \Delta y = Gu + f = u(t)*Psi(x) + f	in Q=\Omega x (0,T)
 *                       		   y = g   		on \Sigma_D=\Gamma_D x (0,T)
 *         		    \nabla y \cdot n = j   		on \Sigma_N=\Gamma_N x (0,T)
 *                        		y(0) = y0 		
 * 				           u in D 
 *
 * with D={u: |Du|([0,T])<= \sigma}
 * 
 * and subproblems solved by ADMM
 */

template<class GridView, class Constraints, class Problem, int k>
class OCPDmax : public OCPMaster<GridView, Constraints, Problem, k,1>
{
public: 
	typedef OCPMaster<GridView, Constraints, Problem, k,1> BaseClass; 
	typedef Dmax BaseSwitchSet; 
	using SubType=typename BaseClass::SubType;
	using GV=typename BaseClass::GV;
	using GFS=typename BaseClass::GFS;
	using vtype=typename BaseClass::vtype;
	using MType=typename BaseClass::MType;
	using CGProblem=typename BaseClass::CGProblem;
	using Matrix=typename BaseClass::Matrix; 
	using ControlVector=typename BaseClass::ControlVector; 

	using Vector = BlockVector<FieldVector<double,1>>;
	enum{n=BaseClass::n};

protected:
	int sigma;
	using BaseClass::rfactor;
	using BaseClass::heat;
	using BaseClass::tgrid;  
	using BaseClass::alpha; 
	using BaseClass::beta;
	using BaseClass::rho; 
	using BaseClass::eabs;
	using BaseClass::erel;
	using BaseClass::tol;
	using BaseClass::reduction; 
	using BaseClass::iter_lin; 
	using BaseClass::iter_outer; 
	using BaseClass::timeouter; 
	
	using BaseClass::primalBound;
	using BaseClass::primalBoundFound; 
	using BaseClass::optEps;
 	using BaseClass::optSense;
	

public:
	/* constructor and destructor of BaseClass*/
	OCPDmax(
	GV gv, 
	HeatProblemInterface<GridView,1>& heat,
 	AdjointProblemInterface<GridView>& desired,
	vector<double> grid, 
	const Dmax& set,
	double alpha,
	double rho,
	double beta,
	vector<double> dt, 
	int iter, double factor,
	double tabs,
	double trel,
	double tol,
	double rfac=0.5,
	double rchange=0.01,
	double optEps=1e-2, 
	int iter_o=100, 
	double limit=1800,  
	int outFreq=1) 
	: sigma(set.getSigma()), 
	BaseClass(gv,heat,desired,grid,set,alpha,rho,beta,dt,iter,factor,tabs,trel,tol,rfac,rchange,optEps,iter_o,limit,outFreq)
	{}

	~OCPDmax() {};

	/* master class construction  with input data */
	OCPDmax(
	GV gv, 
	HeatProblemInterface<GridView,1>& heat,
 	AdjointProblemInterface<GridView>& desired,
	vector<double> grid, 
	const BaseSwitchSet& set, 
	vector<pair<int,double>> tau, vector<double> c, vector<pair<int,int>> index,
	ControlVector u,
	Matrix A, Vector b,Vector v, Vector mu,
	Vector w, Vector nu,
	double alpha, 
	double rho,
	double beta,
	vector<double> dt,
	int iter, double factor,
	double tabs,
	double trel,
	double tol,
	double rfac,
	double rchange,
	double optEps, 
	int iter_o, 
	double limit, int outFreq, 
	int count) 
	: sigma(set.getSigma()),
	BaseClass(gv,heat,desired,grid,set,tau, c,index, u, A, b, v, mu, w, nu, alpha,rho, beta, dt,iter,factor,tabs,trel,tol,rfac,rchange,optEps,iter_o,limit,outFreq,count)
	{}

	Subproblem* getfirstSub() {return this->firstSub;}

	void writeParameterSetting()
	{
		stringstream filename;
		filename << BNB_OUTPUT_PATH << OUTPUTDIR << "/settings.txt";
		fstream stream; 
		stream.precision(10);
		stream.open(filename.str(), ios::out | ios::app);
		stream << "smax=" << sigma << endl;
		stream.close();
		BaseClass::writeParameterSetting(); 
	}

	int separate(Subproblem* sub,const vector<double>& dt,const ControlVector& u, double* cut, double& rhs);
	/* compute a violated constraint, if any 
	 * parameters:
	 * - given subproblem with fixations
	 * - temporal grid 
	 * - current switching pattern
	 * - coefficients of violated constraint, if any
	 * - right-hand side of violated constraint, if any
  	 * return value:
  	 * - 0: no violated constraint found
  	 * - 1: violated inequality found
 	 */

	void heuristic(Subproblem* sub,const GFS& gfs,const vector<double>& dt, const ControlVector& u, ControlVector& v);
	/* compute a heuristic solution
	 * parameters:
	 * - given subproblem with fixations
	 * - spatial grid function space
	 * - temporal grid 
	 * - relaxed control
	 * - heuristic control
	 */

	int updateFixControls(Subproblem* sub,const vector<double>& dt, vector<pair<int,int>>& fixIndices, ControlVector& u);
	/* update fixations of control variables 
	 * parameters: 
	 * - given subproblem 
	 * - temporal grid
	 * - vector with fixed indices
         * - control u_old
	 * 
	 * after adding fixed indices
	 *	- change values of u_old for fixed indices and call setControl of subproblem 
	 *	- call updateData in master problem to calculate G*S*Su_old
	 * return value: 
	 * -0: no fixed indices added
	 * -1: fixed indices added
	 */


};

# include <src/bnb/OCPDmax.cc>

#endif
