// Copyright © Alexandra Grütering, see file LICENSE.md in module root
// License: GLP-3.0 or later
#ifndef OCPMASTER_HH
#define OCPMASTER_HH

#include <type_traits>
#include <sys/stat.h>

#include <dune/grid/io/file/gmshwriter.hh>
#include <dune/istl/matrixmarket.hh>

#include <dune/bnb/master.h>
#include <dune/bnb/OCPSub.hh>
#include <dune/OCP/hfunctions.hh>
#include <dune/OCP/heatdriver.hh>
#include <dune/OCP/adjointdriver.hh>
#include <dune/OCP/switchconstr.hh> 
#include <dune/bnb/errorfunctions.hh>
#include <dune/OCP/ADMMsolver.hh>

using namespace Dune; 
using namespace BranchAndBound;
using namespace std;


/* parameter structure of branch and bound algorithm for
 *
 *    min J(y,u)=1/2|y-y_d|_{L^2(Q)}^2 + \apha/2 |u-u_d|_{L^2(0,T)}^2
 *    s.t.  \partial_t y - \Delta y = Gu + f = \sum_j u_j(t)*Psi_j(x) + f	in Q=\Omega x (0,T)
 *                       		   y = g   				on \Sigma_D=\Gamma_D x (0,T)
 *         		    \nabla y \cdot n = j   				on \Sigma_N=\Gamma_N x (0,T)
 *                        		y(0) = y0 		
 * 				           u in D 
 *
 *
 * subproblems solved by ADMM method
*/
template<class GridView, class Constraints, class Problem, int k, int number>
class OCPMaster : public Master 
{
public: 
	enum{n=number}; 
	enum {dim=GridView::dimension}; // dimension of spatial grid
	enum {degree = k}; // degree of spatial ansatz and test functions
	typedef GridView GV; // spatial grid view 
	typedef Constraints CON; // constraints (e.g. conforming dirichlet constraints)

	typedef OCPMaster<GV,CON,Problem,k,number> MType;
	typedef OCPsub<MType> SubType;

	typedef typename GV::ctype DF; // coordinate type of the spatial grid
	
	typedef typename PDELab::PkLocalFiniteElementMap<GV,DF,double,degree> FEM; // finite element map (for simplices in case dim>=2)
  	typedef PDELab::ISTL::VectorBackend<> VBE; // ISTL backend to storage Dofs
  	typedef typename PDELab::GridFunctionSpace<GV,FEM,CON,VBE> GFS; // grid function space for spatial ansatz and test functions (Galerkin method assumed)

	typedef typename PDELab::Backend::Vector<GFS,double> vtype; // type of Dof vectors
	typedef typename PDELab::DiscreteGridFunction<GFS,vtype> DGF; // discrete grid function
	typedef typename DGF::Traits::RangeType RangeType; // range type discrete grid function

     	typedef PDELab::LocalFunctionSpace<GFS> LFS; // local function space 
	typedef PDELab::LFSIndexCache<LFS> LFSCache; // local function space cache
	typedef typename GFS::Traits::FiniteElementType::Traits::LocalBasisType LocalBasis; // local basis 
	typedef typename PDELab::LocalBasisCache<LocalBasis> BasisCache; // local basis cache

	typedef BlockVector<FieldVector<double,n>> ControlVector; // type of control vector
	typedef BCRSMatrix<FieldMatrix<double,1,1>>Matrix; // matrix type for switching constraints

	typedef Problem CGProblem; // problem interface for SGv
	typedef AdjointProblem<GFS> AdjProblem; // problem interface for S*SGv 

	typedef SwitchPoly BaseSwitchSet; 

	using Vector = BlockVector<FieldVector<double,1>>; 

protected:	
	HeatProblemInterface<GridView,n>& heat; //  problem class for Sf (with boundary data)
	AdjointProblemInterface<GridView>& YdProblem; // problem class for yd (with boundary data)
	vector<double> tgrid; // temporal grid of yd and ud
	double alpha; // Thikonov term
	double rho; // penalty term cutting planes
	double beta; // penalty term box constraints

	int iter_lin; // maximal number of iterations for linear system 
	double reduction; // error reduction factor
	double eabs; // absolute tolerance for ADMM method
	double erel; // relative tolerance for ADMM Method
	double tol; // error tolerance for ADMM */

	double rfactor; // refinement factor
	double rchange; // minimal relative change of objective value

	int iter_outer; // maximal number of outer approximation iterations for each subproblem
	double timeouter; // time limit for outer approximation of each subproblem

	vector<double> primalTimeGrid; // time grid of primal solution

	const BaseSwitchSet& D; // switching polytope

public: 

	/* constructor and destructor */
	OCPMaster(GV gv, 
	HeatProblemInterface<GridView,n>& heat_,
 	AdjointProblemInterface<GridView>& desired,
	vector<double> grid, 
	const BaseSwitchSet& set, 
	double alpha_, 
	double rho_, 
	double beta_,
	vector<double> dt, 
	int iter, double factor,
	double tabs,
	double trel,
	double tol_,
	double rfac=0.5,
	double relchange=0.01,
	double optEps=1e-2,
	int iter_o=100, 
	double limit=1800, 
	int outFreq=1) 
	: heat(heat_), YdProblem(desired), tgrid(grid),D(set), 
	alpha(alpha_), rho(rho_), beta(beta_),
	iter_lin(iter), reduction(factor), eabs(tabs), erel(trel), tol(tol_),
	rfactor(rfac), rchange(relchange), iter_outer(iter_o), timeouter(limit), 	
	Master(Minimize,optEps,outFreq) 
	{ 
		// default branching strategy
		branchingStrategy=new OCPBranch<SubType>(); 
		// create first subproblem
	 	delete firstSub;
		firstSub=new SubType(this, gv, dt);
		// set up desired temperature and inhomogeneuous PDE data for first subproblem
		initialize(firstSub);
	}

	~OCPMaster() {};


	/* master class construction  with input data */
	OCPMaster(GV gv, 
	HeatProblemInterface<GridView,n>& heat_,
 	AdjointProblemInterface<GridView>& desired, 
	vector<double> grid, 
	const BaseSwitchSet& set, 
	vector<pair<int,double>> tau, vector<double> d, vector<pair<int,int>> index,
	ControlVector u, 
	Matrix A, Vector b, Vector v, Vector mu,
	Vector w, Vector nu,
	double alpha_,
	double rho_,
	double beta_, 
	vector<double> dt,
	int iter, double factor,
	double tabs,
	double trel,
	double tol_,
	double rfac,
	double relchange,
	double optEps,
	int iter_o,
	double limit, int outFreq, 
	int count) 
	: heat(heat_), YdProblem(desired), tgrid(grid), D(set), 
	alpha(alpha_), rho(rho_), beta(beta_),  
	iter_lin(iter), reduction(factor), eabs(tabs), erel(trel), tol(tol_),
	rfactor(rfac),rchange(relchange), iter_outer(iter_o), timeouter(limit), 	
	Master(Minimize,optEps,outFreq)
	{ 
		this->subCount=count;
		// default branching strategy
		branchingStrategy=new OCPBranch<SubType>(); 
		// create first subproblem
	 	delete firstSub;
		firstSub=new SubType(this, gv, dt, tau, d, index, u, A, b,v, mu, w, nu);
		GFS gfs=((SubType*)firstSub)->getGFS();
		// set up desired temperature and inhomogeneuous PDE data for first subproblem
		initialize(firstSub);
		updateFixData(firstSub,index,gfs,dt,u);
		((SubType*) firstSub)->GridChanged();
	}

	double Rho(){return rho;}
	
	/* set up desired temperature and inhomogeneuous PDE data for root node*/
	void initialize(Subproblem* sub);

	void updateFixData(Subproblem* sub,const vector<pair<int,int>>& index, const GFS& gfs, const vector<double>& dt,const ControlVector& v);
	/* update subproblem data due to changed fixed control vector 
	 * parameters:
 	 * - given subproblem
	 * - indices of fixed control variables
	 * - spatial grid function space
	 * - temporal grid
	 * - control v
	 */

	void updateData(vector<pair<int,int>>&,const GFS&,const vector<double>&, vector<double>&,vector<vtype>&,vector<vtype>&, ControlVector&, ControlVector&,Matrix&, Vector&, Vector&, ControlVector&);
	/* method to update data after grid refinement 
	 * parameters:
	 * - indices of fixed control variables
	 * - spatial grid function space
	 * - temporal grid
	 * - old temporal grid 
	 * - Sf of subproblem 
 	 * - S*yd of subproblem
	 * - G*S*(Sf-yd) 
	 * - G*S*SGu_fix
         * - matrix of cutting planes
	 * - auxiliary vector for box constraints
	 * - lagrange multiplicator of box constraints 
	 * - current control v
	 */

	 void PrimalBound(const GFS&,const vector<double>&,const vector<vtype>&,const vector<vtype>&, const ControlVector&); 
	/* primal bound calculation 
	 * parameters:
	 * - spatial grid function space
	 * - temporal grid of v 
	 * - Sf of subproblem 
 	 * - S*yd of subproblem
	 * - control v
	 */

	double objective(const GFS&, const vector<double>&,const vector<vtype>&, const ControlVector &); 
	/* exacter objective value calculation 
	 * parameters:
	 * - spatial grid function space
	 * - temporal grid of v
	 * - Sf of subproblem 
	 * - control v 
	 */


	double J(const GFS&,const vector<double>&,const ControlVector&,const vector<vtype>&);
	/* calculation of J(y,u)
         * parameters:
	 * - spatial grid function space
	 * - temporal grid 
	 * - control u 
	 * - y=S(Gu+f)
         * return value: objective value
	 */
	
	void call_outer(Subproblem*,const GFS&,const vector<double>&,const vector<pair<int,int>>&,const vector<vtype>&, const ControlVector&,const ControlVector&,Matrix&, Vector&, Vector&, Vector&, Vector&, Vector&, double&,ControlVector&,const int, OutputData&);
	/* outer approximation of subproblem 
	 * parameters:
 	 * - given subproblem with fixations
	 * - spatial grid function space
	 * - temporal grid 
	 * - indices of fixed control variables
	 * - Sf
	 * - G*S*(Sf-yd) 
	 * - G*S*SGu_fix
         * - matrix of cutting planes
	 * - right-hand side of cutting planes
	 * - auxiliary vector for cutting planes
         * - lagrange muliplicator of cutting planes
	 * - auxiliary vector for box constraints
	 * - lagrange multiplicator of box constraints 
	 * - objective value 
	 * - current control v
	 * - subproblem count
         * - structure for output data
	 */

	/* bounds */ 
	int getBounds(Subproblem*,int);
	void setPrimalBound(double& bound, vector<double> solution, vector<double> dt){
	    if(!primalBoundFound
	       	|| (optSense == Maximize && bound > primalBound)
	       	|| (optSense == Minimize && bound < primalBound)){
		  	primalBound = bound;
		  	primalBoundFound = true;
		 	primalSolution = solution;
			primalTimeGrid=dt;
			
	     }
	 }
	
	/* write methods */
	void writePrimalInfo(int);
	void writeParameterSetting();
	
	virtual int separate(Subproblem*,const vector<double>&,const ControlVector&, double*, double&)=0;
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
  
	virtual void heuristic(Subproblem*,const GFS&,const vector<double>&,const ControlVector&, ControlVector&)=0;
	/* compute a heuristic solution
	 * parameters:
	 * - given subproblem with fixations
	 * - spatial grid function space
	 * - temporal grid 
	 * - relaxed control
	 * - heuristic control
	 */

	virtual int updateFixControls(Subproblem* , const vector<double>&,vector<pair<int,int>>&, ControlVector&  )=0; 
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
	
	int refine(Subproblem*,GFS&,vector<double>&,vector<double>&);
	/* refinement strategy of grids 
	 * parameters:
         * - given subproblem
	 * - spatial grid function space 
	 * - temporal grid 
	 * - cellwise error vector 
	 * return value:
	 * -1: grid not refined since minimal grid wide reached
	 *  1: grid refined
	 */

	void Error(const GFS&,const vector<double>&, const vector<vtype>&,const vector<vtype>&, const vector<pair<int,int>>& , const Matrix&, const Vector&, const ControlVector&, vector<double>&);
	/* cellwise error calculation in optimal control
	 * parameter: 
	 * - spatial grid function space 
	 * - temporal grid
	 * - Sf (of subproblem) 
	 * - S*yd (of subproblem) 
	 * - indices of fixed control variables
	 * - matrix of cutting planes
         * - lagrange muliplicator of cutting planes
	 * - control u 
 	 * - cellwise error vector 
	 */

};

#include <src/bnb/OCPMaster.cc>

#endif


