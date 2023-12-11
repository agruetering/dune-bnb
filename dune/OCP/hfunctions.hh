// Copyright © Alexandra Grütering, see file LICENSE.md in module root
// License: GLP-3.0 or later
#ifndef HELPFUNCTIONS_HH
#define HELPFUNCTIONS_HH

using namespace Dune; 
using namespace std;

#include<dune/OCP/probleminterface.hh>


/************************************************************************** 
 *
 * contribution G*v =(\int Psi_j(x)v(t,x) dx )_{1<=j<=n} 
 *
***************************************************************************/
template<class GFS, class Problem,  class ControlVector>
void Gstar(const GFS& gfs, // grid functions space 
const Problem& problem, // problem class containing form functions
const vector<typename PDELab::Backend::Vector<GFS,double>>& v, // Dof vector v(t) at every time point
ControlVector& y) // vector to store result
{
    typedef typename GFS::Traits::GridView GV; // spatial grid view
    typedef typename GFS::Traits::FiniteElementMap FEM; // finite element map
    typedef typename PDELab::Backend::Vector<GFS,double> Y; // type of Dof vector
    typedef typename PDELab::DiscreteGridFunction<GFS,Y> DGF; // discrete grid function

  
    typedef typename DGF::Traits::RangeType RangeType; // range type discrete grid function
    typedef PDELab::LocalFunctionSpace<GFS> LFS; // local function space 
    typedef PDELab::LFSIndexCache<LFS> LFSCache; // local function space  cache
    typedef typename GFS::Traits::FiniteElementType::Traits::LocalBasisType LocalBasis; // local basis 
    typedef typename PDELab::LocalBasisCache<LocalBasis> BasisCache; // local basis cache
	
    // local function space and basis
    LFS lfs(gfs);
    LFSCache lfs_cache(lfs);
    BasisCache cache;

    // clear y
    y.resize(v.size()-1);
    y=0;

    // loop over all elements of leaf grid
    for (const auto& eg : elements(gfs.gridView())){

	// local function space on current element
	lfs.bind(eg); 
       	lfs_cache.update();
        auto geo = eg.geometry();
	auto ref = referenceElement(geo);

        // choose appropiate quadrature rule
        const int order = -1+2*gfs.finiteElementMap().maxLocalSize();
        auto rule = PDELab::quadratureRule(geo,order);

	// loop over quadrature points
	for (const auto& ip : rule)
	{
		decltype(ip.weight()) factor = ip.weight()* geo.integrationElement(ip.position());

		// evaluate basis functions
		auto phihat = cache.evaluateFunction(ip.position(),
               lfs.finiteElement().localBasis());

            	// evaluate Psi
        	auto psi=problem.Psi(eg,ip.position());

            	for (int i=0; i < v.size()-1; i++){
                	RangeType val(0.0); 
                	// evaluate v(t_i,x)
                	for (size_t j=0; j<phihat.size(); j++) {
				val+=v[i][lfs_cache.containerIndex(j)]*phihat[j]; 
			} 
			val*=factor;
                	y[i].axpy(val,psi);
            	}
        }
   }
}



/**************************************************************************
 *
 * calculation of 1/2 ||y-y_d||^2_{L_2(Q)} 
 *
***************************************************************************/
template<class GFS>
void L2Deviation(const GFS& gfs, // spatial grid function space
const vector<double>& dt, // temporal grid 
const vector<typename PDELab::Backend::Vector<GFS,double>>& y, // Dof vector 
AdjointProblemInterface<typename GFS::Traits::GridView>& YdProblem, // problem class for yd
double& bound)
{
	typedef typename GFS::Traits::GridView GV; // spatial grid view
    	typedef typename GFS::Traits::FiniteElementMap FEM; // finite element map
    	typedef typename PDELab::Backend::Vector<GFS,double> Y; // type of Dof vector
   	 typedef typename PDELab::DiscreteGridFunction<GFS,Y> DGF;  // discrete grid function

    	typedef typename DGF::Traits::RangeType RangeType; // range type discrete grid function
    	typedef PDELab::LocalFunctionSpace<GFS> LFS; // local function space 
    	typedef PDELab::LFSIndexCache<LFS> LFSCache; // local function space cache
    	typedef typename GFS::Traits::FiniteElementType::Traits::LocalBasisType LocalBasis; // local basis 
    	typedef typename PDELab::LocalBasisCache<LocalBasis> BasisCache; // local basis cache

	YdProblem.setdeltaT(dt);

    	LFS lfs(gfs);
    	LFSCache lfs_cache(lfs);
    	BasisCache cache;
	
    	bound=0; 
    	// 1/2 int_\Omega \int_[0,T] (y-y_d)^2 dt dx
    	// loop over all elements of leaf grid
    	for (const auto& eg : elements(gfs.gridView())){
		
        	lfs.bind(eg); 
        	lfs_cache.update();
       	 	auto geo = eg.geometry();

        	// choose appropiate quadrature rule
       		const int order =-1+2*gfs.finiteElementMap().maxLocalSize();
        	auto rule = PDELab::quadratureRule(geo,order);

        	// loop over quadrature points
    		for (const auto& ip : rule)
     		{
            		decltype(ip.weight()) factor = ip.weight()* geo.integrationElement(ip.position());

            		// evaluate basis functions
            		auto phihat = cache.evaluateFunction(ip.position(),
                             lfs.finiteElement().localBasis());

                    	// loop over all time points
           		for (int i=0; i < dt.size(); i++){
                        	// evaluate desired temperature
				double t=accumulate(dt.begin(),dt.begin()+i+1,0.0);
				double des=YdProblem.rhs(eg,ip.position(),t);
               	 		RangeType val;
                		// evaluate y(t_i,x)
                		for (size_t j=0; j<lfs.size(); j++){
				        val+=y[i+1][lfs_cache.containerIndex(j)]*phihat[j];
                		}
				
                		double term; 
				term=pow(val-des,2);
                		term*=factor; 
				term*=dt[i];
                        	bound+=term; 
            		}
        	}
	}  
	bound*=1.0/2.0;

}


/***************************************************************************
 *
 * extend control
 *
 ***************************************************************************/
template<class ControlVector>
void extendVector(const ControlVector& u, // control
const vector<double>& dt, // temporal grid of control
const vector<double>& s, // temporal grid of yd
ControlVector& v, // extended vector 
vector<double>& grid) // extended grid
{
	/* time points of temporal grid */
	vector<double> t(dt.size());
	t[0]=dt[0]; 
	for(int i=1; i < dt.size(); i++)
		t[i]=t[i-1]+dt[i];

	/* different time points of both temporal grids */
	vector<double> t_new=t; 
	for(int i=0; i < s.size(); i++){
		double time=accumulate(s.begin(),s.begin()+i+1,0.0);
		if(find_if(t_new.begin(),t_new.end(),[time](const double& a){return abs(a-time)<=1e-10;}) == t_new.end())
			t_new.push_back(time);
	}
	sort(t_new.begin(),t_new.end());	

	v.resize(t_new.size());
	for(int i=0; i < t_new.size(); i++){
		// first time point of temporal grid of control that is not less than t_new[i]
		auto start=lower_bound(t.begin(),t.end(), t_new[i],[](const double& a, double value){return a < value-1e-10;});
		int index=distance(t.begin(),start);
		// set value of new control vector on (t_new[i-1],t_new[i]]
		v[i]=u[index]; 	
	}

	/* temporal grid of all different time points */ 
	grid.resize(t_new.size()); 	
	grid[0]=t_new[0];
	for(int i=1; i<t_new.size(); i++)
		grid[i]=t_new[i]-t_new[i-1]; 
}

/***************************************************************************
 *
 * extend primal solution 
 *
 ***************************************************************************/
template<class ControlVector>
void extendVector(const vector<double>& sol, // primal solution
const vector<double>& dt, // temporal grid of primal solution
const vector<double>& s, // temporal grid of yd
ControlVector& u, // extended vector 
vector<double>& grid) // extended grid
{
	const int n=ControlVector::block_type::dimension;
	
	/* time points of temporal grid */
	vector<double> t(dt.size());
	t[0]=dt[0]; 
	for(int i=1; i < dt.size(); i++)
		t[i]=t[i-1]+dt[i];

	/* different time points of both temporal grids */
	vector<double> t_new=t; 
	for(int i=0; i < s.size(); i++){
		double time=accumulate(s.begin(),s.begin()+i+1,0.0);
		if(find_if(t_new.begin(),t_new.end(),[time](const double& a){return abs(a-time)<=1e-10;}) == t_new.end())
			t_new.push_back(time);
	}
	sort(t_new.begin(),t_new.end());

	u.resize(t_new.size());
	for(int i=0; i < t_new.size(); i++){
		// first time point of temporal grid that is not less than t_new[i]
		auto start=lower_bound(t.begin(),t.end(), t_new[i],[](const double& a, double value){return a < value-1e-10;});
		int index=distance(t.begin(),start);
		// set value of new control vector on (t_new[i-1],t_new[i]]
		for(int j=0; j<n; j++)
			u[i][j]=sol[index*n+j]; 	
	}

	/* temporal grid of all different time points */ 
	grid.resize(t_new.size()); 	
	grid[0]=t_new[0];
	for(int i=1; i<t_new.size(); i++)
		grid[i]=t_new[i]-t_new[i-1]; 
}

/***************************************************************************
 *
 * extend state
 *
 ***************************************************************************/
template<class GFS>
void extendVector(
const vector<typename PDELab::Backend::Vector<GFS,double>>& y, // state
const vector<double>& dt, // temporal grid of state
const vector<double>& dt_new, // new temporal grid for state
vector<typename PDELab::Backend::Vector<GFS,double>>& w) // extended state
{

	/* time points of new temporal grid */
	vector<double> t_new(dt_new.size()+1);
	t_new[0]=0; 
	for(int i=0; i < dt_new.size(); i++)
		t_new[i+1]=t_new[i]+dt_new[i];

	/* time points of old temporal grid */
	vector<double> t(dt.size()+1);
	t[0]=0; 
	for(int i=0; i < dt.size(); i++)
		t[i+1]=t[i]+dt[i];

	for(int i=0; i < dt_new.size()+1; i++){
		// first time point of old temporal grid that is not less than t_new[i]
		auto start=lower_bound(t.begin(),t.end(), t_new[i],[](const double& a, double value){return a < value-1e-10;});
		int index=distance(t.begin(),start);
		// set value of extended state on (t_new[i-1],t_new[i]]
		w[i]=y[index]; 	
	}

}

/***************************************************************************
 *
 * calculation of G*S*(Sf-yd) and Sf
 *
 ***************************************************************************/
template<class GFS, class ControlVector>
void RHScontribution(
const GFS& gfs, // grid function space
const vector<double>& dt, // temporal grid
HeatProblemInterface<typename GFS::Traits::GridView,ControlVector::block_type::dimension>& heat, // for G*S*Sf
AdjointProblemInterface<typename GFS::Traits::GridView>& YdProblem, // for G*S*yd
ControlVector& Gstar_eta, // G*S*(Sf-yd)
vector<typename PDELab::Backend::Vector<GFS,double>>& g1,  // Sf
vector<typename PDELab::Backend::Vector<GFS,double>>& g2) // S*yd
{
	typedef typename PDELab::Backend::Vector<GFS,double> vtype;

	g1=vector<vtype>(dt.size()+1,vtype(gfs));
	g2=vector<vtype>(dt.size()+1,vtype(gfs));
	Gstar_eta.resize(dt.size());
	Gstar_eta=0;

	// inhomogeneous PDE data
	if(!heat.isHomogeneous())
	{
		// calculation of Sf
		heatdriver(gfs,heat,dt,g1);
		AdjointProblem<GFS> adjHeat(dt,g1,gfs);
		// calculation q=S*Sf 
		vector<vtype> q; 
		adjointdriver(gfs,adjHeat,q);  
	      
    		// contibution of G*S*Sf=\int Psi(x)q(t,x) dx 
    		Gstar(gfs,heat,q,Gstar_eta);
	}

	// desired temperature contribution
	ControlVector Gstar_yd;
	YdProblem.setdeltaT(dt);
	{
		// calculation eta=S*yd
		adjointdriver(gfs,YdProblem,g2);
        	
		// contibution of G*S*yd=\int Psi(x)eta(t,x) dx 
		Gstar(gfs,heat,g2,Gstar_yd);
	}
	Gstar_eta-=Gstar_yd;
}


/***************************************************************************
 *
 * calculation of G*S*SGu_fix
 *
 ***************************************************************************/
template<class GFS, class ControlVector>
void RHScontribution(
const GFS& gfs, // grid function space
const vector<double>& dt, // temporal grid
CGProblemInterface<typename GFS::Traits::GridView,ControlVector::block_type::dimension>& cgproblem, // for SGu
ControlVector& Gstar_q)  // G*S*SGu_fix
{
	typedef typename PDELab::Backend::Vector<GFS,double> vtype;

	Gstar_q.resize(dt.size());
	Gstar_q=0;

	// q=SGu
	vector<vtype> y;
	heatdriver(gfs,cgproblem,dt,y);
	AdjointProblem<GFS> adjointProblem(dt,y,gfs);
	// w=S*SGu
	vector<vtype> w;  
	adjointdriver(gfs,adjointProblem,w);
	// Gstar_q= G*S*SGu
	Gstar(gfs,cgproblem,w,Gstar_q); 
}

#endif //ifned HELPFUNCTIONS_HH
