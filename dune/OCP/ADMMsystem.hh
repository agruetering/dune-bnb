// dune istl includes
#include <dune/istl/operators.hh> 
#include<dune/istl/vbvector.hh>
#include <dune/istl/solvers.hh> 
#include <dune/istl/umfpack.hh>
#include <dune/istl/matrix.hh>


using namespace Dune; 
using namespace std; 

/* Linear Operator without cutting planes for ADMM method
*/
template<class MType>
class LinOperator : public LinearOperator< BlockVector<FieldVector<double,1>>, BlockVector<FieldVector<double,1>> >
{
protected:
	typedef BlockVector<FieldVector<double,1>> X;
	
	typedef MType Master;
	typedef typename Master::GFS GFS; 
    	typedef typename Master::vtype vtype; // type of dof vector
	typedef typename Master::ControlVector ControlVector; // type of control vector
	typedef typename Master::Matrix Matrix; // matrix type

	typedef typename Master::AdjProblem AdjProblem;
	typedef typename Master::CGProblem CGProblem;
	// using Vector = BlockVector<FieldVector<double,1>>;
	enum{n=Master::n};

public:
	LinOperator(const GFS& gfs_,const vector<double>& dt_,double alpha_,double beta_,const vector<pair<int,int>>& I_) : gfs(gfs_), dt(dt_),alpha(alpha_),beta(beta_), I(I_) {};

	/* y=Cx */
	void apply(const X& x, X& y) const 
	{
		// clear y
		y=0.0;
		
		// fill x with zeros
		ControlVector z(dt.size(),0);
		for (int i=0;i<I.size(); i++ )
			z[I[i].second][I[i].first]=x[i];
 
		// calculation of w=SGx;
		CGProblem cgproblem(z,dt); 
		vector<vtype> w; 
		heatdriver(gfs,cgproblem,dt,w);
    
    	
		// calculation of v=S*SGx;
    		AdjProblem adjointProblem=AdjProblem(dt,w,gfs); 
    		vector<vtype> v; 
		adjointdriver(gfs, adjointProblem,v);  

    		
    		// y= (alpha+beta) x + \Chi_I(G*S*SGx) 
		y=x;
        	y*=(alpha+beta);
    		// contibution of \Chi_I(G*S*SGx) = \Chi_I(\int Psi(x)v(t,x) dx )
		ControlVector Gstar_z(dt.size()+1); 
    		Gstar(gfs,cgproblem,v,Gstar_z);
		for (int i=0;i<I.size(); i++ ){
			y[i]+=Gstar_z[I[i].second][I[i].first];
			y[i]*=dt[I[i].second];
		}


	}

	/* y= y+tau*Cx */
	void applyscaleadd(field_type tau, const X& x, X& y) const
	{
		X z=y; // auxiliary variable to save result Cx
		apply(x,z); 
		
		z*=tau; 
		y+=z;	
	}

	SolverCategory::Category category() const{
		return SolverCategory::sequential;
	}


private: 
	const GFS& gfs; // spatial grid function space
	const vector<double>& dt; // temporal grid
	double alpha; // Thikonov term
	double beta; // penaly term box constraints
	const vector<pair<int,int>>& I; // free control variables
};



/* Linear Operator with cutting planes for ADMM method
*/
template<class MType>
class LinCutOperator : public LinearOperator< BlockVector<FieldVector<double,1>>, BlockVector<FieldVector<double,1>> >
{
protected:
	typedef BlockVector<FieldVector<double,1>> X;
	
	typedef MType Master;
	typedef typename Master::GFS GFS; 
    	typedef typename Master::vtype vtype; // type of dof vector
	typedef typename Master::ControlVector ControlVector; // type of control vector
	typedef typename Master::Matrix Matrix; // matrix type

	typedef typename Master::AdjProblem AdjProblem;
	typedef typename Master::CGProblem CGProblem;
	//using Vector = BlockVector<FieldVector<double,1>>;
	enum{n=Master::n};

public:
	LinCutOperator(const GFS& gfs_,const Matrix& A_,const vector<double>& dt_,double alpha_,double rho_,double beta_, const vector<pair<int,int>>& I_) : gfs(gfs_), A(A_),dt(dt_),alpha(alpha_),rho(rho_),beta(beta_), I(I_) {};

	/* y=Cx */
	void apply(const X& x, X& y) const 
	{
		// clear y
		y=0.0;
		
		// fill x with zeros
		ControlVector z(dt.size(),0);
		for (int i=0;i<I.size(); i++ )
			z[I[i].second][I[i].first]=x[i];
 
		// calculation of w=SGx;
		CGProblem cgproblem(z,dt); 
		vector<vtype> w; 
		heatdriver(gfs,cgproblem,dt,w);
    
    	
		// calculation of v=S*SGx;
    		AdjProblem adjointProblem=AdjProblem(dt,w,gfs); 
    		vector<vtype> v; 
		adjointdriver(gfs, adjointProblem,v);  

    		
    		// y= (alpha+beta) x + \Chi_I(G*S*SGx)+\Chi_I(rho*A^TAx) 
		y=x;
        	y*=(alpha+beta);
    		// contibution of \Chi_I(G*S*SGx) = \Chi_I(\int Psi(x)v(t,x) dx )
		ControlVector Gstar_z(dt.size()+1); 
    		Gstar(gfs,cgproblem,v,Gstar_z);
		for (int i=0;i<I.size(); i++ ){
			y[i]+=Gstar_z[I[i].second][I[i].first];
			y[i]*=dt[I[i].second];
		}

		// contribution of \Chi_I(rho*A^TAx) 
		// store x filled with zeros appropiately 
		X h(dt.size()*n); 	
		for(int i=0; i < I.size(); i++) 
			h[I[i].second*n+I[i].first]=x[i];
		// calculate r=Ax 
		X r(A.N()); 
		A.mv(h,r);
		// calculate h=rho*A^TAx
		h=0; 
		A.usmtv(rho,r,h); 
		for (int i=0;i<I.size(); i++ ){
			y[i]+=h[I[i].second*n+I[i].first];	
		}
	}
	

	/* y= y+tau*Cx */
	void applyscaleadd(field_type tau, const X& x, X& y) const
	{
		X z=y; // auxiliary variable to save result Cx
		apply(x,z); 
		
		z*=tau; 
		y+=z;	
	}

	SolverCategory::Category category() const{
		return SolverCategory::sequential;
	}


private: 
	const GFS& gfs; // spatial grid function space
	const Matrix& A; // cutting plane matrix
	const vector<double>& dt; // temporal grid
	double alpha; // Thikonov term
	const vector<pair<int,int>>& I; // free control variables
	double rho; // penalty term cutting plane
	double beta; // penaly term box constraints
};




/* Inverse Operator for Preconditioning */ 
template<class M,class X,class Y> 
class InverseP : public InverseOperator<X,Y>{
  public:	
    	typedef M matrix_type;
 	typedef X domain_type;
    	typedef Y range_type;

	// Constructor 
	InverseP(const M& P_): P(P_) {}; 

	// apply inverse operator
	void apply(X& x,Y& b, InverseOperatorResult& res){
		// call UMFPack
		UMFPack<M> Solver(P); 
		Solver.apply(x,b,res);
	}

	// apply inverse operator, with given convergence criteria
	void apply(X& x,Y& b,double reduction, InverseOperatorResult& res){
		// call UMFPack
		UMFPack<M> umfpackSolver(P); 
		umfpackSolver.apply(x,b,reduction,res);
	}

	// Solver category 
	SolverCategory::Category category() const{
		return SolverCategory::sequential;
	}
  private: 
	const M& P; // Preconditioning matrix
};




