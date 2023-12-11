// Copyright © Alexandra Grütering, see file LICENSE.md in module root
// License: GLP-3.0 or later
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

// C, C++ includes
#include<ctime> 
#include<math.h>
#include<iostream>
#include<vector>
#include <stdlib.h> 

// dune-common includes
#include<dune/common/parametertreeparser.hh>
#include<dune/common/fvector.hh> 
#include<dune/common/fmatrix.hh> 

// dune-grid includes
#include<dune/grid/utility/structuredgridfactory.hh>
#include<dune/grid/onedgrid.hh>
// dune-istl includes
#include<dune/istl/bvector.hh> 
#include<dune/istl/bcrsmatrix.hh>
// dune pdelab includes
#include<dune/pdelab/finiteelementmap/pkfem.hh>
#include<dune/pdelab/constraints/common/constraints.hh>
#include<dune/pdelab/constraints/common/constraintsparameters.hh>
#include<dune/pdelab/constraints/conforming.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include<dune/pdelab/backend/istl.hh>
#include<dune/pdelab/common/quadraturerules.hh>
#include<dune/pdelab/finiteelement/localbasiscache.hh>

#include <dune/OCP/probleminterface.hh>
#include <dune/OCP/Dmax.hh>
#include <dune/bnb/OCPbranch.hh>
#include <dune/bnb/OCPDmax.hh>
#include <dune/bnb/breadthfirst.h>

using namespace Dune; 
using namespace std;
using namespace BranchAndBound;

template<class GridView>
class HeatProblem : public HeatProblemInterface<GridView,1>
{
protected:

	typedef HeatProblemInterface<GridView,1> BaseType;

	using typename BaseType::E;
	using typename BaseType::I; 
	using typename BaseType::X; 
	
	using BaseType::t;

public:

	HeatProblem() : BaseType() {}

	FieldVector<double,1> Psi (const E& e, const X& x) const
  	{
		int dimDomain = x.size(); 
		auto global = e.geometry().global(x);   		
		double value = exp(global[0])*sin(M_PI*global[0]);
		return value; 
  	}
	
};


/* problem for calculating Sd in CG algorithm */
template<class GridView>
class CGProblem : public CGProblemInterface<GridView,1>
{
protected:

	typedef CGProblemInterface<GridView,1> BaseType;

	using typename BaseType::ControlVector;
	using typename BaseType::E;
	using typename BaseType::I; 
	using typename BaseType::X; 
	
	using BaseType::t;

public:

	CGProblem(const ControlVector& w, double dt_, int timesteps) : BaseType(w,dt_,timesteps) {}

	CGProblem(const ControlVector& w, vector<double> dt_) : BaseType(w,dt_) {}

	FieldVector<double,1> Psi (const E& e, const X& x) const
  	{
		int dimDomain = x.size(); 
		auto global = e.geometry().global(x);   		
		double value = exp(global[0])*sin(M_PI*global[0])+0.5;
		return value; 
  	}

	
};

template<class GFS>
class YDProblem : public AdjointProblem<GFS>
{
protected:
	typedef   AdjointProblem<GFS> BaseType;
	typedef BlockVector<FieldVector<double,1>> Vector;

	using typename  BaseType::Y;
	using typename BaseType::E;
	using typename BaseType::I; 
	using typename BaseType::X; 

	using BaseType::t;
	using BaseType::T;
	using BaseType::dt;
	using BaseType::gfs;
	using BaseType::y;
	

public: 
	YDProblem(double endTime_, int timesteps, vector<Y>& y_,const GFS& gfs_) : BaseType(endTime_,timesteps,y_,gfs_) {}

	YDProblem(vector<double> dt_, vector<Y>& y_,const GFS& gfs_) : BaseType(dt_,y_,gfs_) {}
	

	double rhs(const E& e, const X& x, double time) const
	{
		typedef typename PDELab::DiscreteGridFunction<GFS,Y> DGF;
    		typedef typename DGF::Traits::RangeType RangeType;

		int l=0;
		double clock=0; 
		while(clock < time-1e-10){
			clock+=0.003125;
			l++; 
		}
		DGF ydgf(gfs,y[l]);
		RangeType value; 
    		ydgf.evaluate(e,x,value); 
		return value;	
	}

};


/* Main program with grid setup */
int main(int argc, char** argv)
{
	// Maybe initialize Mpi
    	Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);
    	if(Dune::MPIHelper::isFake)
      		std::cout<< "This is a sequential program." << std::endl;
   	else
      		std::cout << "Parallel code run on "
                	<< helper.size() << " process(es)" << std::endl;

   	//===============================================================
       // Make parameter structure 
       //===============================================================
    	// open ini file (contains grid information, filenames for output)
    	Dune::ParameterTree ptree;
    	Dune::ParameterTreeParser ptreeparser;
        std::stringstream file;
        file << BNB_SOURCE_PATH << "BnBTest.ini"; 
    	ptreeparser.readINITree(file.str(),ptree);
    	ptreeparser.readOptions(argc,argv,ptree);


        //===============================================================
        // make 1D grid
        //===============================================================
	  const int dim=1;
        typedef Dune::OneDGrid Grid;	
        Dune::FieldVector<double,dim> L;
	for(int i=0; i < dim; i++){
		stringstream name; 
		name << "grid.structured.L" << i; 
		L[i] = ptree.get<double>(name.str(),(double)0.0);
	}
	Dune::FieldVector<double,dim> U;
	for(int i=0; i < dim; i++){
		stringstream name; 
		name << "grid.structured.U" << i; 
		U[i] = ptree.get<double>(name.str(),(double)1.0);
	}
        std::array<unsigned int,dim> N;
	for(int i=0; i < dim; i++){
		stringstream name; 
		name << "grid.structured.N" << i; 
		N[i] = ptree.get<unsigned int>(name.str(),(unsigned int) 100);
	} 
	std::shared_ptr<Grid> gridp = StructuredGridFactory<Grid>::createSimplexGrid(L,U,N);
	
   	typedef Grid::LeafGridView GV;
    	GV gv=gridp->leafGridView();

    	typedef Dune::PDELab::ConformingDirichletConstraints CON; // conforming dirichlet constraints
	typedef PDELab::ISTL::VectorBackend<> VBE; // ISTL backend to storage Dofs
	typedef typename GV::ctype DF;
	typedef typename PDELab::PkLocalFiniteElementMap<GV,DF,double,1> FEM; 
	HeatProblem<GV> heat; // heat problem with boundary data
	typedef typename PDELab::GridFunctionSpace<GV,FEM,CON,VBE> GFS; // grid function space
	typedef typename PDELab::Backend::Vector<GFS,double> vtype;
	typedef AdjointProblem<GFS> AdjProblem;
	typedef BlockVector<FieldVector<double,1>> Vector; // type of control vector
	
	
	double T=ptree.get<double>("problem.T",(double)1.0); 
	int timesteps= ptree.get<int>("problem.timesteps",(int)100);
	vector<double> dt(timesteps,(double) T/timesteps);
	int jumps=ptree.get<int>("problem.jumps",(int) 5);; 
	int smax=ptree.get<int>("problem.smax",(int) 2);

	// random generated jump points
	vector<int> jpoints(jumps); 
	bool stop;
	// set up seed
	unsigned int seed=time(0);
        srand(seed);
	
	// generate random jump points
	for (int i = 0; i < jumps; i++){
        	bool same;
        	do
       		{
           		same = false;
            		jpoints[i] = rand() % (319);
            		// Check if the newly generated number is a duplicate
            		for (int j = 0; j < i;j++){
               			if (jpoints[i] == jpoints[j]){
                    			same = true;
                    			break;
                		}
            		}
        	} while (same);
    	}
	sort(jpoints.begin(),jpoints.end()); 

	// calcuate ud
	vector<double> tp(jumps);
	for(int i=0; i < jumps; i++)
		tp[i]=(jpoints[i]+1)*(T/320.0);
	
	double h=T/320.0;
	// create piecewise constant function with jump points tp
	Vector ud(320,0);
	int j=0;
	for(int i=0; i < jumps; i++){
		while( (j+1)*h < tp[i]+1e-10){
			ud[j]+=i%2;
			j++;
		}
	}
	for(int k=j; k < 320; k++){
		ud[k]+=jumps%2;
	}

	// output jump times
        std::cout << "Jump times" << std::endl;
        for(int i=0; i < jumps; i++)
		cout << "t_" << i << ":" << tp[i]<< endl;
	
	// calculation of desired temperature yd(t,x)
	FEM fem(gv); 
	GFS gfs(gv,fem);
	vector<vtype> yd; 
	vector<double> tgrid(320,h);
	CGProblem<GV> problem(ud,tgrid);
	heatdriver(gfs,problem,tgrid,yd);
	YDProblem<GFS> AdjDesired(dt,yd,gfs);
        
	stringstream filename;
	filename << BNB_OUTPUT_PATH << OUTPUTDIR << "/settings.txt";
	fstream stream; 
	stream.precision(15);
	stream.open(filename.str(), ios::out | ios::app);
	stream << "seed=" << seed << endl;
	stream << "jumps=" << jumps << endl;
	stream.close();
        
	try{
		Dmax D(smax);
		OCPDmax<GV,CON,CGProblem<GV>,1> master(gv, /* spatial grid */
						       heat, /* heat problem */ 
						       AdjDesired, /* desired temperature problem */
						       tgrid, /* temporal grid of yd */
						       D, /* switching polytope */
						       alpha, /* Tikhonov term */
						       (sqrt(5)+1)/2, /* penalty parameter cutting planes in ADMM */
						       alpha, /* penalty parameter box constraints in ADMM */
						       dt, /* start temporal grid */
	                 	    		       40, /* maximum iterations of linear solver */
						       1e-10, /* reduction factor for linear solver */
						       1e-3, /* absolute tolerance in ADMM */
						       1e-3, /* relative tolerance in ADMM */
						       1e-5, /* error tolerance in ADMM */
						       0.5,  /* grid refinement factor */
						       0.02, /* objective reduction factor */
						       0.02);
		master.writeParameterSetting(); 
		EnumerationStrategy* strategy= new BreadthFirst();
		master.setEnumerationStrategy(strategy);
		master.optimize();
	}
	catch (Dune::Exception &e){
		std::cerr << "Dune reported error: " << e << std::endl;
  	}
  	catch(std::exception & e){
    		std::cerr << "Error: " << e.what() << std::endl;
 	}
  	catch (...){
    		std::cerr << "Unknown exception thrown!" << std::endl;
  	}

	return 0;
}
