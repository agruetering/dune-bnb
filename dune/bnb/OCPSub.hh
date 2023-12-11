// Copyright © Alexandra Grütering, see file LICENSE.md in module root
// License: GLP-3.0 or later
#ifndef OCPSUB_HH
#define OCPSUB_HH

#include<vector>
#include<dune/bnb/subproblem.h>

using namespace BranchAndBound;
using namespace std;

/* structure for output data of subproblem
*/
struct OutputData{
	int status;
	/* -1: default 
	    0: feasible solution for subproblem found
	    1: feasible integer solution for subproblem found
	*/
	double time_feasibility;
	double time_updateFixData;
	double time_separation;
	double time_solve;
	vector<int> iter_admm; // admm iterations 
	vector<double> dual_bound; // lower bounds within outer approximation
	vector<double> time_admm_solve; // admm solver time
	vector<double> time_outer_separation; // separation time
	int iter_cut; // number of outer approximations iterations
	OutputData() {
		status=-1;
		time_feasibility=0;
		time_updateFixData=0;
		time_separation=0;
		time_solve=0;
		iter_admm=vector<int>(); 
		dual_bound=vector<double>();
		time_admm_solve=vector<double>();
		time_outer_separation=vector<double>();
		iter_cut=0;
	}

	void writeInfo(int count,int anzFixPoint, int anzFixIndex,double min_dt, double bound, bool fathom)
	{
		stringstream dir;
		dir << BNB_OUTPUT_PATH << OUTPUTDIR << "sub" << count;
		mkdir(dir.str().c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

		string filename;
		fstream file; 
		file.precision(std::numeric_limits<double>::digits10+2);
		// write grid to file with gmswrite
		filename=dir.str()+"/statistics.txt";
		file.open(filename, ios::out | ios::trunc);

		file << "\n[Input]" << endl; 
		file << "anzFixPoint:" << anzFixPoint << endl;
		file << "anzFixIndex:" << anzFixIndex << endl;

		file << "\n[Output]" << endl; 
		file << "t_fcheck:" << time_feasibility << endl; 
		file << "t_updFix:" << time_updateFixData<<endl; 
		file << "t_sep:" << time_separation << endl; 
		file << "iter_cut:" << iter_cut << endl;
		file << "dual_bound:" << bound << endl;
		file << "TimeResolution:" << min_dt << endl; //smallest grid cell
		
		file << "\nADMMIterations" << endl;
		for(int i=0; i < iter_admm.size(); i++)
			file << "a_" << i << ":" << iter_admm[i] << endl;
		file << "a_total:" << accumulate(iter_admm.begin(),iter_admm.end(), 0)<< endl;
		file << "\nADMMTime" << endl;
		for(int i=0; i < time_admm_solve.size(); i++)
			file << "o_" << i << ":" << time_admm_solve[i] << endl;
		file << "o_total:" << accumulate(time_admm_solve.begin(),time_admm_solve.end(), 0.0)<< endl;

		file << "\nTimeSeparation" << endl;
		for(int i=0; i < time_outer_separation.size(); i++)
			file << "s_" << i << ":" << time_outer_separation[i] << endl;
		file << "s_total:" << accumulate(time_outer_separation.begin(),time_outer_separation.end(), 0.0)<< endl;

		file << "\ntime_total:" << time_solve << endl;

		file << "\nDualBounds" << endl;
		for (int i=0; i< dual_bound.size(); i++)
			file << "d_" << i << ":" << dual_bound[i] << endl;
		
		file << "\n[Status]" << endl;
		file << "status:" <<status << endl;
		file << "fathomed:" << fathom << endl;

		file.close();
	}
		

};

 /* Optimal control subproblem class */
template<class OCPMaster>
class OCPsub : public Subproblem
{
public: 
	typedef OCPMaster Master;
	typedef typename OCPMaster::GV GV;
	typedef typename GV::Grid Grid; 
	typedef typename OCPMaster::FEM FEM; 
	typedef typename OCPMaster::GFS GFS; 
	typedef typename OCPMaster::vtype vtype;
	typedef typename OCPMaster::ControlVector ControlVector;
	typedef typename OCPMaster::Matrix Matrix; 
	
	using Vector = BlockVector<FieldVector<double,1>>;

	enum{n=OCPMaster::n};
private: 
	// fixed control variable indices
	vector<pair<int,int>> fixIndices;

	// spatial discretization
	GV gv; 
	FEM fem; 
	GFS gfs; 
	// temporal discretization
	vector<double> dt;	
	bool grid_changed; // flag if temporal gird was changed 

	// desired temperature and inhomogeneuous PDE data
	vector<vtype> Sf; // encode Sf on given discretization (needed for error estimator)
	vector<vtype> Staryd; // encode S*yd on given discretization (needed for error estimator)
	ControlVector Gstar_eta; // encodes G*S*(Sf-yd)
	

	ControlVector Gstar_q; // encodes G*S*SGu_fix with u_fix extended vector of fixed control variables

	// cutting planes 
	Matrix A; // coefficient matrix
	Vector b; // right-hand side
	Vector v; // auxiliary variable for cutting planes
	Vector mu; // lagrange multiplicator 

	// box constraints
	Vector w; // auxiliary variable for box constraints
	Vector nu; // lagrange multiplicator for box constraints
	
	ControlVector u; // control

public: 
	/* constructor and destructor */
	OCPsub(Master* pMaster_,GV grid,vector<double> dt_) : 
	Subproblem(pMaster_), fixIndices(vector<pair<int,int>>()),
	dt(dt_), gv(grid), fem(FEM(gv)), gfs(GFS(gv,fem)), 
	Sf(vector<vtype>(dt.size()+1,vtype(gfs))),
	Staryd(vector<vtype>(dt.size()+1,vtype(gfs))),
	Gstar_eta( ControlVector(dt.size()) ), 
	Gstar_q(ControlVector(dt.size())), 
	A(Matrix()), 
	b(Vector()),
	v(Vector()), 
	mu(Vector()),
	w(dt.size()*n),
	nu(dt.size()*n),
	u(ControlVector(dt.size(),0)), 	
	grid_changed(0)
	{}

	OCPsub(Subproblem* pParent_,GFS gfs_,vector<pair<int,double>> tau_, vector<double> d_) : 
	Subproblem(pParent_,tau_,d_), gv(gfs_.gridView()), fem(FEM(gv)), gfs(GFS(gv,fem))
	{
		((OCPsub<OCPMaster>*)pParent_)->copyData(this);
	}

	~OCPsub(){}

	/* subproblem construction with input data */
	OCPsub(Master* pMaster_,GV grid,vector<double> dt_ ,vector<pair<int,double>> tau_,vector<double> d_,vector<pair<int,int>>index, ControlVector u_, Matrix A_, Vector b_, Vector v_, Vector mu_, Vector w_, Vector nu_) : 
	Subproblem(pMaster_,tau_,d_), fixIndices(index),
	dt(dt_), gv(grid), fem(FEM(gv)), gfs(GFS(gv,fem)), 
	Sf(vector<vtype>(dt.size()+1,vtype(gfs))),
	Staryd(vector<vtype>(dt.size()+1,vtype(gfs))),
	Gstar_eta( ControlVector(dt.size()) ),
	Gstar_q(ControlVector(dt.size())), 
	A(A_), 
	b(b_), 
	v(v_),
	mu(mu_),
	w(w_),
	nu(nu_),
	u(u_), 
	grid_changed(0)
	{}

	/* inheritance from father node */
	void copyData(OCPsub<OCPMaster>* sub)
	{
		sub->setData(fixIndices,gfs,dt,Sf,Staryd,Gstar_eta,Gstar_q, A,b,v,mu,w,nu,u);
	} 

	/* set all data of subproblem */
	void setData(vector<pair<int,int>>& index, GFS& gfs_, vector<double>& dt_, vector<vtype>& Sf_, vector<vtype>& g, ControlVector& Gstar_eta_, ControlVector& Gstar_q_, Matrix& A_,Vector&b_, Vector& v_, Vector& mu_, Vector& w_, Vector& nu_, ControlVector&u_)
	{
		fixIndices=index;

		// discretization
		gv=gfs_.gridView(); 
		fem=FEM(gv);
		gfs=GFS(gv,fem);
		dt=dt_;

		// desired temperature and inhomogeneuous PDE data
		Sf=vector<vtype>(dt.size()+1,vtype(gfs));
		Sf=Sf_;
		Staryd=vector<vtype>(dt.size()+1,vtype(gfs));
		Staryd=g;
		Gstar_eta=Gstar_eta_;
		Gstar_q=Gstar_q_;

		// cutting planes
		A=A_;
		b=b_; 
		v=v_;
		mu=mu_; 

		// box constraints
		w=w_; 
		nu=nu_;

		// previous control
		u=u_;

	
	}
	/* set solution of subproblem after outer approximation */ 
	void setData(Matrix& A_,Vector&b_, Vector& v_, Vector& mu_, Vector& w_, Vector& nu_, ControlVector&u_)
	{
		// cutting planes
		A=A_;
		b=b_; 
		v=v_;
		mu=mu_; 

		// box constraints
		w=w_; 
		nu=nu_;

		// control
		u=u_;

	}
	/* set suproblem data after grid refinement */ 
	void setData(vector<pair<int,int>>& index, GFS& gfs_, vector<double>& dt_, vector<vtype>& Sf_,vector<vtype>& g, ControlVector& Gstar_eta_, ControlVector& Gstar_q_, Matrix& A_, Vector& w_, Vector& nu_, ControlVector&u_)
	{
		fixIndices=index;

		// discretization
		gv=gfs_.gridView(); 
		fem=FEM(gv);
		gfs=GFS(gv,fem);
		dt=dt_;
		grid_changed=1;

		// desired temperature and inhomogeneuous PDE data
		Sf=vector<vtype>(dt.size()+1,vtype(gfs));
		Sf=Sf_;
		Staryd=vector<vtype>(dt.size()+1,vtype(gfs));
		Staryd=g;
		Gstar_eta=Gstar_eta_;
		Gstar_q=Gstar_q_;


		// cutting plane matrix
		A=A_; 

		// box constraints
		w=w_; 
		nu=nu_;

		// previous control
		u=u_;

	}
	/* set data of root node */
	void setData( vector<vtype>& Sf_,vector<vtype>& g, ControlVector& Gstar_eta_)
	{
		// desired temperature and inhomogeneuous PDE data
		Gstar_eta=Gstar_eta_;
		Sf=vector<vtype>(dt.size()+1,vtype(gfs));
		Sf=Sf_;
		Staryd=vector<vtype>(dt.size()+1,vtype(gfs));
		Staryd=g;
	}


	/* Setter */
	void setFixData(ControlVector& Gstar_z){Gstar_q=Gstar_z;}
	void setControl(ControlVector& u_) { u=u_;}
	void setFixIndices(vector<pair<int,int>>& index) { fixIndices=index;}
	void GridChanged() {grid_changed=1;}

	/* Getter */
	GFS getGFS() {return gfs;}
	vector<vtype> getStaryd() {return Staryd;}
	vector<double> getdeltaT() {return dt;}
	ControlVector getControl() {return u;}
	vector<pair<int,int>> getFixIndices() {return fixIndices;}
	ControlVector getGeta() {return Gstar_eta;}
	bool getGridChanged() {return grid_changed;}
	void getLambda(vector<pair<int,int>>& index,vector<vtype>& Sf_, ControlVector& Gy, ControlVector& Gq, Matrix& A_, Vector& b_,Vector& v_, Vector& mu_, Vector& w_, Vector& nu_)
	{
		// fixed control variable indices
		index=fixIndices;

		// desired temperature and inhomogeneuous PDE data
		Gy=Gstar_eta; 
		Gq=Gstar_q;
		Sf_=vector<vtype>(dt.size()+1,vtype(gfs));
		Sf_=Sf;

		// cutting planes 
		A_=Matrix(A);
		b_=b; 
		v_=v;
		mu_=mu; 

		// box constraints 
		w_=w; 
		nu_=nu;
	}

	/* extension of cutting planes*/
	void extendCut(double* cut, double& s){
	
		// add row to A
		Matrix A_new=Matrix(A.N()+1,dt.size()*n,Matrix::row_wise); 
		for(auto row=A_new.createbegin(); row!=A_new.createend(); ++row){
                	if (row.index()< A.N()){
				for( auto it=A[row.index()].begin(); it!=A[row.index()].end(); ++it) 
					row.insert(it.index());
                	}
               		else{
                    		for(int k=0; k<n*dt.size(); k++){
                        		if(cut[k]!=0) 
                            			row.insert(k);  
                		}
            		} 
		}
		for(int i=0; i< A.N(); i++) 
                	A_new[i]=A[i];
            	int end = A.N();
            	for(int k=0; k<n*dt.size(); k++){
			if(cut[k]!=0)
                    		A_new[end][k]=cut[k]; 
		}
		A=A_new;

		// add value s to b
		b.resize(end+1);
		b[end]=s;

		// extend lagrange multiplicator
		mu.resize(end+1);

		// extend auxiliary variable 
		v.resize(end+1);
	}
		
	void resetBox()
	{
		w.resize(dt.size()*n-fixIndices.size());
		w=0;
		nu.resize(dt.size()*n-fixIndices.size());
		nu=0;
	}

	// write info of sub problem to file
	void writeInfo(int count)
	{
		stringstream dir;
		dir << BNB_OUTPUT_PATH << OUTPUTDIR << "sub" << count;
		mkdir(dir.str().c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

		vector<pair<int,double>> tau=this->getTau(); 
		vector<double> d=this->getFix();

		string filename;
		fstream file; 
		file.precision(std::numeric_limits<double>::digits10+2);
		#ifdef BNBTREE
		// write grid to file with gmswrite
		filename=dir.str()+"/grid.msh"; 
		GmshWriter<GV> writer(gv);
  		writer.setPrecision(10);
  		writer.write(filename);
		// write dt to file 
		filename=dir.str()+"/dt.txt";
		file.open(filename, ios::out | ios::trunc);
		for(int i=0; i<dt.size(); i++)
			file << dt[i] <<endl; 
		file.close();
		// write u to file 
		filename=dir.str()+"/u.txt";
		storeMatrixMarket<ControlVector>(u,filename);
		// write A to file 
		filename=dir.str()+"/A.txt";
		storeMatrixMarket<Matrix>(A,filename);
		// write b to file 
		filename=dir.str()+"/b.txt";
		storeMatrixMarket<Vector>(b,filename);
		// write mu to file 
		filename=dir.str()+"/mu.txt";
		storeMatrixMarket<Vector>(mu,filename);
		// write v to file
		filename=dir.str()+"/v.txt";
		storeMatrixMarket<Vector>(v,filename);
		// write w to file
		filename=dir.str()+"/w.txt";
		storeMatrixMarket<Vector>(w,filename);
		// write nu to file
		filename=dir.str()+"/nu.txt";
		storeMatrixMarket<Vector>(nu,filename);
		// write fixIndices to file 
		filename=dir.str()+"/fixIndices.txt";
		file.open(filename, ios::out | ios::trunc);
		for(int i=0; i<fixIndices.size(); i++)
			file << fixIndices[i].first << " " << fixIndices[i].second << endl;  
		file.close();
		#endif
		// write tau to file 
		filename=dir.str()+"/tau.txt";
		file.open(filename, ios::out | ios::trunc);
		for(int i=0; i<tau.size(); i++)
			file  << tau[i].first << " " << tau[i].second  << endl; 
		file.close();
		// write c to file 
		filename=dir.str()+"/d.txt";
		file.open(filename, ios::out | ios::trunc);
		for(int i=0; i<d.size(); i++)
			file << d[i] << endl; 
		file.close();
	}
};

#endif
