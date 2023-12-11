#include <dune/bnb/OCPMaster.hh>
#include <dune/istl/io.hh>

/**************************************************************************
 *
 * set up desired temperature and inhomogeneuous PDE data
 *
 **************************************************************************/
template<class GridView, class Constraints, class Problem, int k, int number>
void OCPMaster<GridView,Constraints,Problem,k,number>::initialize(Subproblem* sub)
{
	GFS gfs=((SubType*)sub)->getGFS();
	vector<double> dt=((SubType*)sub)->getdeltaT();
	/* G*S*(Sf-yd) and Sf */
	ControlVector Gstar_eta;
	vector<vtype> Sf; 
	vector<vtype> Staryd;
	RHScontribution(gfs,dt,this->heat,this->YdProblem,Gstar_eta,Sf,Staryd);

	((SubType*)sub)->setData(Sf,Staryd,Gstar_eta);
}

/************************************************************************** 
 *
 * update subproblem data due to changed fixed control vector
 *
 **************************************************************************/
template<class GridView, class Constraints, class Problem, int k, int number>
void OCPMaster<GridView,Constraints,Problem,k,number>::updateFixData(Subproblem* sub,const vector<pair<int,int>>& index, const GFS& gfs, const vector<double>& dt, const ControlVector& v)
{
	/* calculate G*S*Su_fix */
	// set up extended vector of fixed control variables
	ControlVector ufix(dt.size());
	bool null=true; // indicate whether fixed values are all zero
	for(int i=0; i < index.size(); i++){
		ufix[index[i].second][index[i].first]=v[index[i].second][index[i].first];
		if(ufix[index[i].second][index[i].first]!=0)
			null=false;
	}
	// z= G*S*Su_fix
	ControlVector Gstar_z(dt.size());
	vector<vtype> z;
	if(!null){
		CGProblem cgproblem(ufix,dt); 
		RHScontribution(gfs,dt,cgproblem,Gstar_z); 
	}
	((SubType*)sub)->setFixData(Gstar_z);	 

	
}

/**************************************************************************
 *
 * method to update data after grid refinement 
 *
 **************************************************************************/ 
template<class GridView, class Constraints, class Problem, int k, int number>
void OCPMaster<GridView,Constraints,Problem,k,number>::updateData(vector<pair<int,int>>& fixIndex,
const GFS& gfs,const vector<double>& dt_new, vector<double>& dt, 
vector<vtype>& Sf, vector<vtype>& Staryd, ControlVector& Gstar_eta, ControlVector& Gstar_q,Matrix& A, Vector& w, Vector& nu, ControlVector& u)
{
	/* time points of new temporal grid */
	vector<double> t_new(dt_new.size());
	t_new[0]=dt_new[0]; 
	for(int i=1; i < dt_new.size(); i++)
		t_new[i]=t_new[i-1]+dt_new[i];

	/* time points of old temporal grid */
	vector<double> t(dt.size());
	t[0]=dt[0]; 
	for(int i=1; i < dt.size(); i++)
		t[i]=t[i-1]+dt[i];


	/* not fixed variables for old temporal grid */
	vector<pair<int,int>> I;
	for(int j=0; j < n; j++){
		for(int i=0; i < dt.size(); i++){
			if (find(fixIndex.begin(),fixIndex.end(),make_pair(j,i))==fixIndex.end())
				I.push_back(make_pair(j,i)); 
		}
	}


	/* new control vector, (un)fixed indices and ufix */
	ControlVector u_new(dt_new.size()); 
	ControlVector ufix(dt_new.size());
	bool null=true; // indicate whether fixed values are all zero
	vector<pair<int,int>> newfixIndex; 
	vector<pair<int,int>> I_new;
	for(int i=0; i < dt_new.size(); i++){
		// first time point of old temporal grid that is not less than t_new[i]
		auto start=lower_bound(t.begin(),t.end(), t_new[i],[](const double& a, double value){return a < value-1e-10;});
		int index=distance(t.begin(),start);
		// set value of new control vector on (t_new[i-1],t_new[i]]
		u_new[i]=u[index]; 	
		// determine new (un)fixed control variables on (t_new[i-1],t_new[i]] and ufix
		for(int j=0; j < n; j++){
			if(find(fixIndex.begin(),fixIndex.end(),make_pair(j,index))!=fixIndex.end()){ // u_j[index] fixed
				newfixIndex.push_back(make_pair(j,i)); 
				ufix[i][j]=u_new[i][j];
				if(ufix[i][j]!=0)
					null=false;
			}
			else
				I_new.push_back(make_pair(j,i));
		}
	}


	/* new box constraints variables */
	ControlVector w_new(I_new.size()); 
	ControlVector nu_new(I_new.size()); 
	for(int i=0; i < I_new.size(); i++){
		auto start=lower_bound(t.begin(),t.end(), t_new[I_new[i].second],[](const double& a, double value){return a < value-1e-10;});
		int index=distance(t.begin(),start);
		int pos=find(I.begin(),I.end(),make_pair(I_new[i].first,index))-I.begin();
		nu_new[i]=nu[pos];
		w_new[i]=w[pos];
	}


	/* new vectors G*S*(Sf-yd) and G*S*SGu_fix */
	// G*S*(Sf-yd)
	ControlVector newGstar_eta(dt_new.size());
	vector<vtype> newg1(dt_new.size()+1,vtype(gfs));
	vector<vtype> newg2(dt_new.size()+1,vtype(gfs));
	RHScontribution(gfs,dt_new,this->heat,this->YdProblem,newGstar_eta,newg1,newg2);
	// calculate G*S*SGu_fix 
	ControlVector newGstar_q(dt_new.size());
	if(!null){
		CGProblem cgproblem(ufix,dt_new); 
		RHScontribution(gfs,dt_new,cgproblem,newGstar_q);
	}

	/* new cutting plane matrix */ 
	if(A.N()>0){
		Matrix A_new=Matrix(A.N(), dt_new.size()*n, Matrix::row_wise);
		// set up sparsity pattern
		for(auto row=A_new.createbegin(); row!=A_new.createend(); ++row){
			for(int i=0; i < dt_new.size(); i++){
				// first time point of old temporal grid that is not less than t_new[i]
				int index=lower_bound(t.begin(), t.end(), t_new[i],[](const double& a, double value){return a < value-1e-10;})-t.begin();
				for(int j=0; j < n; j++){
					if(A.exists(row.index(),index*n+j))
						row.insert(i*n+j);
				}
			}
		}
	
		// set up values 
		for(int i=0; i < dt_new.size(); i++){
			// first time point of old temporal grid that is not less than t_new[i]
			int index=lower_bound(t.begin(), t.end(), t_new[i],[](const double& a, double value){return a < value-1e-10;})-t.begin();
			for(int j=0; j < n; j++){
				for(int l=0; l < A.N(); l++){
					if(A.exists(l,index*n+j))
						A_new[l][i*n+j]=A[l][index*n+j]*(dt_new[i]/dt[index]);
				}
			}
		}	
		A=A_new;
	}
		

	/* adapt all data */ 
	u=u_new; 
	w=w_new; 
	nu=nu_new; 
	fixIndex=newfixIndex;
	Sf=vector(dt_new.size()+1,vtype(gfs));
	Sf=newg1;
	Staryd=vector(dt_new.size()+1,vtype(gfs));
	Staryd=newg2;
	Gstar_eta=newGstar_eta; 
	Gstar_q=newGstar_q;	
}


/************************************************************************** 
 *
 * primal bound calculation 
 *
 *************************************************************************/
template<class GridView, class Constraints, class Problem, int k, int number>
void OCPMaster<GridView,Constraints,Problem,k,number>::PrimalBound(const GFS& gfs,const vector<double>& dt, 
const vector<vtype>& Sf, const vector<vtype>& Staryd, 
const ControlVector& v)
{
	/* calculate new primal bound */
	// y=SGv
	vector<vtype> y;
	CGProblem cgproblem(v,dt); 
	heatdriver(gfs,cgproblem,dt,y);

	vector<vtype> q=y;
	for(int i=0; i < dt.size(); i++) 
		q[i]+=Sf[i];
	

	vector<double> dt_new;
	// extend v
	ControlVector v_ex; 
	extendVector(v,dt,tgrid,v_ex,dt_new);
	// extend y
	vector<vtype> q_ex(dt_new.size()+1,vtype(gfs));
	extendVector<GFS>(q,dt,dt_new,q_ex);
	
	
	// objective value
	double bound=J(gfs,dt_new,v_ex,q_ex); 
	
	// error in finite element 
	vector<double> e(dt.size());
	FEError(gfs,dt,heat,y,Sf,Staryd,v,e);
	double error=accumulate(e.begin(),e.end(),0.0);
	bound+=error;

	// save solution appropiately
	vector<double> sol(dt.size()*n);
	for(int i=0; i < dt.size(); i++)
		for(int j=0; j < n; j++)
			sol[i*n+j]=v[i][j];

	/* update primal bound */
	setPrimalBound(bound,sol,dt);
	
}


/************************************************************************** 
 *
 * extacter objective value calculation 
 *
 *************************************************************************/
template<class GridView, class Constraints, class Problem, int k, int number>
double OCPMaster<GridView,Constraints,Problem,k,number>::objective(const GFS& gfs, const vector<double>& dt, 
const vector<vtype>& Sf, const ControlVector& u)
{
	/* calculate new primal bound */
	// y=SGu
	vector<vtype> y;
	CGProblem cgproblem(u,dt); 
	heatdriver(gfs,cgproblem,dt,y);
	for(int i=0; i < dt.size()+1; i++) 
		y[i]+=Sf[i];
	
	vector<double> dt_new;
	// extend u
	ControlVector u_ex; 
	extendVector(u,dt,tgrid,u_ex,dt_new);
	// extend y
	vector<vtype> y_ex(dt_new.size()+1,vtype(gfs));
	extendVector<GFS>(y,dt,dt_new,y_ex);
	
	// objective value
	double bound=J(gfs,dt_new,u_ex,y_ex); 

	return bound;
} 



/************************************************************************** 
 *
 * objective value calculation J(g,u,y)=1/2||y-yd||^2_{L^2(Q)}
 *					+\alpha/2||u-ud||^2_{L^2(0,T;R^n)}
 *
 *************************************************************************/
template<class GridView, class Constraints, class Problem, int k, int number>
double OCPMaster<GridView,Constraints,Problem,k,number>::J(const GFS& gfs,const vector<double>& dt,
const ControlVector& u,const vector<vtype>& y)
{
	// distribution 1/2 int_\Omega \int_[0,T] (y-y_d)^2 dt dx
	double bound; 
	L2Deviation(gfs,dt,y,YdProblem,bound); 

    	// distribution alpha/2 \int_[0,T] ||u-u_d||^2 dt 
    	for (int i=0; i < dt.size(); i++){
		double time=accumulate(dt.begin(),dt.begin()+i+1,0.0);
		auto diff=u[i]-heat.ud(time); 
       		bound+=diff.two_norm2()*alpha*dt[i]*0.5; 
	}

    	return bound;
}

/**************************************************************************
 *
 * outer approximation of subproblem
 *
 **************************************************************************/ 
template<class GridView, class Constraints, class Problem, int k, int number>
void OCPMaster<GridView,Constraints,Problem,k,number>::call_outer(Subproblem* sub,const GFS&gfs,const vector<double>& dt, 
const vector<pair<int,int>>& index,
const vector<vtype>& Sf,
const ControlVector& Gstar_eta,const ControlVector& Gstar_q,
Matrix& A, Vector&b, Vector & v, Vector& mu, Vector&w, Vector& nu, double& opt, 
ControlVector& u,const int subCount, OutputData& res)
{
	typedef ADMMsolver<MType> Solver;

	ControlVector ud(dt.size()); 
	ControlVector ua(dt.size()); 
	ControlVector ub(dt.size()); 
	for(int i=0; i < dt.size(); i++){
		double time = accumulate(dt.begin(),dt.begin()+i+1, 0.0); 
		ud[i]=heat.ud(time);
		ua[i]=heat.ua(time);
		ub[i]=heat.ub(time);
	}
	ControlVector u_=u;
	Vector v_=v;
	Vector mu_=mu;
	Vector w_=w; 
	Vector nu_=nu; 


	Solver solver=Solver(this->heat, this->YdProblem,gfs,dt,tgrid,alpha,rho,beta,eabs,erel,tol,reduction,iter_lin,index,Sf,Gstar_eta,Gstar_q,ud, ua, ub); 

	int iter =1; 
	int count=0; // frequency of low relative change

   	double time=0.0; // overall runtime in seconds
	double obj; // objective of previous iteration
	double bound=0; // current bound 
	bool stop_outer = false;   
		
	int iterations;
	
	cout << "\n---------------------- Outer Approximation Algorithm ----------------------\n" << endl;

	cout << "\n--- Outer Approximation Iteration " << iter << endl;
	// apply solver 
	if(A.N()==0){
		clock_t cstart=clock();
		solver.apply(u_,w_,nu_,bound,optSense,primalBoundFound,primalBound,optEps,iterations);
		time+=(double)(clock()-cstart)/CLOCKS_PER_SEC;      

		/* accept iteration */
		res.time_admm_solve.push_back((double)(clock()-cstart)/CLOCKS_PER_SEC);
		res.iter_admm.push_back(iterations);
		u=u_;
		nu=nu_;
		w=w_;
		// calculate objective value 
		opt=bound;
	}
	else{
		clock_t cstart=clock();	
		stop_outer=solver.apply(u_,v_,mu_,w_,nu_,A,b,bound,optSense,primalBoundFound,primalBound,optEps,iterations);
		time+=(double)(clock()-cstart)/CLOCKS_PER_SEC;      
			
		/* accept iteration */
		res.time_admm_solve.push_back((double)(clock()-cstart)/CLOCKS_PER_SEC);
		res.iter_admm.push_back(iterations);
		u=u_;
		w=w_; 	
		nu=nu_;
		// save active cutting planes
		vector<int> B; 
		for(int i=0; i < mu_.size(); i++){
			if(mu_[i]>1e-5) 
				B.push_back(i);
		}
		if(B.size()==mu_.size()){
			mu=mu_; 
			v=v_;
		}
		else if(B.size()==0){
			b=Vector();
			mu=Vector();
			v=Vector();
			A=Matrix();
		}
		else{
			b.resize(B.size());
			mu.resize(B.size());
			v.resize(B.size());
			for(int i=0; i < B.size(); i++){
				b[i]=b[B[i]]; 
				mu[i]=mu_[B[i]];
				v[i]=v_[B[i]];
			}

			Matrix A_B;
			A_B.setBuildMode(Matrix::row_wise);
			A_B.setSize(B.size(),A.M());
			for(auto row=A_B.createbegin(); row!=A_B.createend(); ++row)
				for( auto it=A[B[row.index()]].begin(); it!=A[B[row.index()]].end(); ++it) 
					row.insert(it.index()); 

			for(int i=0; i < B.size(); i++) 
				A_B[i]=A[B[i]]; 

			A=A_B;
		}

		opt=bound;

	} 
	res.dual_bound.push_back(bound);

    	while(!stop_outer && iter < iter_outer && time < timeouter && count<3){
		// store current objective 
		obj=opt; 

		// separation of control
		double* cut=new double[dt.size()*n]; 
		double rhs;

		clock_t ostart=clock();
		int status = separate(sub,dt,u,cut,rhs); 
		time+=(double)(clock()-ostart)/CLOCKS_PER_SEC;
		res.time_outer_separation.push_back((double)(clock()-ostart)/CLOCKS_PER_SEC);
        	if(status==0){ // control feasible: stop outer approximation
			cout<<"Control feasible for subproblem" << endl;
            		stop_outer = true; 
			sub->OAstopped(true);
        	}
       		else{ // control infeasible: add cutting plane and start next outer approximation iteration
	       		iter++;
			cout << "\n--- Outer Approximation Iteration " << iter << endl; 
 
			// append cutting plane to matrix
            		Matrix A_new=Matrix(A.N()+1,dt.size()*n,Matrix::row_wise); 
           		for(auto row=A_new.createbegin(); row!=A_new.createend(); ++row){
            			if (row.index()< A.N()){
			        	for( auto it=A[row.index()].begin(); it!=A[row.index()].end(); ++it) 
						row.insert(it.index());
                		}
               			else{
                    			for(int l=0; l<dt.size()*n; l++){
                        			if(cut[l]!=0)
                        			row.insert(l);  					
					}
            			}
			}
			
            		for(int i=0; i< A.N(); i++) 
                		A_new[i]=A[i];
           		int end = A.N();
            		for(int l=0; l<dt.size()*n; l++){
				if(cut[l]!=0)
                   			A_new[end][l]=cut[l]; 
			}
			// equilibration
			A_new[end]*=alpha;

			// extend right hand side
			Vector b_new(end+1);
			for(int i=0; i < end; i++) 
				b_new[i]=b[i]; 
			b_new[end]=rhs;
			b_new[end]*=alpha;
			
			// extend lagrange multiplicator 
			Vector mu_new(end+1);
			for(int i=0; i < end; i++)
            			mu_new[i]=mu[i];

			Vector v_new(end+1); 
			for(int i=0; i< end; i++)
				v_new[i]=v[i];
	
			// start outer approximation iteration
        	    	clock_t cstart=clock();
			stop_outer=solver.apply(u_,v_new,mu_new,w_,nu_,A_new,b_new,bound,optSense,primalBoundFound,primalBound,optEps,iterations);
			time+=(double)(clock()-cstart)/CLOCKS_PER_SEC;

			/* accept iteration */ 
			res.time_admm_solve.push_back((double)(clock()-cstart)/CLOCKS_PER_SEC);
			res.iter_admm.push_back(iterations);
			res.dual_bound.push_back(bound);
			u=u_;
			w=w_; 
			nu=nu_;
			// save active cutting planes
			vector<int> B; 
			for(int i=0; i < mu_new.size(); i++){
				if(mu_new[i]>1e-5) 
					B.push_back(i);
			}		
			if(B.size()==mu_new.size()){
				A=A_new;
				b=b_new; 
				mu=mu_new; 
				v=v_new;
			}
			else if(B.size()==0){
				b=Vector();
				mu=Vector();
				v=Vector();
				A=Matrix();
			}
			else{
				b.resize(B.size());
				mu.resize(B.size());
				v.resize(B.size());
				for(int i=0; i < B.size(); i++){
					b[i]=b_new[B[i]]; 
					mu[i]=mu_new[B[i]];
					v[i]=v_new[B[i]];
				}

				Matrix A_B;
				A_B.setBuildMode(Matrix::row_wise);
				A_B.setSize(B.size(),A_new.M());
				for(auto row=A_B.createbegin(); row!=A_B.createend(); ++row)
					for( auto it=A_new[B[row.index()]].begin(); it!=A_new[B[row.index()]].end(); ++it) 
						row.insert(it.index()); 

				for(int i=0; i < B.size(); i++) 
					A_B[i]=A_new[B[i]]; 

				A=A_B;
			}
		
			opt=bound;
			
			// check relative change
			if((opt-obj)/obj<rchange)
				count++; 
			else
				count=0;
		}		
		delete[] cut;	
	}
	res.iter_cut+=iter;
    	cout << "\n---------------------- End Outer Approximation Algorithm ----------------------\n" << endl;
    	
}

/*********************************************************************
 *
 * get bounds for subproblem 
 *
 *********************************************************************/
template<class GridView, class Constraints, class Problem, int k, int number>
int OCPMaster<GridView,Constraints,Problem,k,number>::getBounds(Subproblem* sub,int count)
{
	// write subproblem info to file
	((SubType*)sub)->writeInfo(count); 
	// write current primal info to file 
	writePrimalInfo(count);

	OutputData result;
	vector<double> dt= ((SubType*)sub)->getdeltaT(); // step wides
	GFS gfs=((SubType*)sub)->getGFS();
	
	/* check feasibility of subproblem */
	clock_t cstart=clock();
	sub->setFeasible(true);	
	if(sub->getNumFixPoints()>1){
		int t=sub->getNumFixPoints();
		double* cut=new double[t]; 
		double rhs;
		// prove feasibility
		int status = this->separate(sub,vector<double>(),ControlVector(),cut,rhs); 

		if(status==1) // control given by fixations infeasible
			sub->setFeasible(false);
		delete[] cut;
	}
	result.time_feasibility=(double)(clock()-cstart)/CLOCKS_PER_SEC; 
	

	/* calculation of dual bound */
	if(sub->isFeasible()){ 
		// get initial control
		ControlVector u=((SubType*)sub)->getControl();
		// get initial fixed indices
		vector<pair<int,int>> fixIndex=((SubType*)sub)->getFixIndices();

		// update fixed control variables
		bool changed=updateFixControls(sub,dt,fixIndex,u);
		if(changed){ // update data of subproblem
			((SubType*)sub)->setControl(u);
			((SubType*)sub)->setFixIndices(fixIndex);
			((SubType*)sub)->resetBox();
			cstart=clock();
			updateFixData(sub,fixIndex,gfs,dt,u);
			result.time_updateFixData=(double)(clock()-cstart)/CLOCKS_PER_SEC;
		}
		else{
			// separate last control (if possible)
			double* cut=new double[dt.size()*n]; 
			double rhs;
			cstart=clock();
			int status = separate(sub,dt,u,cut,rhs); 
			result.time_separation=(double)(clock()-cstart)/CLOCKS_PER_SEC;
			// outer approximation: if last control not feasible for subproblem add new cutting plane
			if(status!=0){	
				// add equilibrated cutting plane to problem
				for(int i=0; i < dt.size()*n; i++)
					cut[i]*=alpha;
				rhs*=alpha;
				((SubType*)sub)->extendCut(cut,rhs);
			}
			delete[] cut;
			bool grid_changed=((SubType*)sub)->getGridChanged();
			// if last control still optimal for subproblem resort to branching
			if(status==0 && count >0 && !(grid_changed)){
				// set dual solution  to last control
				vector<double> sol(dt.size()*n);
				for(int i=0; i < dt.size(); i++)
					for(int j=0; j < n; j++)
						sol[i*n+j]=u[i][j];
				sub->setDualSolution(sol);

				result.writeInfo(count,sub->getNumFixPoints(), fixIndex.size(),*min_element(dt.begin(),dt.end()),sub->getDualBound(),false);

				/* write solution to file */
				stringstream filename;
				filename << BNB_OUTPUT_PATH << OUTPUTDIR << "sub" << count <<"/DualSol.txt";
				storeMatrixMarket<ControlVector>(u,filename.str());
	
				return false;
			}
		}
		
		// get data for outer approximation
		ControlVector Gstar_q,Gstar_eta;
		vector<vtype> Sf;
		Matrix A; 
		Vector v, mu,w, nu, b;
		((SubType*)sub)->getLambda(fixIndex,Sf,Gstar_eta,Gstar_q, A, b,v, mu, w, nu);

	
		/* outer approximation as long as dual bound >= primal bound and grid refinement possible */
		bool stop=false;
		double opt;
		double error; 
		int refined;
		int steps=0; // number of grid refinements
		do{
			// reset/initialize some data
			opt=0;
			refined=0;
			result.status=-1;
			sub->OAstopped(false);

			// call outer approximation if not everything is fixed
			if(fixIndex.size()!=dt.size()){
				cstart=clock();
				call_outer(sub,gfs,dt,fixIndex,Sf,Gstar_eta,Gstar_q,A,b,v,mu,w,nu,opt,u,count,result);
				result.time_solve+=(double)(clock()-cstart)/CLOCKS_PER_SEC;
			} 
			else{
				sub->OAstopped(true);
				opt=objective(gfs,dt,Sf,u); 
			}

			// set dual bound, solution 
			vector<double> sol(dt.size()*n);
			for(int i=0; i < dt.size(); i++)
				for(int j=0; j < n; j++)
					sol[i*n+j]=u[i][j];
			sub->setDualBound(opt);
			sub->setDualSolution(sol);

			// update data of subproblem
			((SubType*)sub)->setData(A, b,v, mu, w, nu, u);

			//check if the calculated solution is feasible and integral
			if(sub->doesOAstopped()){
				result.status=0;
      				double eps = 1e-5;
      				bool Integral = true;
      				for(int i = 0; i < sol.size(); i++){
					double value = sol[i];
	 				if(value - floor(value) > eps && value - floor(value) < 1.-eps){
	    					Integral = false;
						break;
					}
				}
				if(Integral){
					vector<vtype> Staryd=((SubType*)sub)->getStaryd();
					PrimalBound(gfs,dt,Sf,Staryd,u);
					result.status=1;
				}
			}

			// grid refinement 
			vector<double> dt_old=dt;
			if(fathomSubproblem(sub)){
				vector<vtype> Staryd=((SubType*)sub)->getStaryd();
				vector<double> e(dt.size());

				// unscaled dual variables
				BlockVector<FieldVector<double,1>> mu_=mu; 
				mu_*=rho;
				/* a posteriori error  */
				// calculation of y=Su
				vector<vtype> y;
				CGProblem cgproblem(u,dt);
				heatdriver(gfs,cgproblem,dt,y);
				Error(gfs,dt,Sf,Staryd,fixIndex,A,mu_,u,e);
				error=accumulate(e.begin(),e.end(),0.0);
				double obj=objective(gfs,dt,Sf,u); 
				obj+=error;
				sub->setDualBound(obj);
				if(!fathomSubproblem(sub))
					refined=refine(sub, gfs, dt, e);
			}
			if(refined==1 && steps<3){
				cout << "\n=== Refinement of grid" << endl;
				vector<vtype> Staryd(dt.size()+1,vtype(gfs));
				updateData(fixIndex,gfs,dt,dt_old,Sf,Staryd,Gstar_eta,Gstar_q,A,w,nu,u);
				((SubType*)sub)->setData(fixIndex, gfs, dt, Sf,Staryd, Gstar_eta,Gstar_q,A, w, nu, u);
				((SubType*)sub)->GridChanged();
				steps++;
			}
			else{	
				stop=true;
				dt=dt_old;
			}
				
		} while(!stop);

		/* write statistics to file */
		bool fathom=fathomSubproblem(sub); 
		if(refined==1 && steps==3 || refined==-1) // grid not refined 
			fathom=false;
		result.writeInfo(count,sub->getNumFixPoints(), fixIndex.size(),*min_element(dt.begin(),dt.end()),sub->getDualBound(),fathom);

		/* write solution to file */
		stringstream filename;
		filename << BNB_OUTPUT_PATH << OUTPUTDIR << "sub" << count <<"/DualSol.txt";
		storeMatrixMarket<ControlVector>(u,filename.str());

		/* heuristic solution with new control */ 
		ControlVector hsol(dt.size());
		heuristic(sub,gfs,dt,u,hsol);
		vector<vtype> Staryd=((SubType*)sub)->getStaryd();
		// update primal bound
		PrimalBound(gfs,dt,Sf,Staryd,hsol);

		return fathom;

	}
	return true;
}



/*********************************************************************
 *
 * write methods 
 *
 *********************************************************************/
template<class GridView, class Constraints, class Problem, int k, int number>
void OCPMaster<GridView,Constraints,Problem,k,number>::writePrimalInfo(int count)
{
	if(this->primalBoundFound){
		stringstream dir;
		dir << BNB_OUTPUT_PATH << OUTPUTDIR << "sub" << count; 
		string filename=dir.str() + "/primalSol.txt";
		fstream file; 
		file.precision(std::numeric_limits<double>::digits10+2);
		file.open(filename, ios::out | ios::trunc);
		// write primal bound and solution to file 
		file << this->primalBound << endl;
		for(int i=0; i < this->primalSolution.size(); i++)
			file << this->primalSolution[i] << endl;
		file.close();

		filename=dir.str() + "/primalTimeGrid.txt";
		file.open(filename, ios::out | ios::trunc);
		// write primal solution to file 
		for(int i=0; i < primalTimeGrid.size(); i++)
			file << primalTimeGrid[i] << endl;
		file.close();
	}
}

template<class GridView, class Constraints, class Problem, int k, int number>
void OCPMaster<GridView,Constraints,Problem,k,number>::writeParameterSetting()
{
	stringstream filename;
	filename << BNB_OUTPUT_PATH << OUTPUTDIR << "settings.txt";
	fstream stream; 
	stream.precision(std::numeric_limits<double>::digits10+2);
	stream.open(filename.str(), ios::out | ios::app);

	stream << "switches=" << this->n << endl;
	stream << "spatial_degree=" << this->degree << endl;
	stream << "Tikhonov=" << this->alpha << endl; 
	stream << "Rho=" << this-> rho << endl; 
	stream << "Beta=" << this->beta << endl;
	stream << "Iter_Lin=" << this->iter_lin << endl;
	stream << "Reduction=" << this->reduction << endl;
	stream << "ADMM_AbsTol=" << this->eabs << endl;
	stream << "ADMM_RelTol=" << this->erel << endl;
	stream << "ADMM_ETol=" << this->tol << endl;
	stream << "ObjReduction=" <<this ->rchange << endl;
	stream << "Iter_Outer=" << this->iter_outer << endl;
	stream << "TimeLimit=" << this->timeouter << endl;
	stream << "OPTEps=" << this->optEps << endl;
	stream << "rfactor=" << rfactor << endl; 
	stream << "[grid_yd]" << endl;
	for(int i=0; i < this->tgrid.size(); i++)
		stream << "dt_" << i << "="<< this->tgrid[i] << endl;

	stream.close();
}


/*********************************************************************
 *
 * refinement stratedy based on the DWR method
 *
 *********************************************************************/
template<class GridView, class Constraints, class Problem, int k, int number>
int OCPMaster<GridView,Constraints,Problem,k,number>::refine(Subproblem* sub, GFS& gfs, vector<double>& dt, vector<double>& e)
{
	for(int i=0; i < dt.size(); i++) 
		e[i]=abs(e[i]);
	// sort cellwise error in descending order
	vector<int> order(dt.size()); 
	for(int i=0; i < dt.size(); i++) 
		order[i]=i;
	sort(order.begin(),order.end(),[&](int left, int right){return e[left]>e[right];});
	 
	// determine intervals to refine
	double error=accumulate(e.begin(),e.end(),0.0);
	int pos=0; 
	double sum=0.0;
	while(sum<rfactor*error-1e-10){
		sum+=e[order[pos]];
		pos++;
	}

	// sort cells to refine in ascending order
	vector<int> cells(pos); 
	for(int i=0; i < pos; i++) 
		cells[i]=order[i];
	sort(cells.begin(),cells.end());

	// determine new temporal grid 
	vector<double> dt_new=dt; 
	for(int i=0; i < pos; i++){
		dt_new[cells[i]+i]=0.5*dt[cells[i]];
		dt_new.emplace(dt_new.begin()+cells[i]+i+1,0.5*dt[cells[i]]);
	}	
	
	if(dt_new.size()==dt.size()) 
		return -1; 
	else{
		dt=dt_new;
		return 1;
	}

}


/*********************************************************************
 *
 * cellwise error in optimal control
 *
 *********************************************************************/
template<class GridView, class Constraints, class Problem, int k, int number>
void OCPMaster<GridView,Constraints,Problem,k,number>::Error(const GFS& gfs,
const vector<double>&dt, const vector<vtype>& g1,const vector<vtype>& g2, 
const vector<pair<int,int>>& fixIndex, 
const Matrix& A, const Vector& mu, 
const ControlVector& u, vector<double>& e)
{
	vector<vtype> y;
	CGProblem cgproblem(u,dt);
	heatdriver(gfs,cgproblem,dt,y);
	OCPError(gfs,dt,this->heat,this->YdProblem,this->alpha,y,g1,g2,fixIndex,A,mu,u,e);
}
		
