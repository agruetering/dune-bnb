// Copyright © Alexandra Grütering, see file LICENSE.md in module root
// License: GLP-3.0 or later
#ifndef ERRORFUNCTIONS_HH
#define ERRORFUNCTIONS_HH


/**************************************************************************
 * cellwise error
 * 		(y(t_i)-yd(t_i),y(t_i))_{L^2(\Omega)) for all time points
 *
 ***************************************************************************/
template<class GFS>
void cellError(const GFS& gfs,
const AdjointProblemInterface<typename GFS::Traits::GridView>& YdProblem, 
const vector<typename PDELab::Backend::Vector<GFS,double>>& y, 
vector<double>& e) 
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

	LFS lfs(gfs);
  	LFSCache lfs_cache(lfs);
  	BasisCache cache;

	// reset cellwise error 
	e=vector<double>(y.size(),0); 
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
			double t=0.0;
           		for (int i=0; i < y.size(); i++){
            			// evaluate desired temperature
				double des=YdProblem.rhs(eg,ip.position(),t);
               			RangeType val=RangeType(0); 
                		// evaluate y(t_i,x)
                		for (size_t j=0; j<lfs.size(); j++){
				     val+=y[i][lfs_cache.containerIndex(j)]*phihat[j];
               			}
				
                		double term=(val-des)*val;
                		term*=factor; 
				e[i]+=term;
				t+=YdProblem.deltaT(i);
            		}
        	}
	}  
}


/**************************************************************************
 * cellwise error
 * 		(y(t_i),p(t_i))_{L^2(\Omega)) for all time points
 *
 ***************************************************************************/
template<class GFS>
void cellError(const GFS& gfs,
const vector<typename PDELab::Backend::Vector<GFS,double>>& y, 
const vector<typename PDELab::Backend::Vector<GFS,double>>& p, 
vector<double>& e)
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

	LFS lfs(gfs);
  	LFSCache lfs_cache(lfs);
  	BasisCache cache;

	// reset cellwise error 
	e=vector<double>(y.size(),0); 
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
           		for (int i=0; i < y.size(); i++){
               			RangeType val0=RangeType(0); 
				RangeType val1=RangeType(0); 
                		// evaluate y(t_i,x) and p(t_i,x)
                		for (size_t j=0; j<lfs.size(); j++){
					val0+=y[i][lfs_cache.containerIndex(j)]*phihat[j];
					val1+=p[i][lfs_cache.containerIndex(j)]*phihat[j];
                		}
				
                		double term=val0*val1;
                		term*=factor; 
				e[i]+=term;
            		}
        	}
	}  
}


/**************************************************************************
 * cellwise error
 * 	(\nabla y(t_i),\nabla p(t_i))_{L^2(\Omega)) for all time points
 *
 ***************************************************************************/
template<class GFS>
void cellErrorLaplace(const GFS& gfs,
const vector<typename PDELab::Backend::Vector<GFS,double>>& y, 
const vector<typename PDELab::Backend::Vector<GFS,double>>& p, 
vector<double>& e)
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

	const int dim=GV::dimension;
	LFS lfs(gfs);
  	LFSCache lfs_cache(lfs);
  	BasisCache cache;

	// reset cellwise error 
	e=vector<double>(y.size(),0);  
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

        		// evaluate gradient of shape functions
        		auto& gradphihat = cache.evaluateJacobian(ip.position(),
                		lfs.finiteElement().localBasis());
	
        		// transform gradients from reference element to grid element
        		const auto S = geo.jacobianInverseTransposed(ip.position());
        		auto gradphi = makeJacobianContainer(lfs);
        		for (size_t i=0; i<lfs.size(); i++)
          			S.mv(gradphihat[i][0],gradphi[i][0]);

            		// loop over all time points
           		for (int i=0; i < y.size(); i++){
              			// evaluate \nabla y(t_i) and \nabla p(t_i)
               			FieldVector<RangeType,dim> grady(0);
				FieldVector<RangeType,dim> gradp(0);
                		for (size_t j=0; j<lfs.size(); j++){
				     	grady.axpy(y[i][lfs_cache.containerIndex(j)],gradphi[j][0]);
					gradp.axpy(p[i][lfs_cache.containerIndex(j)],gradphi[j][0]);
                		}
				
                		double term=0.0; 
				for(int k=0; k < dim; k++)
					term+=grady[k]*gradp[k];
                		term*=factor; 
				e[i]+=term;
            		}
        	}
	}  
}


/**************************************************************************
 * cellwise error
 * 	(\sum u_j(t_i)\Psi_j,p(t_i))_{L^2(\Omega)) for all time points
 *
 ***************************************************************************/
template<class GFS, class ControlVector>
void cellError(const GFS& gfs, 
const HeatProblemInterface<typename GFS::Traits::GridView,ControlVector::block_type::dimension>& heat, 
const ControlVector& u, 
const vector<typename PDELab::Backend::Vector<GFS,double>>& p, 
vector<double>& e)
{

	const int n=ControlVector::block_type::dimension;

	typedef typename GFS::Traits::GridView GV; // spatial grid view
    	typedef typename GFS::Traits::FiniteElementMap FEM; // finite element map
    	typedef typename PDELab::Backend::Vector<GFS,double> Y; // type of Dof vector
   	typedef typename PDELab::DiscreteGridFunction<GFS,Y> DGF;  // discrete grid function

    	typedef typename DGF::Traits::RangeType RangeType; // range type discrete grid function
   	typedef PDELab::LocalFunctionSpace<GFS> LFS; // local function space 
  	typedef PDELab::LFSIndexCache<LFS> LFSCache; // local function space cache
  	typedef typename GFS::Traits::FiniteElementType::Traits::LocalBasisType LocalBasis; // local basis 
  	typedef typename PDELab::LocalBasisCache<LocalBasis> BasisCache; // local basis cache

	LFS lfs(gfs);
  	LFSCache lfs_cache(lfs);
  	BasisCache cache;

	// reset cellwise error 
	e=vector<double>(u.size(),0); 
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

			// evaluate Psi
        		auto psi=heat.Psi(eg,ip.position());

            		// loop over all time points
           		for (int i=0; i < u.size(); i++){
               			RangeType val0=RangeType(0);
                		// evaluate \sum u_j(t_i)\Psi_j(x) and p(t_i,x)
                		for (size_t j=0; j<lfs.size(); j++){
				    val0+=p[i][lfs_cache.containerIndex(j)]*phihat[j];
                		}
				double val1=0.0; 
				for(int j=0; j < n; j++) 
					val1+=u[i][j]*psi[j];
				
                		double term=val0*val1;
                		term*=factor; 
				e[i]+=term;
            		}
        	}
	}  
}


/**************************************************************************
 * cellwise error
 * 	(f(t_i),p(t_i))_{L^2(\Omega)) for all time points
 *
 ***************************************************************************/
template<class GFS>
void cellError(const GFS& gfs,
const vector<double>& dt, 
ProblemInterface<typename GFS::Traits::GridView>& heat, 
const vector<typename PDELab::Backend::Vector<GFS,double>>& p, 
vector<double>& e)
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

	LFS lfs(gfs);
  	LFSCache lfs_cache(lfs);
  	BasisCache cache;

	// reset cellwise error 
	e=vector<double>(dt.size()+1,0); 
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
           		for (int i=0; i < p.size(); i++){
               			RangeType val0=RangeType(0);
               	 		// evaluate f(t_i,x) and p(t_i,x)
                		for (size_t j=0; j<lfs.size(); j++){
					val0+=p[i][lfs_cache.containerIndex(j)]*phihat[j];
                		}
				double time=accumulate(dt.begin(),dt.begin()+i,0.0); 
				heat.setTime(time); 
				double val1=heat.q(eg,ip.position());
				
                		double term=val0*val1;
                		term*=factor; 
				e[i]+=term;
            		}
        	}
	}  
}



/*****************************************************************************************
 * cellwise error
 * 	L_u'(X)(\tilde u)+ error terms of box constraints for all time points
 *
 *****************************************************************************************/
template<class GFS, class ControlVector>
void cellErrorApprox(const GFS& gfs,
const vector<double>& dt,
const vector<pair<int,int>>& I, // not fixed variables 
const HeatProblemInterface<typename GFS::Traits::GridView,ControlVector::block_type::dimension>& heat, 
const double& alpha,
const BCRSMatrix<FieldMatrix<double,1,1>>& A,  
const BlockVector<FieldVector<double,1>>& mu, 
const ControlVector& u,
const vector<typename PDELab::Backend::Vector<GFS,double>>& p, 
vector<double>& e)
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

	using Vector=BlockVector<FieldVector<double,1>>; 
	

	// G*p(t_i), A^Tmu(t_i) for all time points
	ControlVector Gp(dt.size()); 
	Gstar(gfs,heat,p,Gp);
	const int n=ControlVector::block_type::dimension;
	Vector h(dt.size()*n); 
	A.mtv(mu,h);
	 

	e=vector<double>(dt.size(),0); 
	// evaluation of L_u'(X)(\tilde u)+ error terms of box constraints
	for(int i=0; i < dt.size(); i++){
		double oldtime=accumulate(dt.begin(),dt.begin()+i,0.0);
		double time=accumulate(dt.begin(),dt.begin()+i+1,0.0);
		FieldVector<double,n> h0; // -1/alpha(G*I_k^1p(t_i)+A^Tmu(t_{i+1}))+ud(t_i)
		FieldVector<double,n> h1; // -1/alpha(G*I_k^1p(t_{i+1})+A^Tmu(t_{i+1}))+ud(t_{i+1})
		if(i==0)
			h0=0; 
		else{
			h0.axpy(-1.0/alpha,Gp[i-1]); 
			h0+=heat.ud(oldtime); 
			for(int j=0; j< n; j++)
				h0[j]-=1.0/(alpha*dt[i])*h[i*n+j];
		}
		h1.axpy(-1.0/alpha,Gp[i]); 
		for(int j=0; j< n; j++)
			h1[j]-=1.0/(alpha*dt[i])*h[i*n+j];
		h1+=heat.ud(time); 
		for(int j=0; j<n; j++){
			if(find(I.begin(),I.end(),make_pair(j,i))==I.end()) // (j,i) fixed
				continue;

			double ta; // time point when -1/alpha(G*I_k^1p(t)+A^Tmu(t))+ud(t) becomes ua
			double tb;// time point when -1/alpha(G*I_k^1p(t)+A^Tmu(t))+ud(t) becomes ub
			double value=Gp[i][j]+(1.0/dt[i])*h[i*n+j]+alpha*u[i][j]-alpha*heat.ud(time)[j];
			if(h0[j]==h1[j])
				e[i]+=max(heat.ua(time)[j],min(h1[j],heat.ub(time)[j]))*dt[i]*value;
			else{
				ta=((heat.ua(time)[j]-h0[j])*dt[i])/(h1[j]-h0[j])+oldtime; 
				tb=((-h0[j]+heat.ub(time)[j])*dt[i])/(h1[j]-h0[j])+oldtime;
				if(h0[j]<h1[j] && ta<=time-1e-10 && tb>=oldtime+1e-10){ // monoton increasing and projection not constant
					e[i]+=(min(time,tb)-max(oldtime,ta))*0.5*(max(h0[j],heat.ua(time)[j])+min(h1[j],heat.ub(time)[j]))*value; // monoton increasing part
					e[i]+=max(time-tb,0.0)*heat.ub(time)[j]*value;  // constant ub part
					e[i]+=max(ta-oldtime,0.0)*heat.ua(time)[j]*value;  // constant ua part
					e[i]+=max(time-tb,0.0)*0.5*(u[i][j]-heat.ub(time)[j])*alpha*(h1[j]-heat.ub(time)[j]); // contribution \tilde nu_+(u-ub)
					e[i]+=max(ta-oldtime,0.0)*0.5*(heat.ua(time)[j]-u[i][j])*alpha*(heat.ua(time)[j]-h0[j]); // contribution \tilde nu_-(ua-u)
				}
				else if(h0[j]<h1[j] && tb<oldtime+1e-10){ // monoton increasing and projection constant ub
					e[i]+=heat.ub(time)[j]*value*dt[i]; 
					e[i]+=dt[i]*(u[i][j]-heat.ub(time)[j])*alpha*(0.5*(h0[j]+h1[j])-heat.ub(time)[j]); // contribution \tilde nu_+(u-ub)
				}
				else if(h0[j]<h1[j] && ta>time-1e-10){ // monoton increasing and projection constant ua
					e[i]+=heat.ua(time)[j]*value*dt[i]; 
					e[i]+=dt[i]*(heat.ua(time)[j]-u[i][j])*alpha*(heat.ua(time)[j]-0.5*(h0[j]+h1[j])); // contribution \tilde nu_-(ua-u)
				}
				else if(h0[j]>h1[j] && tb<=time-1e-10 && ta>=oldtime+1e-10){ // mono/home/user/agrueter/Promotion/Projekte/DUNE/dune-bnb/output/1DTest2/sub112/b.txtton decreasing and projection not constant
					e[i]+=(min(time,ta)-max(oldtime,tb))*0.5*(min(h0[j],heat.ub(time)[j])+max(h1[j],heat.ua(time)[j]))*value; // monoton decreasing part
					e[i]+=max(time-ta,0.0)*heat.ua(time)[j]*value;  // constant ua part
					e[i]+=max(tb-oldtime,0.0)*heat.ub(time)[j]*value;  // constant ub part
					e[i]+=max(tb-oldtime,0.0)*0.5*(u[i][j]-heat.ub(time)[j])*alpha*(h0[j]-heat.ub(time)[j]); // contribution \tilde nu_+(u-ub)
					e[i]+=max(time-ta,0.0)*0.5*(heat.ua(time)[j]-u[i][j])*alpha*(heat.ua(time)[j]-h1[j]); // contribution \tilde nu_-(ua-u)
				}
				else if(h0[j]>h1[j] && tb>time-1e-10){ // monoton decreasing and projection constant ub
					e[i]+=heat.ub(time)[j]*value*dt[i];
					e[i]+=dt[i]*(u[i][j]-heat.ub(time)[j])*alpha*(0.5*(h0[j]+h1[j])-heat.ub(time)[j]); // contribution \tilde nu_+(u-ub)
				}
				else if(h0[j]>h1[j] && ta<oldtime+1e-10){ // monoton decreasing and projection constant ua
					e[i]+=heat.ua(time)[j]*value*dt[i]; 
					e[i]+=dt[i]*(heat.ua(time)[j]-u[i][j])*alpha*(heat.ua(time)[j]-0.5*(h0[j]+h1[j])); // contribution \tilde nu_-(ua-u)
				}
			}
		}	
	}
}


/**************************************************************************
 * cellwise error
 * 	L_u'(X)(u) for all time points
 *
 ***************************************************************************/
template<class GFS, class ControlVector>
void cellError(const GFS& gfs,
const vector<double>& dt, 
const vector<pair<int,int>>& I, // not fixed variables
const HeatProblemInterface<typename GFS::Traits::GridView,ControlVector::block_type::dimension>& heat, 
const double& alpha,
const BCRSMatrix<FieldMatrix<double,1,1>>& A, 
const BlockVector<FieldVector<double,1>>& mu, 
const ControlVector& u, 
const vector<typename PDELab::Backend::Vector<GFS,double>>& p, 
vector<double>& e)
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

	using Vector=BlockVector<FieldVector<double,1>>; 
	

	// G*p(t_i), A^Tmu(t_i) for all time points
	ControlVector Gp(dt.size()); 
	Gstar(gfs,heat,p,Gp);
	const int n=ControlVector::block_type::dimension;
	Vector h(dt.size()*n); 
	A.mtv(mu,h);
	 

	e=vector<double>(dt.size(),0); 
	// evaluation of L_u'(X)(u)
	for(int k=0; k < I.size(); k++){
		int i=I[k].second; 
		int j=I[k].first;
		double time=accumulate(dt.begin(),dt.begin()+i+1,0.0); 
		double value=Gp[i][j]+alpha*u[i][j]-alpha*heat.ud(time)[j];
		value*=dt[i]; 
		value+=h[i*n+j]; 
		e[i]+=value*u[i][j];
	}		
}


/*********************************************************************
 *
 * cellwise error in finite element method
 *
 *********************************************************************/
template<class GFS, class ControlVector>
void FEError(const GFS& gfs, // spatial grid function space 
const vector<double>&dt, // temporal grid
HeatProblemInterface<typename GFS::Traits::GridView,ControlVector::block_type::dimension>& heat, // heat problem with form functions
vector<typename PDELab::Backend::Vector<GFS,double>>& y, // y=Su 
const vector<typename PDELab::Backend::Vector<GFS,double>>& g1, // Sf  
const vector<typename PDELab::Backend::Vector<GFS,double>>& g2, // S*yd
const ControlVector& u, // control u
vector<double>& e) // cellwise error vector
{

	typedef typename PDELab::Backend::Vector<GFS,double> vtype;

	// y=S(Gu+f)
	for(int i=0; i < y.size(); i++) 
		y[i]+=g1[i];
	// p=S*(S(Gu+f)-yd)
	vector<vtype> p;
	AdjointProblem<GFS> adjointProblem(dt,y,gfs);
	adjointdriver(gfs,adjointProblem,p);
	for(int i=0; i < p.size(); i++) 
		p[i]-=g2[i];

	vector<double> r; // auxiliary vector to save cellwise errors 
	e=vector<double>(dt.size(),0); 

	/* cellwise error of heat equation */
	vector<vtype> h1(dt.size(),vtype(gfs)); 
	vector<vtype> h2(dt.size(),vtype(gfs)); 
	for(int i=0; i < dt.size(); i++){
		h1[i]=y[i+1]; 
		h1[i]-=y[i];
		h2[i]=p[i+1];
		h2[i]-=p[i];
	}
	cellError(gfs,h1,h2,r); 
	for(int i=0; i < dt.size(); i++) 
		e[i]+=r[i]; 
	for(int i=0; i < dt.size(); i++) 
		h1[i]=y[i+1]; 
	cellErrorLaplace(gfs,h1,h2,r); 
	for(int i=0; i < dt.size(); i++) 
		e[i]+=0.5*r[i]*dt[i];
	cellError(gfs,heat,u,h2,r);
	for(int i=0; i < dt.size(); i++) 
		e[i]-=0.5*dt[i]*r[i]; 
	cellError(gfs,dt,heat,p,r); 
	for(int i=0; i < dt.size(); i++)
		e[i]+=0.5*dt[i]*(r[i]-r[i+1]);
}


/*********************************************************************
 *
 * cellwise error in optimal control
 *
 *********************************************************************/
template<class GFS, class ControlVector>
void OCPError(const GFS& gfs, // spatial grid
const vector<double>&dt, // temporal grid
HeatProblemInterface<typename GFS::Traits::GridView,ControlVector::block_type::dimension>& heat, // heat problem with form functions
AdjointProblemInterface<typename GFS::Traits::GridView>& YdProblem, // problem with desired state yd
double alpha, // Tikhonov term
vector<typename PDELab::Backend::Vector<GFS,double>>& y, // y=Su 
const vector<typename PDELab::Backend::Vector<GFS,double>>& g1, // Sf
const vector<typename PDELab::Backend::Vector<GFS,double>>& g2, // S*yd
const vector<pair<int,int>>& fixIndex, // fixed Indizes
const BCRSMatrix<FieldMatrix<double,1,1>>& A, // cutting plane matrix
const BlockVector<FieldVector<double,1>>& mu, // Lagrange multiplicator of cutting planes
const ControlVector& u,  // control u
vector<double>& e) // cellwise error vector ‚
{

	typedef typename PDELab::Backend::Vector<GFS,double> vtype;
	const int n=ControlVector::block_type::dimension;

	// y=S(Gu+f)
	for(int i=0; i < y.size(); i++) 
		y[i]+=g1[i];
	// p=S*(S(Gu+f)-yd)
	vector<vtype> p;
	AdjointProblem<GFS> adjointProblem(dt,y,gfs);
	adjointdriver(gfs,adjointProblem,p);
	for(int i=0; i < p.size(); i++) 
		p[i]-=g2[i];

	vector<double> r; // auxiliary vector to save cellwise errors 
	e=vector<double>(dt.size(),0); 

	/* cellwise error of adjoint equation */ 
	YdProblem.setdeltaT(dt);
	cellError(gfs,YdProblem,y,r); 
	for(int i=0; i < dt.size(); i++) 
		e[i]+=0.5*dt[i]*(r[i]-r[i+1]);
	vector<vtype> h1(dt.size(),vtype(gfs)); 
	vector<vtype> h2(dt.size(),vtype(gfs)); 
	for(int i=0; i < dt.size(); i++){
		h1[i]=y[i+1]; 
		h1[i]-=y[i];
		h2[i]=p[i+1];
	}
	cellErrorLaplace(gfs,h1,h2,r);
	for(int i=0; i < dt.size(); i++) 
		e[i]+=0.5*dt[i]*r[i];


	/* cellwise error of heat equation */ 
	for(int i=0; i < dt.size(); i++){
		h2[i]=p[i+1];
		h2[i]-=p[i];
	}
	cellError(gfs,h1,h2,r); 
	for(int i=0; i < dt.size(); i++) 
		e[i]+=r[i]; 
	for(int i=0; i < dt.size(); i++) 
		h1[i]=y[i+1]; 
	cellErrorLaplace(gfs,h1,h2,r); 
	for(int i=0; i < dt.size(); i++) 
		e[i]+=0.5*r[i]*dt[i];
	cellError(gfs,heat,u,h2,r);
	for(int i=0; i < dt.size(); i++) 
		e[i]-=0.5*dt[i]*r[i]; 
	cellError(gfs,dt,heat,p,r); 
	for(int i=0; i < dt.size(); i++)
		e[i]+=0.5*dt[i]*(r[i]-r[i+1]);

	/* cellwise error of control and box constraints */
	vector<pair<int,int>> I; 
	for(int i=0; i < dt.size(); i++){ 
		for(int j=0; j < n; j++){
			if (find(fixIndex.begin(), fixIndex.end(), make_pair(j,i)) == fixIndex.end()) // (i,j) not fixed
					I.push_back(make_pair(j,i));
		}
	}
	cellErrorApprox(gfs,dt,I,heat,alpha,A,mu,u,p,r); 
	for(int i=0; i < dt.size(); i++)
		e[i]+=r[i];
	cellError(gfs,dt,I,heat,alpha,A,mu,u,p,r);
	for(int i=0; i < dt.size(); i++)
		e[i]-=r[i];
	
	
	for(int i=0; i < dt.size(); i++)
		e[i]*=0.5;
}


#endif 
