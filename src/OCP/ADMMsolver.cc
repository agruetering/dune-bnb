#include <dune/OCP/ADMMsolver.hh> 


/* ADMM algorithm in scaled version 
 *
 * augmented Lagrangian: L(u,v,mu,w,nu)=f(u)+I_{(-inf,b]}(v)+I_{[0,1]}(w)
 * 					+ beta/2||Au-v+mu|| + rho/2*||u-w+nu||
 * with 
 * - mu scaled dual variable for cutting plane constraints (unscaled dual variable = beta*mu)
 * - nu scaled dual variable for bx constraints (unscaled dual variable = beta*nu)
*/

/* for Branch and bound without switching constraints */
template<class MType>
bool ADMMsolver<MType>::apply(ControlVector& u,Vector& w, Vector& nu,double& opt,const OptSense& optSense,const bool& primalBoundFound, const double& primalBound,const double& optEps, int &iterations) 
{
	bool fathom; 
	cout << "=== ADMM Method" << endl;

	Vector q(I.size());
	for(int i=0; i < I.size(); i++)
		q[i]=u[I[i].second][I[i].first];

	bool stop=false;
	int count=0;
	double C=0.0;
	for(int i=0; i < I.size(); i++) 
		C+=dt[I[i].second];

	while(!stop){
		cout << "\nADMM Iteration " << count << "\n" << endl; 
		
		/* calculate q(k+1)=(G*S*SG+(alpha+beta)I)^{-1}
		 *	\Chi_I[beta*(w(k)-nu(k))+G*S*(yd-Sf-SGu_fix)+alpha ud]-rho*A^TAu_fix
		 */
		Vector r(I.size()); 
		for(int i=0; i < I.size(); i++){
			r[i]=-Gstar_eta[I[i].second][I[i].first]+alpha*ud[I[i].second][I[i].first]-Gstar_q[I[i].second][I[i].first]; 
			r[i]+=beta*(w[i]-nu[i]); 
			r[i]*=dt[I[i].second];// r=\Chi_I beta*(w(k)-nu(k))+G*S*(yd-Sf-SGu_fix)+alpha ud  
		}
		
		// linear operator
		typedef LinOperator<Master> LinearOperator;
		LinearOperator op(gfs,dt,alpha,beta,I);
	
		ScalarProduct<Vector> sp;

		// CG solver
		Richardson<Vector,Vector> prec(1.0); // no preconditioning
		CGSolver<Vector> solver(op,sp,prec,reduction,iter_lin, 2);
	
		// storage for results
 		Dune::InverseOperatorResult result;
		// call CG solver
		solver.apply(q,r,result);
		
#ifdef DEBUG
		for(int i=0; i < dt.size(); i++){ 
			for(int j=0; j < n; j++){
				auto it=find(I.begin(), I.end(), make_pair(j,i)); 
				if(it==I.end())
					cout << u[i][j] << endl;
				else{
					int pos=it-I.begin();
					cout << q[pos] << endl;
				}
			}
		}
#endif
		

		/* w(k+1)=max( min(q(k+1)+nu(k),ub), ua) */
		Vector w_old=w; // save previous w(k) for stopping criteria
		for(int i=0; i < I.size(); i++) 
			w[i]=max(min((double) (q[i]+nu[i]),ub[I[i].second][I[i].first]), ua[I[i].second][I[i].first]); 

		/* nu(k+1)=nu(k)+q(k+1)-w(k+1) */
		nu+=q; 
		nu-=w;

		/* primal error estimation */	
		/* contribution --beta*nu(k+1)^T[q(k+1)-w(k+1)] */
		double sum=0;
		for(int i=0; i < I.size(); i++)
			sum-=beta*nu[i]*(q[i]-w[i])*dt[I[i].second];
		/* dual residual: calculate norm[beta(w(k)-w(k+1))]*/
		Vector h=w_old; 
		h-=w; 
		h*=beta; // h=beta*(w(k)-w(k+1))
		double d_residual=0.0;
    		for (int i=0; i < I.size(); i++){
       			d_residual+=pow(h[i],2)*dt[I[i].second];
		}
		/* contribution dual residual */
		sum+=sqrt(C*d_residual);

		ControlVector z=u;
		for(int i=0; i < I.size(); i++) 
			z[I[i].second][I[i].first]=q[i]; 
		primalObj(z,opt); 


		/* primal residual*/
		double p_residual=0.0; 
		for(int i=0; i < I.size(); i++) 	
			p_residual+=pow(q[i]-w[i],2)*dt[I[i].second];
		p_residual=sqrt(p_residual); 
		
		/* norm beta*nu(k+1) */
		double dnorm=0.0; 
		for(int i=0; i < I.size(); i++){
			dnorm+=pow(beta*nu[i],2)*dt[I[i].second]; 
		}
		dnorm=sqrt(dnorm);
		/* norm q(k+1) */
		double pnorm=0.0; 
		for(int i=0; i < I.size(); i++)
			pnorm+=pow(q[i],2)*dt[I[i].second]; 
		pnorm=sqrt(pnorm);
		/* norm w(k+1) */
		double norm=0.0; 
		for(int i=0; i < I.size(); i++)
			norm+=pow(w[i],2)*dt[I[i].second]; 
		norm=sqrt(norm);
			
		opt-=abs(sum); // dual bound
		if(p_residual<=eabs+erel*max(pnorm,norm) && d_residual<=eabs+erel*dnorm && sum<=tol)
			stop=true;
		
		// check if dual bound exceeds primal bound
		if(optSense == Minimize && primalBoundFound){
      			if((primalBound-opt)/primalBound <= optEps){
	 			stop=true;
				fathom=true;
			}
		}
   		if(optSense == Maximize && primalBoundFound){
      			if((opt-primalBound)/primalBound <= optEps){
	 			stop=true;
				fathom=true;
			}
		}

		count ++;
		
	}
	cout << "\n===" << "End ADMM Method \n" << endl;
	
	iterations=count;
	for(int i=0; i < I.size(); i++) 
		u[I[i].second][I[i].first]=q[i];

	return fathom;
}

/* for root node relaxation without switching constraints */
template<class MType>
void ADMMsolver<MType>::apply(ControlVector& u,Vector& w, Vector& nu,double& opt, int &iterations) 
{ 
	cout << "=== ADMM Method" << endl;

	Vector q(I.size());
	for(int i=0; i < I.size(); i++)
		q[i]=u[I[i].second][I[i].first];

	bool stop=false;
	int count=0;
	double C=0.0;
	for(int i=0; i < I.size(); i++) 
		C+=dt[I[i].second];

	while(!stop){
		cout << "\nADMM Iteration " << count << "\n" << endl; 
		
		/* calculate q(k+1)=(G*S*SG+(alpha+beta)I)^{-1}
		 *	\Chi_I[beta*(w(k)-nu(k))+G*S*(yd-Sf-SGu_fix)+alpha ud]-rho*A^TAu_fix
		 */
		Vector r(I.size()); 
		for(int i=0; i < I.size(); i++){
			r[i]=-Gstar_eta[I[i].second][I[i].first]+alpha*ud[I[i].second][I[i].first]-Gstar_q[I[i].second][I[i].first]; 
			r[i]+=beta*(w[i]-nu[i]); 
			r[i]*=dt[I[i].second];// r=\Chi_I beta*(w(k)-nu(k))+G*S*(yd-Sf-SGu_fix)+alpha ud  
		}
		
		// linear operator
		typedef LinOperator<Master> LinearOperator;
		LinearOperator op(gfs,dt,alpha,beta,I);
	
		ScalarProduct<Vector> sp;

		// CG solver
		Richardson<Vector,Vector> prec(1.0); // no preconditioning
		CGSolver<Vector> solver(op,sp,prec,reduction,iter_lin, 2);
	
		// storage for results
 		Dune::InverseOperatorResult result;
		// call CG solver
		solver.apply(q,r,result);
		
#ifdef DEBUG
		for(int i=0; i < dt.size(); i++){ 
			for(int j=0; j < n; j++){
				auto it=find(I.begin(), I.end(), make_pair(j,i)); 
				if(it==I.end())
					cout << u[i][j] << endl;
				else{
					int pos=it-I.begin();
					cout << q[pos] << endl;
				}
			}
		}
#endif
		

		/* w(k+1)=max( min(q(k+1)+nu(k),ub), ua) */
		Vector w_old=w; // save previous w(k) for stopping criteria
		for(int i=0; i < I.size(); i++) 
			w[i]=max(min((double) (q[i]+nu[i]),ub[I[i].second][I[i].first]), ua[I[i].second][I[i].first]); 

		/* nu(k+1)=nu(k)+q(k+1)-w(k+1) */
		nu+=q; 
		nu-=w;

		/* primal error estimation */	
		/* contribution --beta*nu(k+1)^T[q(k+1)-w(k+1)] */
		double sum=0;
		for(int i=0; i < I.size(); i++)
			sum-=beta*nu[i]*(q[i]-w[i])*dt[I[i].second];
		/* dual residual: calculate norm[beta(w(k)-w(k+1))]*/
		Vector h=w_old; 
		h-=w; 
		h*=beta; // h=beta*(w(k)-w(k+1))
		double d_residual=0.0;
    		for (int i=0; i < I.size(); i++){
       			d_residual+=pow(h[i],2)*dt[I[i].second];
		}
		/* contribution dual residual */
		sum+=sqrt(C*d_residual);

		ControlVector z=u;
		for(int i=0; i < I.size(); i++) 
			z[I[i].second][I[i].first]=q[i]; 
		primalObj(z,opt); 


		/* primal residual*/
		double p_residual=0.0; 
		for(int i=0; i < I.size(); i++) 	
			p_residual+=pow(q[i]-w[i],2)*dt[I[i].second];
		p_residual=sqrt(p_residual); 
		
		/* norm beta*nu(k+1) */
		double dnorm=0.0; 
		for(int i=0; i < I.size(); i++){
			dnorm+=pow(beta*nu[i],2)*dt[I[i].second]; 
		}
		dnorm=sqrt(dnorm);
		/* norm q(k+1) */
		double pnorm=0.0; 
		for(int i=0; i < I.size(); i++)
			pnorm+=pow(q[i],2)*dt[I[i].second]; 
		pnorm=sqrt(pnorm);
		/* norm w(k+1) */
		double norm=0.0; 
		for(int i=0; i < I.size(); i++)
			norm+=pow(w[i],2)*dt[I[i].second]; 
		norm=sqrt(norm);
			
		opt-=abs(sum); // dual bound
		if(p_residual<=eabs+erel*max(pnorm,norm) && d_residual<=eabs+erel*dnorm && sum<=tol)
			stop=true;

		count ++;
		
	}
	cout << "\n===" << "End ADMM Method \n" << endl;
	
	iterations=count;
	for(int i=0; i < I.size(); i++) 
		u[I[i].second][I[i].first]=q[i];
}

/* for Branch and bound with switching constraints */
template<class MType>
bool ADMMsolver<MType>::apply(ControlVector& u, Vector& v, Vector& mu,Vector& w, Vector& nu, const Matrix &A, const Vector& b,double& opt,const OptSense& optSense,const bool& primalBoundFound, const double& primalBound,const double& optEps, int &iterations) 
{
	bool fathom; 
	cout << "=== ADMM Method" << endl;
	// Preconditioner
	Dune::Matrix<double> ATA(I.size(),I.size()); 
	for(int i=0; i < I.size(); i++){
		for(int j=0; j < I.size(); j++){
			ATA[i][j]=0; 
			int k=I[i].second*n+I[i].first;
			int h=I[j].second*n+I[j].first;
			for(int l=0; l < A.N(); l++){
				if(A.exists(l,k) && A.exists(l,h))
					ATA[i][j]+=A[l][k]*A[l][h];
			}
		}
	}
	Matrix P(I.size(),I.size(),Matrix::row_wise);
	for(auto row=P.createbegin(); row!=P.createend(); ++row){
		row.insert(row.index());
		for(int j=0; j < I.size(); j++){
			if(ATA[row.index()][j]!=0)
				row.insert(j);
		}
	}
	
	for(int i=0; i < I.size(); i++){
		for(int j=0; j < I.size(); j++){
			if(ATA[i][j]!=0)
				P[i][j]=rho*ATA[i][j]; 
			if(j==i)
				P[i][j]+=(alpha+beta)*dt[I[i].second]; 
		}
	}
	InverseP<Matrix,Vector,Vector> invOperator(P);
	Precond prec(invOperator);

	Vector q(I.size());
	for(int i=0; i < I.size(); i++)
		q[i]=u[I[i].second][I[i].first];

	bool stop=false;
	int count=0;

	Vector d(A.N());
	if(I.size()!=dt.size()*n){ // fixed control variables existent
		Vector h(dt.size()*n); 	
		for(int i=0; i < dt.size(); i++){
			for(int j=0; j < n; j++){
				if (find(I.begin(), I.end(), make_pair(j,i)) == I.end()) // (j,i) fixed
					h[i*n+j]=u[i][j];
			}
		}
		A.mv(h,d);		
	}

	double C=0.0;
	for(int i=0; i < I.size(); i++) 
		C+=dt[I[i].second];

	while(!stop){
		cout << "\nADMM Iteration " << count << "\n" << endl; 
		
		/* calculate q(k+1)=(G*S*SG+(alpha+beta)I+rho*A^TA)^{-1}
		 *	\Chi_I[rho*A^T(v(k)-mu(k))+beta*(w(k)-nu(k))+G*S*(yd-Sf-SGu_fix)+alpha ud]
		 */
		Vector h1(A.M()); // auxiliary variable to save A^T(-mu(k)+v(k))
		Vector h2(A.N()); 

		h2-=mu;
		h2+=v; // h2=-mu(k)+v(k)
		A.usmtv(rho,h2,h1); // h1=rho*A^T(-mu(k)+v(k))

		Vector r(I.size()); 
		for(int i=0; i < I.size(); i++){
			r[i]=-Gstar_eta[I[i].second][I[i].first]+alpha*ud[I[i].second][I[i].first]-Gstar_q[I[i].second][I[i].first]; 
			r[i]+=beta*(w[i]-nu[i]); 
			r[i]*=dt[I[i].second];// r=\Chi_I beta*(w(k)-nu(k))+G*S*(yd-Sf-SGu_fix)+alpha ud  
			r[i]+=h1[I[i].second*n+I[i].first]; // r+=\Chi_I rho*A^T(-mu(k)+v(k))
		}
		
		// linear operator
		typedef LinCutOperator<Master> LinearOperator;
		LinearOperator op(gfs,A,dt,alpha,rho,beta,I);
	
		ScalarProduct<Vector> sp;

		// CG solver
		CGSolver<Vector> solver(op,sp,prec,reduction,iter_lin, 2);
	
		// storage for results
 		Dune::InverseOperatorResult result;
		// call CG solver
		solver.apply(q,r,result);
		
#ifdef DEBUG
		for(int i=0; i < dt.size(); i++){ 
			for(int j=0; j < n; j++){
				auto it=find(I.begin(), I.end(), make_pair(j,i)); 
				if(it==I.end())
					cout << u[i][j] << endl;
				else{
					int pos=it-I.begin();
					cout << q[pos] << endl;
				}
			}
		}
#endif
		
		/* v(k+1)=min(Aq(k+1)+p(k),b-Au_fix) */
		Vector v_old=v; // save previous v(k) for stopping criteria
		
		h1=0; 
		for(int i=0; i < I.size(); i++) 
			h1[I[i].second*n+I[i].first]=q[i];
			
		A.mv(h1,h2); // h2=Aq(k+1)
		for(int i=0; i < A.N(); i++)
			v[i]=min(h2[i]+mu[i],b[i]-d[i]);
		
		/* mu(k+1)=mu(k)+Aq(k+1)-v */
		mu+=h2; 
		mu-=v;	

		/* w(k+1)=max( min(q(k+1)+nu(k),ub), ua) */
		Vector w_old=w; // save previous w(k) for stopping criteria
		for(int i=0; i < I.size(); i++) 
			w[i]=max(min((double) (q[i]+nu[i]),ub[I[i].second][I[i].first]), ua[I[i].second][I[i].first]); 

		/* nu(k+1)=nu(k)+q(k+1)-w(k+1) */
		nu+=q; 
		nu-=w;

		/* primal error estimation */	
		/* contribution -rho*mu(k+1)^T[Aq(k+1)-v(k+1)]-beta*nu(k+1)^T[q(k+1)-w(k+1)] */
		double sum=0;
		sum-=mu.dot(h2);
		sum+=mu.dot(v);
		sum*=rho;
		for(int i=0; i < I.size(); i++)
			sum-=beta*nu[i]*(q[i]-w[i])*dt[I[i].second];
		/* dual residual: calculate norm \Chi_I rho*A^T(v(k)-v(k+1)) + beta(w(k)-w(k+1))*/
		Vector h=v_old; 
		h-=v; 
		h1=0;
		A.usmtv(rho,h,h1); // h1=rho*A^T(v(k)-v(k+1))
		h.resize(I.size());
		h=w_old; 
		h-=w; 
		h*=beta; // h=beta*(w(k)-w(k+1))
		double d_residual=0.0;
    		for (int i=0; i < I.size(); i++){
			d_residual+=pow(h1[I[i].second*n+I[i].first],2)/dt[I[i].second];
       			d_residual+=pow(h[i],2)*dt[I[i].second];
			d_residual+=2*h1[I[i].second*n+I[i].first]*h[i];
		}
		/* contribution dual residual */
		sum+=sqrt(C*d_residual);


		ControlVector z=u;
		for(int i=0; i < I.size(); i++) 
			z[I[i].second][I[i].first]=q[i]; 
		primalObj(z,opt); 
		
		/* primal residual*/
		double p_residual=0.0; 
		h.resize(A.N()); 
		h=h2; 
		h-=v; 
		p_residual+=h.two_norm(); 
		double val=0.0;
		for(int i=0; i < I.size(); i++) 	
			val+=pow(q[i]-w[i],2)*dt[I[i].second];
		p_residual+=sqrt(val); 
		
		/* norm rho*A^Tmu(k+1)+beta*nu(k+1) */
		double dnorm=0.0; 
		A.usmtv(rho,mu,h1);
		for(int i=0; i < I.size(); i++){
			dnorm+=pow(beta*nu[i],2)*dt[I[i].second]; 
			dnorm+=pow(h1[I[i].second*n+I[i].first],2)/dt[I[i].second];
		}
		dnorm+=sqrt(dnorm);
		/* norm [Aq(k+1),q(k+1)] */
		double pnorm=0.0; 
		pnorm+=h2.two_norm();
		val=0.0;
		for(int i=0; i < I.size(); i++)
			val+=pow(q[i],2)*dt[I[i].second]; 
		pnorm+=sqrt(val);
		/* norm [v(k+1),w(k+1)] */
		double norm=0.0; 
		norm+=v.two_norm(); 
		val=0.0;
		for(int i=0; i < I.size(); i++)
			val+=pow(w[i],2)*dt[I[i].second]; 
		norm+=sqrt(val);
			
		opt-=abs(sum); // dual bound
		if(p_residual<=(sqrt(A.N())+1)*eabs+erel*max(pnorm,norm) && d_residual<=eabs+erel*dnorm && sum<=tol)
			stop=true;
		
		// check if dual bound exceeds primal bound
		if(optSense == Minimize && primalBoundFound){
      			if((primalBound-opt)/primalBound <= optEps){
	 			stop=true;
				fathom=true;
			}
		}
   		if(optSense == Maximize && primalBoundFound){
      			if((opt-primalBound)/primalBound <= optEps){
	 			stop=true;
				fathom=true;
			}
		}

		count ++;
		
	}
	cout << "\n===" << "End ADMM Method \n" << endl;
	
	iterations=count;
	for(int i=0; i < I.size(); i++) 
		u[I[i].second][I[i].first]=q[i];

	return fathom;
}


/* for root node relaxation with switching constraints */
template<class MType>
bool ADMMsolver<MType>::apply(ControlVector& u, Vector& v, Vector& mu,Vector& w, Vector& nu, const Matrix &A, const Vector& b,double& opt, int &iterations, double timelimit) 
{
	bool fathom;
	cout << "=== ADMM Method" << endl;
	// Preconditioner
	Dune::Matrix<double> ATA(I.size(),I.size()); 
	for(int i=0; i < I.size(); i++){
		for(int j=0; j < I.size(); j++){
			ATA[i][j]=0; 
			int k=I[i].second*n+I[i].first;
			int h=I[j].second*n+I[j].first;
			for(int l=0; l < A.N(); l++){
				if(A.exists(l,k) && A.exists(l,h))
					ATA[i][j]+=A[l][k]*A[l][h];
			}
		}
	}
	Matrix P(I.size(),I.size(),Matrix::row_wise);
	for(auto row=P.createbegin(); row!=P.createend(); ++row){
		row.insert(row.index());
		for(int j=0; j < I.size(); j++){
			if(ATA[row.index()][j]!=0)
				row.insert(j);
		}
	}
	
	for(int i=0; i < I.size(); i++){
		for(int j=0; j < I.size(); j++){
			if(ATA[i][j]!=0)
				P[i][j]=rho*ATA[i][j]; 
			if(j==i)
				P[i][j]+=(alpha+beta)*dt[I[i].second]; 
		}
	}
	InverseP<Matrix,Vector,Vector> invOperator(P);
	Precond prec(invOperator);

	Vector q(I.size());
	for(int i=0; i < I.size(); i++)
		q[i]=u[I[i].second][I[i].first];

	bool stop=false;
	int count=0;
	double time=0.0;

	Vector d(A.N());
	if(I.size()!=dt.size()*n){ // fixed control variables existent
		Vector h(dt.size()*n); 	
		for(int i=0; i < dt.size(); i++){
			for(int j=0; j < n; j++){
				if (find(I.begin(), I.end(), make_pair(j,i)) == I.end()) // (j,i) fixed
					h[i*n+j]=u[i][j];
			}
		}
		A.mv(h,d);		
	}

	double C=0.0;
	for(int i=0; i < I.size(); i++) 
		C+=dt[I[i].second];

	while(!stop && time<timelimit){
		clock_t cstart=clock();
		cout << "\nADMM Iteration " << count << "\n" << endl; 
		
		/* calculate q(k+1)=(G*S*SG+(alpha+beta)I+rho*A^TA)^{-1}
		 *	\Chi_I[rho*A^T(v(k)-mu(k))+beta*(w(k)-nu(k))+G*S*(yd-Sf-SGu_fix)+alpha ud]
		 */
		Vector h1(A.M()); // auxiliary variable to save A^T(-mu(k)+v(k))
		Vector h2(A.N()); 

		h2-=mu;
		h2+=v; // h2=-mu(k)+v(k)
		A.usmtv(rho,h2,h1); // h1=rho*A^T(-mu(k)+v(k))

		Vector r(I.size()); 
		for(int i=0; i < I.size(); i++){
			r[i]=-Gstar_eta[I[i].second][I[i].first]+alpha*ud[I[i].second][I[i].first]-Gstar_q[I[i].second][I[i].first]; 
			r[i]+=beta*(w[i]-nu[i]); 
			r[i]*=dt[I[i].second];// r=\Chi_I beta*(w(k)-nu(k))+G*S*(yd-Sf-SGu_fix)+alpha ud  
			r[i]+=h1[I[i].second*n+I[i].first]; // r+=\Chi_I rho*A^T(-mu(k)+v(k))
		}
		
		// linear operator
		typedef LinCutOperator<Master> LinearOperator;
		LinearOperator op(gfs,A,dt,alpha,rho,beta,I);
	
		ScalarProduct<Vector> sp;

		// CG solver
		CGSolver<Vector> solver(op,sp,prec,reduction,iter_lin, 2);
	
		// storage for results
 		Dune::InverseOperatorResult result;
		// call CG solver
		solver.apply(q,r,result);
		
#ifdef DEBUG
		for(int i=0; i < dt.size(); i++){ 
			for(int j=0; j < n; j++){
				auto it=find(I.begin(), I.end(), make_pair(j,i)); 
				if(it==I.end())
					cout << u[i][j] << endl;
				else{
					int pos=it-I.begin();
					cout << q[pos] << endl;
				}
			}
		}
#endif
		
		/* v(k+1)=min(Aq(k+1)+p(k),b-Au_fix) */
		Vector v_old=v; // save previous v(k) for stopping criteria
		
		h1=0; 
		for(int i=0; i < I.size(); i++) 
			h1[I[i].second*n+I[i].first]=q[i];
			
		A.mv(h1,h2); // h2=Aq(k+1)
		for(int i=0; i < A.N(); i++)
			v[i]=min(h2[i]+mu[i],b[i]-d[i]);
		
		/* mu(k+1)=mu(k)+Aq(k+1)-v */
		mu+=h2; 
		mu-=v;	

		/* w(k+1)=max( min(q(k+1)+nu(k),ub), ua) */
		Vector w_old=w; // save previous w(k) for stopping criteria
		for(int i=0; i < I.size(); i++) 
			w[i]=max(min((double) (q[i]+nu[i]),ub[I[i].second][I[i].first]), ua[I[i].second][I[i].first]); 

		/* nu(k+1)=nu(k)+q(k+1)-w(k+1) */
		nu+=q; 
		nu-=w;

		/* primal error estimation */	
		/* contribution -rho*mu(k+1)^T[Aq(k+1)-v(k+1)]-beta*nu(k+1)^T[q(k+1)-w(k+1)] */
		double sum=0;
		sum-=mu.dot(h2);
		sum+=mu.dot(v);
		sum*=rho;
		for(int i=0; i < I.size(); i++)
			sum-=beta*nu[i]*(q[i]-w[i])*dt[I[i].second];
		/* dual residual: calculate norm \Chi_I rho*A^T(v(k)-v(k+1)) + beta(w(k)-w(k+1))*/
		Vector h=v_old; 
		h-=v; 
		h1=0;
		A.usmtv(rho,h,h1); // h1=rho*A^T(v(k)-v(k+1))
		h.resize(I.size());
		h=w_old; 
		h-=w; 
		h*=beta; // h=beta*(w(k)-w(k+1))
		double d_residual=0.0;
    		for (int i=0; i < I.size(); i++){
			d_residual+=pow(h1[I[i].second*n+I[i].first],2)/dt[I[i].second];
       			d_residual+=pow(h[i],2)*dt[I[i].second];
			d_residual+=2*h1[I[i].second*n+I[i].first]*h[i];
		}
		/* contribution dual residual */
		sum+=sqrt(C*d_residual);


		ControlVector z=u;
		for(int i=0; i < I.size(); i++) 
			z[I[i].second][I[i].first]=q[i]; 
		primalObj(z,opt); 
		
		/* primal residual*/
		double p_residual=0.0; 
		h.resize(A.N()); 
		h=h2; 
		h-=v; 
		p_residual+=h.two_norm(); 
		double val=0.0;
		for(int i=0; i < I.size(); i++) 	
			val+=pow(q[i]-w[i],2)*dt[I[i].second];
		p_residual+=sqrt(val); 
		
		/* norm rho*A^Tmu(k+1)+beta*nu(k+1) */
		double dnorm=0.0; 
		A.usmtv(rho,mu,h1);
		for(int i=0; i < I.size(); i++){
			dnorm+=pow(beta*nu[i],2)*dt[I[i].second]; 
			dnorm+=pow(h1[I[i].second*n+I[i].first],2)/dt[I[i].second];
		}
		dnorm+=sqrt(dnorm);
		/* norm [Aq(k+1),q(k+1)] */
		double pnorm=0.0; 
		pnorm+=h2.two_norm();
		val=0.0;
		for(int i=0; i < I.size(); i++)
			val+=pow(q[i],2)*dt[I[i].second]; 
		pnorm+=sqrt(val);
		/* norm [v(k+1),w(k+1)] */
		double norm=0.0; 
		norm+=v.two_norm(); 
		val=0.0;
		for(int i=0; i < I.size(); i++)
			val+=pow(w[i],2)*dt[I[i].second]; 
		norm+=sqrt(val);
			
		opt-=abs(sum); // dual bound
		if(p_residual<=(sqrt(A.N())+1)*eabs+erel*max(pnorm,norm) && d_residual<=eabs+erel*dnorm && sum<=tol)
			stop=true;

		count ++;
		time+=(clock()-cstart)/CLOCKS_PER_SEC;
		
	}
	cout << "\n===" << "End ADMM Method \n" << endl;
	
	iterations=count;
	for(int i=0; i < I.size(); i++) 
		u[I[i].second][I[i].first]=q[i];

	if(!stop && time > timelimit)
		return 1; 
	else 
		return 0;
}



// primal objective value (if switching constraints existent) 
template<class MType>
void ADMMsolver<MType>::primalObj(const ControlVector& u,double& opt)
{
	/* calculate objective value*/
	// y=S(Gu+f)
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
	opt=0;
	L2Deviation(gfs,dt_new,y_ex,YdProblem,opt); 
    	// distribution alpha/2 \int_[0,T] ||u-u_d||^2 dt 
    	for (int i=0; i < dt_new.size(); i++){
		double time=accumulate(dt_new.begin(),dt_new.begin()+i+1,0.0);
		auto diff=u_ex[i]-heat.ud(time); 
       		opt+=diff.two_norm2()*alpha*dt_new[i]*0.5; 
	}
	
	
}
	

 
