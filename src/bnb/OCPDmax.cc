#include <dune/bnb/OCPDmax.hh>

/* separation of current control for subproblem */
template<class GridView, class Constraints, class Problem, int k>
int OCPDmax<GridView,Constraints,Problem,k>::separate(Subproblem* sub,const vector<double>& dt,const ControlVector& u, double* cut, double& rhs) 
{
	int status;
	// store current control 
	double* v=new double[dt.size()];
	for(int i=0; i<dt.size(); i++){
		v[i]=u[i];
	}
	if(sub->getNumFixPoints()==0)
		status = this->D.separate(dt.size(),v,cut,rhs); 
	else{
		// extend control vector by fixations
		int nFix=sub->getNumFixPoints();
		double* Ev=new double[dt.size()+nFix]; 
		vector<pair<int,double>> tau=sub->getTau(); // fixation time points
		double endTime=accumulate(dt.begin(),dt.end(),0.0);
		tau.push_back(pair(0,endTime+1));

		vector<double> c=sub->getFix(); // fixation values

		int j=0; // index of fixations
		int i=0; // index of v
		double time=0;
		for(int l=0; l < nFix+dt.size(); l++){
			if(i<dt.size() && time+dt[i]<=tau[j].second+1e-10){
				Ev[l]=v[i];
				time+=dt[i];
				i++; 
			}
			else{
				Ev[l]=c[j];
				j++;
			}
		}

		// separate extended vector 
		double* Ecut=new double[dt.size()+nFix];
		status = this->D.separate(dt.size()+nFix, Ev, Ecut, rhs); 
		if(status!=0){
			// reconstruct cutting plane for v
			j=0; 
			i=0;
			time=0;
			for(int l=0; l < nFix+dt.size(); l++){
				if(i<dt.size() && time+dt[i]<=tau[j].second+1e-10){
					cut[i]=Ecut[l];
					time+=dt[i];
					i++; 
				}
				else{
					rhs-=Ecut[l]*c[j];
					j++;
				}
			}
		}

		delete[] Ev; 
		delete[] Ecut;

		
	}
	delete[] v;

	return status;
}


/* heuristic solution*/
template<class GridView, class Constraints, class Problem, int k>
void OCPDmax<GridView,Constraints,Problem,k>::heuristic(Subproblem* sub,const GFS& gfs,const vector<double>& dt,const ControlVector& u,ControlVector& v)
{
	// fixations
	int nFix=sub->getNumFixPoints(); 
	vector<pair<int,double>> tau=sub->getTau();
	vector<double> d=sub->getFix();
	double endTime=accumulate(dt.begin(),dt.end(),0.0);
	tau.push_back(pair(0,endTime+1));
	
	int Nt=dt.size();

	// set up objective function for extended vector
	double* c=new double[Nt+nFix];
	int i=0; // time index
	int j=0; // index of fixations 
	double time=0; // time 
	for(int l=0; l < Nt+nFix; l++){
		if(i<Nt && time+dt[i]<=tau[j].second+1e-10){
			c[l]=0.5-u[0][i];
			c[l]*=dt[i];
			time+=dt[i];
			i++; 
		}
		else{
			// force fixations in heuristic solution
			if(d[j]==0)
				c[l]=(nFix+1)*endTime; 
			else 
				c[l]=-1*endTime; 
			j++;
		}
	}
	// calculate extended heuristic solution 
	double* w=new double[Nt+nFix];
	this->D.optimize(Nt+nFix,c,w);
	// reconstruct heuristic solution 
	j=0; 
	i=0;
	time=0;
	for(int l=0; l < nFix+Nt; l++){
		if(i<Nt && time+dt[i]<=tau[j].second+1e-10){
			v[i]=w[l];
			time+=dt[i];
			i++; 
		}
		else
			j++;
	}

	delete[] c; 
	delete[] w;
}

/* update fixations of control variables */
template<class GridView, class Constraints, class Problem, int k>
int OCPDmax<GridView,Constraints,Problem,k>::updateFixControls(Subproblem* sub,const vector<double>& dt,vector<pair<int,int>>& fixIndices, ControlVector& u)
{
	
	vector<double> c=sub->getFix();
	int numPoints=sub->getNumFixPoints();
	if(numPoints==0)
		return 0;

	 // store total variation of fixations 
	double sum=c[0];
	for(int i=1; i < numPoints; i++) 
		sum+=abs(c[i]-c[i-1]);

	// update fixations
	vector<pair<int,int>> fixations;
	// no fixations possible		
	if(sum < sigma-1){ 
		if(fixations!=fixIndices){ // fixed indices compared to father node changed
			fixIndices=fixations;
			return 1;
		}
		else 
			return 0;
	}

	vector<pair<int,double>> tau=sub-> getTau();
	// store time grid points
	vector<double> tpoints(dt.size());
	tpoints[0]=dt[0];
	for(int i=1; i < dt.size(); i++)
		tpoints[i]=tpoints[i-1]+dt[i]; 	
		
	for(int i=0; i < numPoints-1; i++){
		if(c[i+1]==c[i]){ // no jump in (tau[i],tau[i+1]] possible
		
			// start element that is not less than tau[i]
			auto start=lower_bound(tpoints.begin(), tpoints.end(), tau[i].second,[](const double& a, double value){return a < value+1e-10;}); 
			// corresponding start index
			int index_start=distance(tpoints.begin(), start);
			// last element that is not greater than tau[i+1]
			auto end=lower_bound(tpoints.begin(), tpoints.end(), tau[i+1].second,[](const double& a, double value){return a < value+1e-10;}); 
			// corresponding end index
			int index_end=distance(tpoints.begin(), end);

			// update fixed control indices and values
			for(int j=index_start; j< index_end; j++){
				fixations.push_back(pair(0,j));
				u[j]=c[i];
			}
		}
	}
	// no jump in (0,tau[i]] if c[0]=0
	if(c[0]==0){
		auto it=lower_bound(tpoints.begin(), tpoints.end(), tau[0].second,[](const double& a, double value){return a < value+1e-10;}); 
		int index = distance(tpoints.begin(), it);
		for(int j=index-1; j>=0; j--){ 
			fixations.insert(fixations.begin(),pair(0,j));
			u[j]=c[0];
		}
	}
	if(sum>=sigma){ // no jump in (tau[end],T] possible
		// no jump in (tau[end],T]
		auto it=lower_bound(tpoints.begin(), tpoints.end(), tau[numPoints-1].second,[](const double& a, double value){return a < value+1e-10;}); 
		int index = distance(tpoints.begin(), it);
		for(int j=index; j < dt.size(); j++){ 
			fixations.push_back(pair(0,j));
			u[j]=c[numPoints-1];
		}
	}

	if(fixations!=fixIndices){ // fixed indices compared to father node changed
		fixIndices=fixations;
		return 1; 
	}
	else 
		return 0;
		
}


	
