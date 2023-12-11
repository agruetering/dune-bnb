#ifndef Dmax_HH
#define Dmax_HH

#include <dune/OCP/switchconstr.hh> 

/* switching interface for an upper bound on the number of switchings for a single switch*/
class Dmax : public SwitchPoly
{
protected: 
	const int sigma; // upper bound on the number of switchings 
public: 
	Dmax(const int ub) : sigma(ub) {}

	int separate(int,const double*, double*, double&) const; 
	/* compute a violated constraint, if any 
	 * parameters: 
	 * - current switching pattern
	 * - coefficients of violated constraint, if any
	 * - right-hand side of violated constraint, if any
	 * return value: 
	 * - 0: no violated constraint found 
	 * - 1: violated inequality found
	 */

	int optimize(int,const double*, double*) const; 
	/* linear optimization over switching polytope
	 * parameters: 
	 * - objective coefficients 
	 * - computed solution
	 * return value:
 	 * - status of the optimization
  	 *   -1: error
  	 *    0: problem solved
	 */

	int getSigma() const {return sigma;}
}; 

#endif
