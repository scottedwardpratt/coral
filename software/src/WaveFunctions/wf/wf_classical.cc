#include "wavefunction.h"
#include "constants.h"

CWaveFunction_classical::CWaveFunction_classical(){
}

double CWaveFunction_classical::CalcPsiSquared(double q,double r,double ctheta,double m1,double m2,int q1q2){
	double psisquared=1.0;
	if(r>1000.0){
		psisquared=1.0;
	}
	else{
		double mu,eratio,root,sign,stheta,stheta0,ctheta0,Jacobian;
		mu=m1*m2/(m1+m2);
		eratio=2.0*mu*q1q2*ALPHA*HBARC/(r*q*q);  // ratio of PE to E
		//printf("mu=%g, eratio=%g, q=%g, r=%g, q1q2=%d, ALPHA=%g, HBARC=%g\n",mu,eratio,q,r,q1q2,ALPHA,HBARC);
		root=1.0-2.0*eratio/(1.0+ctheta);
		if(root>0.0){
			root=sqrt(root);
			stheta=sqrt(1.0-ctheta*ctheta);
			sign=1.0;
			stheta0=0.5*stheta*(1.0+sign*root);
			ctheta0=sqrt(1.0-stheta0*stheta0);
			Jacobian=(0.5/ctheta0)*(ctheta+sign*(ctheta-eratio)/root);
			Jacobian*=0.5*(1.0+sign*root);	
			psisquared=fabs(Jacobian);
			sign=-1.0;
			stheta0=0.5*stheta*(1.0+sign*root);
			ctheta0=sqrt(1.0-stheta0*stheta0);
			Jacobian=(0.5/ctheta0)*(ctheta+sign*(ctheta-eratio)/root);
			Jacobian*=0.5*(1.0+sign*root);
			psisquared+=fabs(Jacobian);
			psisquared*=sqrt(1.0-eratio);
			//printf("psisquared=%g\n",psisquared);
		}
		else
			psisquared=0.0;
	}
	return psisquared;
}

double CWaveFunction_classical::CalcPsiSquared(int iq,double r,double ctheta,double m1,double m2,int q1q2){
	return CalcPsiSquared((iq+0.5)*delq,r,ctheta,m1,m2,q1q2);
}
