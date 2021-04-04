#include "wavefunction.h"
#include "constants.h"

CWaveFunction_classical::CWaveFunction_classical(){
}

double CWaveFunction_classical::CalcPsiSquared(double q,double r,double ctheta,double m1,double m2,int q1q2){
	double mu,psisquared=0.0,eratio,root,sign;
	mu=m1*m2/(m1+m2);
	eratio=2.0*mu*q1q2*ALPHA/(HBARC*r*q*q);  // ratio of PE to E
	if(eratio<1.0){
		if(fabs(eratio)<1.0E-7){
			psisquared=1.0;
		}
		else{
			root=1.0-2.0*eratio/(1.0+ctheta);
			if(root>0.0){
				root=sqrt(root);
				sign=1.0;
				psisquared =fabs(1.0+sign*eratio*eratio/((1.0+sign*root)*(1.0+sign*root)*root*(1.0+ctheta)*(1.0+ctheta)));
				sign=-1.0;
				psisquared+=fabs(1.0+sign*eratio*eratio/((1.0+sign*root)*(1.0+sign*root)*root*(1.0+ctheta)*(1.0+ctheta)));
			}
		}
	}
	return psisquared;
}

double CWaveFunction_classical::CalcPsiSquared(int iq,double r,double ctheta,double m1,double m2,int q1q2){
	return CalcPsiSquared((iq+0.5)*delq,r,ctheta,m1,m2,q1q2);
}
