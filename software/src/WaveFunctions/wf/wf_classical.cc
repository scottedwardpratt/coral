#include "msu_coral/wavefunction.h"
#include "msu_commonutils/constants.h"
#include "msu_commonutils/log.h"
using namespace NMSUPratt;

CWaveFunction_classical::CWaveFunction_classical(){
}

double CWaveFunction_classical::CalcPsiSquared(double q,double r,double ctheta,double m1,double m2,int q1q2){
	double psisquared=1.0;
	double mu,eratio,root,sign,Jacobian1,Jacobian2,q0ratio;
	mu=m1*m2/(m1+m2);
	eratio=2.0*mu*q1q2*ALPHA*HBARC/(r*q*q);  // ratio of PE to E
	if(q1q2==0 || fabs(eratio)<1.0E-8){
		psisquared=1.0;
	}
	else if(eratio!=eratio){
		psisquared=0.0;
	}
	else{
		root=1.0-2.0*eratio/(1.0+ctheta);
		if(root<=0.0 || root!=root){
			psisquared=0.0;
		}
		else{
			root=sqrt(root);
			q0ratio=sqrt(1.0-eratio);
			
			sign=1.0;
			
			Jacobian1=1.0 +sign*pow(eratio/((1.0+sign*root)*(1.0+ctheta)),2)/root;
			if(Jacobian1!=Jacobian1){
				snprintf(message,CLog::CHARLENGTH,"J1=%g, q0ratio=%g, eratio=%g\n",Jacobian1,q0ratio,eratio);
				CLog::Info(message);
				snprintf(message,CLog::CHARLENGTH,"ctheta=%g, root=%g\n",ctheta,root);
				CLog::Info(message);
				snprintf(message,CLog::CHARLENGTH,"q=%g, r=%g, q1q2=%d, ctheta=%g\n",q,r,q1q2,ctheta);
				CLog::Fatal(message);
			}
			//Avoid dividing 0/0, instead use limit for ctheta0=0
			sign=-1.0;
			Jacobian2=1.0 +sign*pow(eratio/((1.0+sign*root)*(1.0+ctheta)),2)/root;
			if(Jacobian2!=Jacobian2){
				snprintf(message,CLog::CHARLENGTH,"J2=%g, q0ratio=%g, eratio=%g\n",Jacobian2,q0ratio,eratio);
				CLog::Info(message);
				snprintf(message,CLog::CHARLENGTH,"ctheta=%g, root=%g\n",ctheta,root);
				CLog::Info(message);
				snprintf(message,CLog::CHARLENGTH,"q=%g, r=%g, q1q2=%d, ctheta=%g\n",q,r,q1q2,ctheta);
				CLog::Fatal(message);
			}
			
			psisquared=fabs(Jacobian1)+fabs(Jacobian2);
		}
		// there is an integrable singularity as root->0. To reduce noise, psisquared is cut off at 1000.0
		psisquared=psisquared/sqrt(1.0+0.000001*pow(psisquared-1.0,2));
	}
	// there is an integrable singularity as root->0. To reduce noise, psisquared is cut off at 100.0
	//psisquared=psisquared/sqrt(1.0+0.0001*pow(psisquared-1.0,2));
	return psisquared;
}

double CWaveFunction_classical::CalcPsiSquared(int iq,double r,double ctheta,double m1,double m2,int q1q2){
	return CalcPsiSquared((iq+0.5)*delq,r,ctheta,m1,m2,q1q2);
}
