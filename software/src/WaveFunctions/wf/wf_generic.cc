#include "msu_coral/wavefunction.h"
using namespace NMSUPratt;

CWaveFunction_generic::CWaveFunction_generic(string  parsfilename,int q1q2set,double m1set,
double m2set,double symmweightset): CWaveFunction(){
	generic=1;
	ParsInit(parsfilename);
	m1=m1set;
	m2=m2set;
  muscale=m1*m2/(m1+m2);
  mu=muscale;
  symmweight=symmweightset;
  q1q2scale=q1q2set;
  q1q2=q1q2scale;
	if(q1q2!=0)
		COULOMB=true;
	else
		COULOMB=false;
	STRONG=false;
	KILLSYM=false;
  nchannels=0;
  ellmax=0;
  InitArrays();
  InitWaves();
}

void CWaveFunction_generic::reset(int q1q2set,double m1set,double m2set,
				  double symmweightset){
  m1=m1set;
  m2=m2set;
  mu=m1*m2/(m1+m2);
  if(q1q2*q1q2set<0){
	  snprintf(message,CLog::CHARLENGTH,"Illegal: Trying to reset q1q2 to opposite charge\n");
	  CLog::Fatal(message);
  }
  q1q2=q1q2set;
  if(q1q2!=0)
		COULOMB=true;
	else
		COULOMB=false;
  symmweight=symmweightset;
}

double CWaveFunction_generic::CalcPsiSquared(int iq,double r,double ctheta){
  double psisquared,asymmweight;
  complex<double> psi1,psi2,psisymm,psiasymm;
  const double ROOT2=sqrt(2.0);
	
  if(iq>=nqmax){
    snprintf(message,CLog::CHARLENGTH,"iq too large =%d, nqmax=%d\n",iq,nqmax);
	 CLog::Info(message);
    psisquared=1.0;
  }
  else{
    psi1=planewave[iq]->planewave(r,ctheta);
		if(KILLSYM)
			psisquared=real(psi1*conj(psi1));
		else{
			psi2=conj(psi1);
			asymmweight=1.0-symmweight;
			psisymm=(psi1+psi2)/ROOT2;
			psiasymm=(psi1-psi2)/ROOT2;
	    psisquared=symmweight*real(psisymm*conj(psisymm))
	      +asymmweight*real(psiasymm*conj(psiasymm));
		}
    KILLSYM=false;
  }
	if(q1q2!=0) psisquared*=RelativisticCorrection(r,iq);
  return psisquared;

}
