#include "msu_coral/wavefunction.h"
#include "msu_commonutils/constants.h"
#include "msu_commonutils/sf.h"
using namespace std;

CWaveFunction_pd_sqwell::CWaveFunction_pd_sqwell(string parsfilename) : CWaveFunction(){
  ParsInit(parsfilename);

  m1=MPROTON; 
  m2=2*MPROTON-2.224;
  IDENTICAL=0;

  q1q2=1;
  nchannels=2;

  ellmax=0;
  InitArrays();
  CLog::Info("Arrays Initialized\n");

  ell[0]=0;
  ell[1]=0;
  InitWaves();

  nwells=new int[nchannels];
  nwells[0]=3;
  nwells[1]=1;

  SquareWell_MakeArrays();
	
  a[0][0]=3.23445; a[0][1]=5.57680; a[0][2]=7.71351; // S11, I=1, J=1/2
  a[1][0]=4.99965; // S13, I=1, J=3/2

  V0[0][0]=-28.005; V0[0][1]=37.283; V0[0][2]=-5.829;
  V0[1][0]=108.203;
	//V0[0][0]=V0[0][1]=V0[0][2]=V0[1][0]=0.0;

  SquareWell_Init();
}

CWaveFunction_pd_sqwell::~CWaveFunction_pd_sqwell(){
  SquareWell_DeleteArrays();
}


double CWaveFunction_pd_sqwell::CalcPsiSquared(int iq,double r,double ctheta){
  double psisquared,theta=acos(ctheta),x;
  double q=GetQ(iq);
  complex<double> psi,psia,Xlm00;
  psia=planewave[iq]->planewave(r,ctheta);
  complex<double> DelPhi[2];

  SquareWell_GetDelPhi(iq,r,DelPhi);
	x=q*r/HBARC;
  Xlm00=sqrt(4.0*PI)*SpherHarmonics::Ylm(0,0,theta,0.0)/x;
	// for S=1/2
  psi=psia+Xlm00*DelPhi[0];
  psisquared=real(psi*conj(psi))/3.0;
	// for S=3/2
  psi=psia+Xlm00*DelPhi[1];
  psisquared+=2.0*real(psi*conj(psi))/3.0;

	//psisquared*=RelativisticCorrection(r,iq);
  return psisquared;

}
