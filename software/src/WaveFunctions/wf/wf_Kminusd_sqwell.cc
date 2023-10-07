#include "msu_coral/wavefunction.h"
#include "msu_commonutils/constants.h"
#include "msu_commonutils/sf.h"
using namespace std;
using namespace NMSUPratt;

CWaveFunction_Kminusd_sqwell::CWaveFunction_Kminusd_sqwell(string parsfilename) : CWaveFunction(){
	CLog::Info("Beware! The K^-d wavefunction was tuned to match only the scattering length and ignores the inelastic process.\n The K-d system is inelastic at threshold due to the chance of producing a Lambda.\n The deuteron breakup channel also comes into play for q>41 MeV. So this treatment is very questionable!\n");
	// scattering length comes from M.Doring and U.-G. Meissner, PLB 704, p. 663-666 (2011).
  ParsInit(parsfilename);

  m1=PionMass; 
  m2=ProtonMass+NeutronMass-2.224;
  IDENTICAL=0;

  q1q2=-1;
  nchannels=1;

  ellmax=0;
  InitArrays();
  CLog::Info("Arrays Initialized\n");

  ell[0]=0;
  InitWaves();

  nwells=new int[nchannels];
  nwells[0]=1;

  SquareWell_MakeArrays();
	
  a[0][0]=1.0; 

  V0[0][0]=-78.882;

  SquareWell_Init();
}

CWaveFunction_Kminusd_sqwell::~CWaveFunction_Kminusd_sqwell(){
  SquareWell_DeleteArrays();
}


double CWaveFunction_Kminusd_sqwell::CalcPsiSquared(int iq,double r,double ctheta){
  double psisquared,theta=acos(ctheta),x;
  double q=GetQ(iq);
  complex<double> psi,psia,Xlm00;
  psia=planewave[iq]->planewave(r,ctheta);
  complex<double> DelPhi[2];

  SquareWell_GetDelPhi(iq,r,DelPhi);
	x=q*r/HBARC;
  Xlm00=sqrt(4.0*PI)*SpherHarmonics::Ylm(0,0,theta,0.0)/x;

  psi=psia+Xlm00*DelPhi[0];
  psisquared=real(psi*conj(psi));

	//psisquared*=RelativisticCorrection(r,iq);
  return psisquared;

}
