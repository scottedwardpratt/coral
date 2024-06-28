#include "msu_coral/coral.h"
#include <vector>
using namespace std;
using namespace NMSUPratt;

int main(){
  CWaveFunction_pd_tune *wf;
  double q,r,ctheta=1.0,Rx,Ry,Rz,Rinv,offset,phi,stheta;
	double x,y,z,qx,qy,qz,root2=sqrt(2.0);
  int iq,Nq;
  vector<double> cf;
  int imc,NMC;
  Crandy randy(-time(NULL));  // random number generator (in msu_commonutils)
  
  wf=new CWaveFunction_pd_tune("parameters/wfparameters.dat"); // This is for pp, solves Schrod. eq. 
  // parameter file includes information for solving and storing wave functions
  
	
	// Data PRC 82, 034002 (2010) singlet/triplet
	//delta(KE_proton=2.25 MeV) = -39.1, -34.5.    q= 58.11 MeV/c
	//delta(KE_proton=3.15 MeV) = -48.7, -42.9     q= 68.76 MeV/c

  //CparameterMap parmap;
  //parmap.ReadParsFromFile("parameters/coralpars.dat");  // Reads
  
	return 0;
}

