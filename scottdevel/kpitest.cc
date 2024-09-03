#include "msu_coral/coral.h"
#include <vector>
using namespace std;
using namespace NMSUPratt;

int main(){
  CWaveFunction_kpluspiminus_sqwell *wf;
  double q,r,ctheta=1.0,Rx,Ry,Rz,Rinv,offset,phi,stheta;
	double x,y,z,qx,qy,qz,root2=sqrt(2.0);
  int iq,Nq;
  vector<double> cf;
  int imc,NMC;
  Crandy randy(-time(NULL));  // random number generator (in msu_commonutils)
  
  wf=new CWaveFunction_kpluspiminus_sqwell("parameters/wfparameters.dat"); 

  CparameterMap parmap;  
  parmap.ReadParsFromFile("parameters/coralpars.dat");  // Reads
  NMC=parmap.getD("NMC",10000);  // Number of Monte Carlo points to sample CF
  Rinv=parmap.getD("RINV",3.0);
	Rx=Ry=Rz=Rinv;
  offset=parmap.getD("OFFSET",0.0);  // separation between two particles (along outward direction)
	

	Nq=wf->GetNQMAX(); // Number of q values
	cf.resize(Nq);
	for(iq=0;iq<Nq;iq++){
		q=wf->GetQ(iq);
		ctheta=1.0-2.0*randy.ran();
		stheta=sqrt(1.0-ctheta*ctheta);
		phi=2.0*PI*randy.ran();
		qz=q*ctheta;
		qx=q*stheta*cos(phi);
		qy=q*stheta*sin(phi);
		
		cf[iq]=0.0;
		for(imc=0;imc<NMC;imc++){
			x=Rx*root2*randy.ran_gauss();
			y=Ry*root2*randy.ran_gauss();
			z=Rz*root2*randy.ran_gauss();
			x=x+offset;
			r=sqrt(x*x+y*y+z*z);
			ctheta=(qx*x+qy*y+qz*z)/(q*r);
			//cf[iqdir][iq]+=wf->GetPsiSquared(q,r,ctheta);
			cf[iq]+=wf->CalcPsiSquared(iq,r,ctheta);
		}
		cf[iq]=cf[iq]/double(NMC);
	}

	char filename[120];
	snprintf(filename,120,"results/Rinv%g.txt",Rinv);
	FILE *fptr=fopen(filename,"w");
	printf(" q      C(qinv)\n");
	fprintf(fptr," q      C(qinv)\n");
	for(iq=0;iq<Nq;iq++){
		q=wf->GetQ(iq);
		printf("%5.1f %8.5f\n",q,cf[iq]);		
		fprintf(fptr,"%5.1f %8.5f\n",q,cf[iq]);		
	}
	fclose(fptr);
	return 0;
}

