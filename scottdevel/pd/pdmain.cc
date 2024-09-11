#include "msu_coral/coral.h"
#include <vector>
using namespace std;
using namespace NMSUPratt;

int main(){
	CWaveFunction_pd_sqwell *wf;
	double q,r,ctheta=1.0,Rx,Ry,Rz,Rinv=2.0,offset=0.0,phi,stheta;
	double x,y,z,qx,qy,qz,root2=sqrt(2.0);
	int iq,Nq;
	vector<double> cf;
	int imc,NMC=4000000;
	Crandy randy(-time(NULL));  // random number generator (in msu_commonutils)
  
	wf=new CWaveFunction_pd_sqwell("parameters/wfparameters.dat"); // This is for pp, solves Schrod. eq. 
	// parameter file includes information for solving and storing wave functions
	
	// Data PRC 82, 034002 (2010) singlet/triplet
	//delta(KE_proton=2.25 MeV) = -39.1, -34.5.    q= 58.11 MeV/c
	//delta(KE_proton=3.15 MeV) = -48.7, -42.9     q= 68.76 MeV/c

	//CparameterMap parmap;
	//parmap.ReadParsFromFile("parameters/coralpars.dat");  // Reads
	
	//wf->PrintPhaseShifts();

	printf("Enter Rinv: ");
	scanf("%lf",&Rinv);
	
	
	Rx=Ry=Rz=Rinv;
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
			/*do{
			x=Rx*root2*randy.ran_gauss();
			y=Ry*root2*randy.ran_gauss();
			z=Rz*root2*randy.ran_gauss();
			x=x+offset;
			r=sqrt(x*x+y*y+z*z);
		}while(r>40);*/
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
