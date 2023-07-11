#include "msu_coral/coral.h"
#include <vector>
using namespace std;

int main(){
  CWaveFunction *wf;
  double q,r,ctheta=1.0,Rx,Ry,Rz,offset;
	double x,y,z,qx,qy,qz,root2=sqrt(2.0);
  int iq,Nq,iqdir;
  vector<vector<double>> cf;
  int imc,NMC;
  Crandy randy(-time(NULL));  // random number generator (in msu_commonutils)
  
  wf=new CWaveFunction_pHe3_sqwell("parameters/wfparameters.dat"); // This is for pp, solves Schrod. eq. 
  // parameter file includes information for solving and storing wave functions
  
  CparameterMap parmap;  
  parmap.ReadParsFromFile("parameters/coralpars.dat");  // Reads
  NMC=parmap.getD("NMC",10000);  // Number of Monte Carlo points to sample CF
  Rz=parmap.getD("Rlong",4.0);
  Ry=parmap.getD("Rside",3.0);
  Rx=parmap.getD("Rout",3.5);
  offset=parmap.getD("OFFSET",0.0);  // separation between two particles (along outward direction)

  Nq=wf->GetNQMAX(); // Number of q values
  cf.resize(3); // Do out/long/side
  for(iqdir=0;iqdir<3;iqdir++){
	cf[iqdir].resize(Nq);
	for(iq=0;iq<Nq;iq++){
		q=wf->GetQ(iq);
		if(iqdir==0){
			qx=q; qy=qz=0.0;
		}
		else if(iqdir==1){
			qy=q; qx=qz=0.0;
		}
		else{
			qz=q; qx=qy=0.0;
		}
		cf[iqdir][iq]=0.0;
		for(imc=0;imc<NMC;imc++){
			x=Rx*root2*randy.ran_gauss();
			y=Ry*root2*randy.ran_gauss();
			z=Rz*root2*randy.ran_gauss();
			x=x+offset;
			r=sqrt(x*x+y*y+z*z);
			ctheta=(qx*x+qy*y+qz*z)/(q*r);
			//cf[iqdir][iq]+=wf->GetPsiSquared(q,r,ctheta);
			cf[iqdir][iq]+=wf->CalcPsiSquared(iq,r,ctheta);
		}
		cf[iqdir][iq]=cf[iqdir][iq]/double(NMC);
	}
	printf("--- finished %d / 3 directions ---\n",iqdir+1);
}
char filename[120];
snprintf(filename,120,"results/Rx%g_Ry%g_Rz%g_offset%g.dat",Rx,Ry,Rz,offset);
FILE *fptr=fopen(filename,"w");
printf(" q      Cout      Cside     Clong\n");
fprintf(fptr," q      Cout      Cside     Clong\n");
	for(iq=0;iq<Nq;iq++){
		q=wf->GetQ(iq);
		printf("%5.1f %8.5f  %8.5f  %8.5f\n",q,cf[0][iq],cf[1][iq],cf[2][iq]);		
		fprintf(fptr,"%5.1f %8.5f  %8.5f  %8.5f\n",
		q,cf[0][iq],cf[1][iq],cf[2][iq]);		
	}
	fclose(fptr);
	return 0;
}

