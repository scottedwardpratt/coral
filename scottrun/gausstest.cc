#include "coral.h"
#include <vector>

using namespace std;

int main(){
  CWaveFunction *wf;
  double q,r,ctheta=1.0,Rx,Ry,Rz,offset;
	double x,y,z,x1,y1,z1,x2,y2,z2,qx,qy,qz,root2=sqrt(2.0);
  int iq,Nq,iqdir;
  vector<vector<double>> cf;
  int imc,NMC;
  CRandy randy(-time(NULL));
	string wf_parsfilename="parameters/wfparameters.dat";
  wf=new CWaveFunction_pp_schrod(wf_parsfilename);
	CparameterMap parmap;
	string parsfilename="coralpars.dat";
	parmap.ReadParsFromFile(parsfilename);
	NMC=parmap.getD("NMC",10000);
	Rz=parmap.getD("Rlong",4.0);
	Ry=parmap.getD("Rside",3.0);
	Rx=parmap.getD("Rout",3.5);
	offset=parmap.getD("OFFSET",0.0);  // separation between two particles (along outward direction)

	Nq=wf->GetNQMAX();
	cf.resize(3);
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
				x1=Rx*root2*randy.ran_gauss();
				y1=Ry*root2*randy.ran_gauss();
				z1=Rz*root2*randy.ran_gauss();
				x2=Rx*root2*randy.ran_gauss();
				y2=Ry*root2*randy.ran_gauss();
				z2=Rz*root2*randy.ran_gauss();
				x=x1-x2;
				y=y1-y2;
				z=z1-z2+offset;
				r=sqrt(x*x+y*y+z*z);
				ctheta=(qx*x+qy*y+qz*z)/(q*r);
				//c[iq]+=wf->GetPsiSquared(q,r,ctheta);
				cf[iqdir][iq]+=wf->CalcPsiSquared(iq,r,ctheta);
			}
			cf[iqdir][iq]=cf[iqdir][iq]/double(NMC);
		}
		printf("--- finished %d / 3 directions ---\n",iqdir+1);
	}
	char filename[120];
	sprintf(filename,"results/Rx%g_Ry%g_Rz%g_offset%g.dat",Rx,Ry,Rz,offset);
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

