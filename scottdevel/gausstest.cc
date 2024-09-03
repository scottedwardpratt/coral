#include "msu_coral/coral.h"
#include <vector>
using namespace std;
using namespace NMSUPratt;

int main(){
	CWaveFunction *wf;
	double q,r,cthetaq,sthetaq,ctheta,phi,Rx,Ry,Rz,offset,wf2;
	double x,y,z,qx,qy,qz,root2=sqrt(2.0);
	int iq,Nq,idir,Ndir=5;
	vector<vector<double>> cf;
	int imc,NMC;
	Crandy randy(-time(NULL));  // random number generator (in msu_commonutils)
	string qdirection[5]={"qinv","qout","qside","qlong","ell=1_m=0(using x axis)"};
  
	wf=new CWaveFunction_pp_schrod("parameters/wfparameters.dat"); // This is for pp, solves Schrod. eq.
	//wf=new CWaveFunction_ppbar_nocoul("parameters/wfparameters.dat"); // This is for pp, solves Schrod. eq. 
	//wf=new CWaveFunction_ppiminus_sqwell("parameters/wfparameters.dat"); // This is for pp, solves Schrod. eq. 
	// parameter file includes information for solving and storing wave functions
	//wf=new CWaveFunction_pd_sqwell("parameters/wfparameters.dat");
  
	CparameterMap parmap;  
	parmap.ReadParsFromFile("parameters/coralpars.dat");  // Reads
	NMC=parmap.getD("NMC",10000);  // Number of Monte Carlo points to sample CF
	Rz=parmap.getD("Rlong",4.0);
	Ry=parmap.getD("Rside",3.0);
	Rx=parmap.getD("Rout",3.5);
	offset=parmap.getD("OFFSET",0.0);  // separation between two particles (along outward direction)
	//offset=0.0;
	//double Rinv=parmap.getD("RINV",3.0);
	//printf("Enter Rinv: ");
	//scanf("%lf",&Rinv);
	//Rx=Ry=Rz=Rinv;
	
	Nq=wf->GetNQMAX(); // Number of q values
	cf.resize(Ndir);
	for(idir=0;idir<Ndir;idir++){
		cf[idir].resize(Nq);
		for(iq=0;iq<Nq;iq++)
			cf[idir][iq]=0.0;
		
	}
	for(idir=0;idir<Ndir;idir++){
		printf("-------  %s --------\n",qdirection[idir].c_str());
		for(iq=0;iq<Nq;iq++){
			q=wf->GetQ(iq);
			for(imc=0;imc<NMC;imc++){
				if(idir==0 || idir==4){
					cthetaq=1.0-2.0*randy.ran();
					sthetaq=sqrt(1.0-cthetaq*cthetaq);
					phi=2.0*PI*randy.ran();
					qz=q*cthetaq;
					qx=q*sthetaq*cos(phi);
					qy=q*sthetaq*sin(phi);				
				}
				else if(idir==1){
					qx=q; qy=qz=0.0;
				}
				else if(idir==2){
					qy=q; qx=qz=0.0;
				}
				else if(idir==3){
					qz=q; qx=qy=0.0;
				}
		
			
				x=Rx*root2*randy.ran_gauss();
				y=Ry*root2*randy.ran_gauss();
				z=Rz*root2*randy.ran_gauss();
				x=x+offset;
				r=sqrt(x*x+y*y+z*z);
				ctheta=(qx*x+qy*y+qz*z)/(q*r);
				wf2=wf->CalcPsiSquared(iq,r,ctheta);
				if(idir!=4){
					cf[idir][iq]+=wf2;
				}
				else{
					cf[idir][iq]+=wf2*qx/q;
				}
			}
		}
		for(iq=0;iq<Nq;iq++)
			cf[idir][iq]=cf[idir][iq]/double(NMC);
	}


	char filename[120];
	snprintf(filename,120,"results/Rx%g_Ry%g_Rz%g_offset%g.txt",Rx,Ry,Rz,offset);
	FILE *fptr=fopen(filename,"w");
	fprintf(fptr,"#  q      inv      out      long     side    ell=1,m=0(along x axis)\n");
	printf("  q      inv      out      long     side    ell=1,m=0(along x axis)\n");
	for(iq=0;iq<Nq;iq++){
		q=wf->GetQ(iq);
		printf("%5.1f ",q);		
		fprintf(fptr,"%5.1f ",q);	
		for(idir=0;idir<Ndir;idir++){
			printf("%8.5f ",cf[idir][iq]);
			fprintf(fptr,"%8.5f ",cf[idir][iq]);
		}
		printf("\n");
		fprintf(fptr,"\n");
	}
	fclose(fptr);
	return 0;
}

