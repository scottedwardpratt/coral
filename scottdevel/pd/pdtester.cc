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
	int imc,NMC=1000000;
	Crandy randy(-time(NULL));  // random number generator (in msu_commonutils)
  
	wf=new CWaveFunction_pd_sqwell("parameters/wfparameters.dat"); // This is for pp, solves Schrod. eq. 
	// parameter file includes information for solving and storing wave functions
	
	// Data PRC 82, 034002 (2010) singlet/triplet
	//delta(KE_proton=2.25 MeV) = -39.1, -34.5.    q= 58.11 MeV/c
	//delta(KE_proton=3.15 MeV) = -48.7, -42.9     q= 68.76 MeV/c

	//CparameterMap parmap;
	//parmap.ReadParsFromFile("parameters/coralpars.dat");  // Reads
	
	//wf->PrintPhaseShifts();

	
	printf("Enter iq:");
	scanf("%d",&iq);
	q=wf->GetQ(iq);
	printf("q=%g\n",q);
	int ichannel;
	double delr=0.1,Integral[6]={0.0},lambda=1000.0;;
	complex<double> DelPhi[6],DelPhiPrime[6];
	double DelPhi2[6];
	
	for(r=0.5*delr;r<15*lambda;r+=delr){ int ic=-1;
	//for(r=2.158;r<2.160;r+=0.0001){ int ic=1;
		//for(r=2.41;r<2.42;r+=0.0001){ int ic=1;
		//for(r=7.899;r<7.90;r+=0.0001){ int ic=1;
		//for(r=10.11;r<10.12;r+=0.0001){ int ic=3;
		wf->SquareWell_CalcDelPhi2(iq,r,DelPhi,DelPhiPrime,DelPhi2);
		for(ichannel=0;ichannel<6;ichannel++){
			Integral[ichannel]+=delr*DelPhi2[ichannel]*exp(-0.5*(r/lambda)*(r/lambda));
		}
		
		if(ic>=0){
			complex<double> dphiprime,dphia[6],dphiprimea[6],dphib[6],dphiprimeb[6];
			double dphi2a[6],dphi2b[6];
			double delta=0.0001;
			wf->SquareWell_CalcDelPhi2(iq,r-0.5*delta,dphia,dphiprimea,dphi2a);
			wf->SquareWell_CalcDelPhi2(iq,r+0.5*delta,dphib,dphiprimeb,dphi2b);
			dphiprime=(dphib[ic]-dphia[ic])*HBARC/delta;
			printf("%7.5f: (%8.5f,%8.5f)  (%8.5f,%8.5f) =? (%8.5f,%8.5f)\n",
			r,real(DelPhi[ic]),imag(DelPhi[ic]),real(DelPhiPrime[ic]),imag(DelPhiPrime[ic]),
			real(dphiprime),imag(dphiprime));
		}
	}

	for(ichannel=0;ichannel<6;ichannel++){
		double delk=5.0/HBARC,guess;
		guess=(wf->delta[ichannel][iq+1]-wf->delta[ichannel][iq-1])/(2.0*delk);
		printf("---- ichannel=%d: Integral=%8.5f, ddelta/dq=%8.5f=?%8.5f, ratio=%g\n",ichannel,Integral[ichannel],2.0*Integral[ichannel],
		guess,2.0*Integral[ichannel]/guess);
	}
	exit(1);
	
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
