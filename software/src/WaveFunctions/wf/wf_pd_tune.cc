#include "msu_coral/wavefunction.h"
#include "msu_commonutils/constants.h"
#include "msu_commonutils/sf.h"
#include "msu_commonutils/randy.h"
using namespace std;
using namespace NMSUPratt;

CWaveFunction_pd_tune::CWaveFunction_pd_tune(string parsfilename) : CWaveFunction(){
	CLog::Info("Beware! The p-d wavefunction was tuned to match phase shifts which were only measured for q<100.\n Also, the p-d system becomes inelastic (deuteron breaks up) above q=52 MeV/c\n So this treatment is pretty questionable for q>50!\n");

	// Interaction fit to phaseshifts from T.C. Black et al., PLB 471, p. 103-107 (1999).
  ParsInit(parsfilename);
	
	if(nqmax!=7)
		CLog::Fatal("nqmax="+to_string(nqmax)+" should be 7\n");

  m1=ProtonMass; 
  m2=ProtonMass+NeutronMass-2.224;
  IDENTICAL=0;

  q1q2=1;
  nchannels=6;
	
  ellmax=2;
  InitArrays();
  CLog::Info("Arrays Initialized\n");

  ell[0]=0;
	ell[1]=0;
	ell[2]=1;
	ell[3]=1;
	ell[4]=2;
	ell[5]=2;
  InitWaves();

  nwells=new int[nchannels];
  nwells[0]=3;
	nwells[1]=3;
	nwells[2]=2;
	nwells[3]=2;
	nwells[4]=2;
	nwells[5]=2;

  SquareWell_MakeArrays();
	
	a[0][0]=4.33842; a[0][1]=4.7281; a[0][2]=9.99999; 
	V0[0][0]=-32.7951; V0[0][1]=38.9506; V0[0][2]=-0.591542;
	
	a[1][0]=2.1591; a[1][1]=2.4132; a[1][2]=7.89978; 
	V0[1][0]=-30.4977; V0[1][1]=38.9727; V0[1][2]=1.10959;
	
	a[2][0]=1.74688; a[2][1]=11.8412; 
	V0[2][0]=-9.99994; V0[2][1]=0.131973; 
	
	a[3][0]=3.52139; a[3][1]=10.1152; 
	V0[3][0]=-9.63551; V0[3][1]=-0.588381;
	
	a[4][0]=3.97389; a[4][1]=13.0859; 
	V0[4][0]=-9.8424; V0[4][1]=-0.0903357;
	
	a[5][0]=1.25107; a[5][1]=11.5236; 
	V0[5][0]=-9.11892; V0[5][1]=0.323206; 
	
	/*
	
	int iq,ichannel,iwell;
	double phaseshift,chisquare,bestchisquare;
	char dumbo[120];
	vector<double> deltaexp,qexp;
	int itry,ntries;
	Crandy randy(123);
	
	double dela=0.001;
	double delV0=0.001; 
	
	ichannel=5;
	
	if(ichannel==0){
		qexp.resize(nqmax);
		deltaexp.resize(nqmax);
		FILE *fptr=fopen("phaseshifts/data/S/delta_2S_12.txt","r");
		for(iq=0;iq<nqmax;iq++){
			fscanf(fptr,"%lf %lf",&qexp[iq],&deltaexp[iq]);
		}
	}
	if(ichannel==1){
		qexp.resize(nqmax);
		deltaexp.resize(nqmax);
		FILE *fptr=fopen("phaseshifts/data/S/delta_4S_32.txt","r");
		for(iq=0;iq<nqmax;iq++){
			fscanf(fptr,"%lf %lf",&qexp[iq],&deltaexp[iq]);
		}
	}
	if(ichannel==2){
		double dumbphase;
		qexp.resize(nqmax);
		deltaexp.resize(nqmax);
		FILE *fptr=fopen("phaseshifts/data/P/deltabar.txt","r");
		fgets(dumbo,120,fptr);
		for(iq=0;iq<nqmax;iq++){
			fscanf(fptr,"%lf %lf %lf",&qexp[iq],&deltaexp[iq],&dumbphase);
		}
	}
	if(ichannel==3){
		double dumbphase;
		qexp.resize(nqmax);
		deltaexp.resize(nqmax);
		FILE *fptr=fopen("phaseshifts/data/P/deltabar.txt","r");
		fgets(dumbo,120,fptr);
		for(iq=0;iq<nqmax;iq++){
			fscanf(fptr,"%lf %lf %lf",&qexp[iq],&dumbphase,&deltaexp[iq]);
		}
	}
	if(ichannel==4){
		double dumbphase;
		qexp.resize(nqmax);
		deltaexp.resize(nqmax);
		FILE *fptr=fopen("phaseshifts/data/D/deltabar.txt","r");
		fgets(dumbo,120,fptr);
		for(iq=0;iq<nqmax;iq++){
			fscanf(fptr,"%lf %lf %lf",&qexp[iq],&deltaexp[iq],&dumbphase);
		}
	}
	if(ichannel==5){
		double dumbphase;
		qexp.resize(nqmax);
		deltaexp.resize(nqmax);
		FILE *fptr=fopen("phaseshifts/data/D/deltabar.txt","r");
		fgets(dumbo,120,fptr);
		for(iq=0;iq<nqmax;iq++){
			fscanf(fptr,"%lf %lf %lf",&qexp[iq],&dumbphase,&deltaexp[iq]);
		}
	}
  
	bestchisquare=1.0E99;
	vector<vector<double>> abest,V0best;
	abest.resize(nchannels);
	V0best.resize(nchannels);
	
	//for(ichannel=0;ichannel<nchannels;ichannel++){
		abest[ichannel].resize(nwells[ichannel]);
		V0best[ichannel].resize(nwells[ichannel]);
		for(iwell=0;iwell<nwells[ichannel];iwell++){
			abest[ichannel][iwell]=a[ichannel][iwell];
			V0best[ichannel][iwell]=V0[ichannel][iwell];
		}
		//}
	
	int nsuccess=0;
	double dely,a0,amax=15;
	
	ntries=10000;
	
	for(itry=0;itry<ntries;itry++){
		SquareWell_Init();
		chisquare=0.0;
		for(iq=0;iq<nqmax;iq++){
			phaseshift=(180.0/PI)*GetDELTA(ichannel,iq);
			if(phaseshift>0.0)
				phaseshift-=180.0;
			//printf("%6.3f %g =? %g\n",qarray[iq],phaseshift,deltaexp[iq]);
			dely=fabs(phaseshift-deltaexp[iq]);
			if(dely>90.0)
				dely-=180.0;
			chisquare+=dely*dely;
		}
		if(chisquare<bestchisquare){
			bestchisquare=chisquare;
			for(iwell=0;iwell<nwells[ichannel];iwell++){
				abest[ichannel][iwell]=a[ichannel][iwell];
				V0best[ichannel][iwell]=V0[ichannel][iwell];
			}
			delV0*=0.999;
			dela*=0.999;
			nsuccess+=1;
			printf("success!!!!! itry=%d, bestchisquare=%g\n",itry,bestchisquare);
		}
		for(iwell=0;iwell<nwells[ichannel];iwell++){
			a0=0.0;
			if(iwell>0)
				a0=a[ichannel][iwell-1];
			do{
				a[ichannel][iwell]=abest[ichannel][iwell]+dela*randy.ran_gauss();
			}while(a[ichannel][iwell] < a0 || a[ichannel][iwell]>amax);
			V0[ichannel][iwell]=V0best[ichannel][iwell]+delV0*randy.ran_gauss();
		}
		
	}
	
	for(iwell=0;iwell<nwells[ichannel];iwell++){
		a[ichannel][iwell]=abest[ichannel][iwell];
		V0[ichannel][iwell]=V0best[ichannel][iwell];
	}
	*/
	
	SquareWell_Init();
	
	/*
	for(iq=0;iq<nqmax;iq++){
		phaseshift=(180.0/PI)*GetDELTA(ichannel,iq);
		if(phaseshift<0.0)
			phaseshift+=180.0;
		printf("%6.3f %10.4f =? %10.4f\n",qarray[iq],phaseshift,deltaexp[iq]);
	}
	
	for(iwell=0;iwell<nwells[ichannel];iwell++){
		printf("a[%d][%d]=%g; ",ichannel,iwell,abest[ichannel][iwell]);
	}
	printf("\n");
	for(iwell=0;iwell<nwells[ichannel];iwell++){
		printf("V0[%d][%d]=%g; ",ichannel,iwell,V0best[ichannel][iwell]);
	}
	printf("\n");
	printf("---- nsuccess=%d -----\n",nsuccess);
	*/
	
}

CWaveFunction_pd_tune::~CWaveFunction_pd_tune(){
  SquareWell_DeleteArrays();
}


double CWaveFunction_pd_tune::CalcPsiSquared(int iq,double r,double ctheta){
  double psisquared,theta=acos(ctheta),x;
  double q=GetQ(iq);
  complex<double> psi,psiplane,Xlm00,Xlm10,Xlm20;
  psiplane=planewave[iq]->planewave(r,ctheta);
  complex<double> DelPhi[6];

  SquareWell_GetDelPhi(iq,r,DelPhi);
	
	x=q*r/HBARC;
  Xlm00=sqrt(4.0*PI)*SpherHarmonics::Ylm(0,0,theta,0.0)/x;
	Xlm10=ci*sqrt(12.0*PI)*SpherHarmonics::Ylm(1,0,theta,0.0)/x;
	Xlm20=-sqrt(20.0*PI)*SpherHarmonics::Ylm(2,0,theta,0.0)/x;
	
	// for S=1/2
  psi=psiplane+Xlm00*DelPhi[0]+Xlm10*DelPhi[2]+Xlm20*DelPhi[4];
  psisquared=(1.0/3.0)*real(psi*conj(psi));
	// for S=3/2
  psi=psiplane+Xlm00*DelPhi[1]+Xlm10*DelPhi[3]+Xlm20*DelPhi[5];
  psisquared+=(2.0/3.0)*real(psi*conj(psi));
	
  return psisquared;

}