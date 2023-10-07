#include "msu_coral/wavefunction.h"
#include "msu_commonutils/constants.h"
#include "msu_commonutils/sf.h"
#include "msu_commonutils/randy.h"
using namespace std;
using namespace NMSUPratt;

CWaveFunction_pHe3_sqwell::CWaveFunction_pHe3_sqwell(string parsfilename) : CWaveFunction(){
	bool TUNINGPHE3=false;
	//CLog::Info("Quick and dirty potentials to replicate low-E phase shifts for p-He3\n");
	// Interaction fit to phaseshifts from T.V. Daniels et al, PRC (2010).
	ParsInit(parsfilename);

	m1=ProtonMass; 
	m2=2.0*ProtonMass+NeutronMass-7.718;
	IDENTICAL=0;

	q1q2=2;
	nchannels=4;

	ellmax=1;
	InitArrays();
	CLog::Info("Arrays Initialized\n");

	ell[0]=0;
	ell[1]=0;
	ell[2]=1;
	ell[3]=1;
	InitWaves();

	nwells=new int[nchannels];
	nwells[0]=2;
	nwells[1]=2;
	nwells[2]=2;
	nwells[3]=2;

	SquareWell_MakeArrays();
	
	// L=0, S=0
	a[0][0]=1.80628; a[0][1]=2.05664;
	V0[0][0]=-2.71804; V0[0][1]=-114.483;
	
	// L=0, S=1
	a[1][0]=0.835825; a[1][1]=1.10527;
	V0[1][0]=-6.16757; V0[1][1]=-155.224;
	
	// L=1, S=0
	a[2][0]=1.10123; a[2][1]=6.55019;
	V0[2][0]=11.5845; V0[2][1]=-1.51854;
	
	// L=1, S=1
	a[3][0]=0.0956635; a[3][1]=5.23552;
	V0[3][0]=31.0429; V0[3][1]=-4.76129;
	
	SquareWell_Init();
	
	if(TUNINGPHE3){
	
		// If you are doing fitting (must use delq=1 MeV/c)
		Crandy *randy=new Crandy(-time(NULL));
		double error,besterror=1.0E99,delta1,delta2,delta3,delta4;
		double delta1target,delta2target,delta3target,delta4target;
		vector<vector<double>> abest,Vbest;
		int ichannel,itry,ntry=1000;
	
		abest.resize(nchannels);
		Vbest.resize(nchannels);
		for(ichannel=0;ichannel<nchannels;ichannel++){
			abest[ichannel].resize(nwells[ichannel]);
			Vbest[ichannel].resize(nwells[ichannel]);
		}
	
		for(ichannel=0;ichannel<nchannels;ichannel++){
			for(int ir=0;ir<nwells[ichannel];ir++){
				abest[ichannel][ir]=a[ichannel][ir];
				Vbest[ichannel][ir]=V0[ichannel][ir];
			}
		}
		
		ichannel=3;
		double dela=0.03;
		double delV=0.5;
		
		// These targets assume delq=2

		if(ichannel==0){
			delta1target=-39.1;
			delta2target=-48.7;
			delta3target=-56.3;
			delta4target=-67.8;
		}
		else if(ichannel==1){
			delta1target=-34.5;
			delta2target=-42.9;
			delta3target=-49.3;
			delta4target=-58.6;
		}
		else if(ichannel==2){ // These are for L=1, S=0
			delta1target=8.0;
			delta2target=13.4;
			delta3target=17.3;
			delta4target=21.2;
		}
		else if(ichannel==3){ // These are averaged over different J for L=1,S=1
			delta1target=15.4;
			delta2target=25.5;
			delta3target=34.1;
			delta4target=46.0;
		}
		
		/* these were actual phase shifts for J=0,1,2, L=1, S=1 (before averaging)
	
		else if(ichannel==4){
		delta1target=5;
		delta2target=9.7;
		delta3target=14.1;
		delta4target=21.3;
		}
		else if(ichannel==4){
		delta1target=17.0;
		delta2target=27.0;
		delta3target=34.9;
		delta4target=45.2;
		}
		else if(ichannel==5){
		delta1target=16.5;
		delta2target=27.7;
		delta3target=37.6;
		delta4target=51.5;
		}*/
			
		// Data PRC 82, 034002 (2010) singlet/triplet
		//delta(KE_proton=2.25 MeV) = -39.1, -34.5.    q=48.706 //32.49 MeV/c
		//delta(KE_proton=3.15 MeV) = -48.7, -42.9     q=57.630 //38.44 MeV/c
		//delta(KE_proton=4.00 MeV) = -56.3, -49.3.    q=69.941 // 43.32 MeV/c
		//delta(KE_proton=5.55 MeV) = -67.8, -58.6     q=76.496 // 51.02 MeV/c
	
		int nmiss=0;
		for(itry=0;itry<ntry;itry++){
			nmiss+=1;
			if(itry>0){
				a[ichannel][0]=fabs(abest[ichannel][0]+dela*randy->ran_gauss());
				a[ichannel][1]=fabs(abest[ichannel][1]+dela*randy->ran_gauss());
				V0[ichannel][0]=Vbest[ichannel][0]+delV*randy->ran_gauss();
				V0[ichannel][1]=Vbest[ichannel][1]+delV*randy->ran_gauss();
			}
			SquareWell_Init();
			delta1=0.794*delta[ichannel][48]+0.206*delta[ichannel][49];
			delta2=0.870*delta[ichannel][57]+0.130*delta[ichannel][58];
			delta3=0.559*delta[ichannel][69]+0.441*delta[ichannel][70];
			delta4=0.004*delta[ichannel][75]+0.996*delta[ichannel][76];
			delta1*=180.0/PI;
			delta2*=180.0/PI;
			delta3*=180.0/PI;
			delta4*=180.0/PI;
			if((ichannel==0 || ichannel==1) && delta1>0)
				delta1-=180.0;
			if((ichannel==0 || ichannel==1) && delta2>0)
				delta2-=180.0;
			if((ichannel==0 || ichannel==1) && delta3>0)
				delta3-=180.0;
			if((ichannel==0 || ichannel==1) && delta4>0)
				delta4-=180.0;
			error=pow(delta1-delta1target,2)+pow(delta2-delta2target,2)
				+pow(delta3-delta3target,2)+pow(delta4-delta4target,2);
			//printf("itry=%d, error=%g, besterror=%g, V[1]=%g\n",itry,error,besterror,V0[ichannel][1]);
			if(error<besterror){
				abest[ichannel][0]=a[ichannel][0]; abest[ichannel][1]=a[ichannel][1];
				Vbest[ichannel][0]=V0[ichannel][0]; Vbest[ichannel][1]=V0[ichannel][1];
				printf("----itry=%d\n",itry);
				printf("a[%d][0]=%g; a[%d][1]=%g;\n",
				ichannel,a[ichannel][0],ichannel,a[ichannel][1]);
				printf("V0[%d][0]=%g; V0[%d][1]=%g;\n",
				ichannel,V0[ichannel][0],ichannel,V0[ichannel][1]);
				//CountNodes();
				besterror=error;
				printf("besterror=%g, delta1=%g =?%g, delta2=%g=?%g, delta3=%g=?%g, delta4=%g=?%g,    dela=%g, delV=%g\n",
				besterror,delta1,delta1target,delta2,delta2target,delta3,delta3target,delta4,delta4target,dela,delV);

			}
			if(nmiss>50){
				dela*=0.8;
				delV*=0.8;
				nmiss=0;
			}
		}
	
		for(ichannel=0;ichannel<nchannels;ichannel++){
			for(int ir=0;ir<nwells[ichannel];ir++){
				a[ichannel][ir]=abest[ichannel][ir];
				V0[ichannel][ir]=Vbest[ichannel][ir];
			}
		}
		//PrintPhaseShifts();
		//
	}
	if(TUNINGPHE3)
		exit(1);
}

CWaveFunction_pHe3_sqwell::~CWaveFunction_pHe3_sqwell(){
	SquareWell_DeleteArrays();
}


double CWaveFunction_pHe3_sqwell::CalcPsiSquared(int iq,double r,double ctheta){
	double psisquared,theta=acos(ctheta),x;
	double q=GetQ(iq);
	complex<double> psi,psia,Xlm00,Xlm10;
	psia=planewave[iq]->planewave(r,ctheta);
	complex<double> DelPhi[nchannels];

	SquareWell_GetDelPhi(iq,r,DelPhi);
	x=q*r/HBARC;
	Xlm00=sqrt(4.0*PI)*SpherHarmonics::Ylm(0,0,theta,0.0)/x;
	Xlm10=sqrt(4.0*PI/3.0)*SpherHarmonics::Ylm(1,0,theta,0.0)/x;
	// for S=0
	psi=psia+Xlm00*DelPhi[0]+Xlm10*3.0*ci*DelPhi[2];
	psisquared=real(psi*conj(psi))/4.0;
	// for S=1
	psi=psia+Xlm00*DelPhi[1]+Xlm10*3.0*ci*DelPhi[3];
	psisquared+=3.0*real(psi*conj(psi))/4.0;
	return psisquared;

}

void CWaveFunction_pHe3_sqwell::CountNodes(){
	int Nnodes=0;
	int ichannel,iq,ir;
	double sign=1.0;
	iq=0;
	ichannel=1;
	for(ir=1;ir<39;ir++){
		sign=real(DelPhiArray[iq][ir][ichannel]*DelPhiArray[iq][ir-1][ichannel]);
		if(sign<0.0){
			printf("ir=%d, sign=%g\n",ir,sign);
			Nnodes+=1;
		}
	}
	if(Nnodes>0)
		CLog::Info("----------- Nnodes="+to_string(Nnodes)+" -------------\n");
}
