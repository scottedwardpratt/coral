#include "wavefunction.h"

// for non-identical particles: square will with complex optical potential

CWaveFunction_optical::CWaveFunction_optical(string parsfilename) : CWaveFunction(){
	ParsInit(parsfilename);
	llmax.resize(nqmax);
	kmax=200;
	etavec.resize(nqmax);
	CL.resize(nqmax);
	a.resize(nqmax);

	A.resize(nqmax);
	B.resize(nqmax);
	qinside.resize(nqmax);
	
	int q1q2set,m1set,m2set,VRset,VIset,Rset;
	q1q2set=parameters.getI("Q1Q2",-1);
	m1set=parameters.getD("M1",938.272);
	m2set=parameters.getD("M2",938.272);
	VRset=parameters.getD("OPTICAL_VR",-8.26);
	VIset=parameters.getD("OPTICAL_VI",25.99);
	Rset=parameters.getD("OPTICAL_A",1.892);
	Reset(q1q2set,m1set,m2set,VRset,VIset,Rset);
  InitArrays();
  InitWaves();
}

CWaveFunction_optical::~CWaveFunction_optical(){
	int iq,ll;
	for(iq=0;iq<nqmax;iq++){
		for(ll=0;ll<llmax[iq];ll++){
			a[iq][ll].resize(0);
		}
		A[iq].resize(0);
		B[iq].resize(0);
		a[iq].resize(0);
		CL[iq].resize(0);
	}
	A.resize(0);
	B.resize(0);
	a.resize(0);
	CL.resize(0);
	llmax.resize(0);
}

void CWaveFunction_optical::Reset(int q1q2set,double m1set,double m2set,double VRset,double VIset,double Rset){
	int iq,L,ll;
	double K,qmag,phiq,delK,qmagR;
	complex<double> q2;
	m1=m1set;
	m2=m2set;
	mu=m1*m2/(m1+m2);
	q1q2=q1q2set;
	VR=VRset;
	VI=VIset;
	R=Rset;
	
	for(iq=0;iq<nqmax;iq++){
		qmagR=delq*iq*R/HBARC;
		if(qmagR>2.5){
			llmax[iq]=4+lrint(qmagR-2);
		}
		else
			llmax[iq]=4;
	}
	
	for(iq=0;iq<nqmax;iq++){
		if(a[iq].size()>0){
			for(L=0;L<=llmax[iq];L++){
				a[iq][L].resize(0);
			}
		}
		a[iq].resize(llmax[iq]+1);
		CL[iq].resize(llmax[iq]+1);
		for(L=0;L<=llmax[iq];L++){
			a[iq][L].resize(kmax);	
		}		
	}

	if(q1q2!=0)
		COULOMB=true;
	else
		COULOMB=false;
	nchannels=0;
	delK=delq;
	for(iq=0;iq<nqmax;iq++){
		K=(iq+0.5)*delK;
		q2=K*K-2.0*mu*(VR-ci*VI);
		qmag=sqrt(real(q2)+imag(q2));
		phiq=atan2(imag(q2),real(q2));
		qinside[iq]=qmag*cos(phiq)+ci*qmag*sin(phiq);
		etavec[iq]=ALPHA*q1q2*mu/qinside[iq];
	}
	
	GetCL();
	GetExpansionCoefficients();
	
	for(iq=0;iq<nqmax;iq++){
		A[iq].resize(llmax[iq]+1);
		B[iq].resize(llmax[iq]+1);
		for(ll=0;ll<=llmax[iq];ll++){
			A[iq][ll]=B[iq][ll]=0.0;
		}
	}
	GetAB();
}

void CWaveFunction_optical::ClearInfo(){
	int iq,L;
	for(iq=0;iq<nqmax;iq++){
		for(L=0;L<=llmax[iq];L++){
			a[iq][L].clear();
		}
		a[iq].clear();
		CL[iq].clear();
	}
	qinside.clear();
	a.clear();
	CL.clear();	
}

void CWaveFunction_optical::GetCL(){
	int iq,ll;
	double prefactor,cmag,phi;
	complex<double> CL2,gamow,etaq;
	for(iq=0;iq<int(etavec.size());iq++){
		etaq=etavec[iq];
		gamow=CoulWave::cgamma(1.0+ci*etaq);
		gamow*=exp(PI*etaq);
		CL2=gamow;
		prefactor=1.0;
		for(ll=0;ll<=llmax[iq];ll++){
			if(ll>0){
				prefactor*=2.0/((2.0*ll)*(2.0*ll-1));
				CL2*=(double(ll*ll)+etaq*etaq);
			}
			cmag=sqrt(real(CL2*conj(CL2)));
			phi=0.5*atan2(imag(CL2),real(CL2));
			CL[iq][ll]=prefactor*sqrt(cmag)*(cos(phi)+ci*sin(phi));
		}
	}
}

void CWaveFunction_optical::GetExpansionCoefficients(){
	int k,ll,iq;
	complex<double> etaq,term1,term2,term3;
	for(iq=0;iq<int(etavec.size());iq++){
		etaq=etavec[iq];
		for(ll=0;ll<=llmax[iq];ll++){
			a[iq][ll][0]=1.0;
			a[iq][ll][1]=0.0;
			a[iq][ll][2]=(ll+1.0-2.0*etaq*etaq)/(4.0*(ll+1.0)*(ll+1.0)*(2.0*ll+3.0));
			a[iq][ll][3]=etaq*((ll+1.0)*(ll+1.0)+etaq*etaq)/(3.0*(ll+1.0)*(ll+1.0)*(ll+1.0)*(ll+2.0)*(2.0*ll+3.0));
			for(k=3;k<kmax;k++){
				term1=etaq*(2.0*k*a[iq][ll][k] -a[iq][ll][k-2]/(ll+1.0));
				term2=a[iq][ll][k-1]*(2.0*etaq*etaq-(2.0*k-1.0)*(ll+1.0)) /(2.0*(ll+1.0));
				term3=a[iq][ll][k-3]/(4.0*(ll+1.0));
				a[iq][ll][k+1]=-(term1+term2+term3)/((ll+1.0)*(k+1.0)*(k+2.0*(ll+1.0)));
			}
		}
	}
}

void CWaveFunction_optical::GetF_Complex(int iq,int L,complex<double> rho,complex<double> &F,complex<double> &Fprime){
	int k;
	complex<double> Phi,Phiprime,sum=1.0,sumprime=0.0,dsum=1.0,dsumprime,etaq,q;
	etaq=etavec[iq];
	Phi=exp( ((etaq*rho)/(L+1.0))-rho*rho/(4.0*(L+1.0)) );
	Phiprime=Phi*( (etaq/(L+1.0))-2.0*rho/(4.0*(L+1.0)) );
	/* printf("------ iq=%d, L=%d ------\n",iq,L);
	printf("rho=(%g,%g), Phi=(%g,%g), Phiprime=(%g,%g)\n",real(rho),imag(rho),
	real(Phi),imag(Phi),real(Phiprime),imag(Phiprime)); */
	k=2;
	while(k<kmax && (k<100 || real(dsum*conj(dsum))>1.0E-40)){
		//printf("a[iq=%d][L=%d][k=%d]=(%g,%g)\n",iq,L,k,real(a[iq][L][k]),imag(a[iq][L][k]));
		dsum=a[iq][L][k]*pow(rho,k);
		dsumprime=a[iq][L][k]*double(k)*pow(rho,k-1);
		sum+=dsum;
		sumprime+=dsumprime;
		k+=1;
	}
	F=Phi*sum*pow(rho,L+1)*CL[iq][L];
	Fprime+=(Phiprime*sum+Phi*sumprime)*pow(rho,L+1)*CL[iq][L];
	Fprime+=F*(L+1.0)/rho;
}

void CWaveFunction_optical::GetAB(){
	int iq,ll;
	double etak,kR,K,qmag2,qmagR,phi;
	vector<double> FL,GL,FLprime,GLprime;
	complex<double> qR;
	complex<double> FLinside,FLprimeinside,HL,HLprime,HLstar,HLprimestar,numer,denom;
	for(iq=0;iq<nqmax;iq++){
		if(int(FL.size())!=llmax[iq]+1){
			FL.resize(llmax[iq]+1);
			GL.resize(llmax[iq]+1);
			FLprime.resize(llmax[iq]+1);
			GLprime.resize(llmax[iq]+1);
		}
		K=(iq+0.5)*delq;
		kR=K*R/HBARC;
		etak=ALPHA*q1q2*mu/K;
		qmag2=2.0*mu*VI+fabs(K*K-2.0*mu*VR);
		qinside[iq]=sqrt(qmag2);
		qmagR=sqrt(qmag2)*R/HBARC;
		phi=0.5*atan2(2.0*mu*VI,K*K-2.0*mu*VR);
		qinside[iq]*=exp(ci*phi);
		qR=qmagR*cos(phi)+ci*qmagR*sin(phi);
		printf("iq=%d: K=%g, kR=%g, qR=(%g,%g)\n",iq,K,kR,real(qR),imag(qR));
		
		CoulWave::GetFGPrimeArray(llmax[iq],kR,etak,FL,GL,FLprime,GLprime);
		
		for(ll=0;ll<=llmax[iq];ll++){
			GetF_Complex(iq,ll,qR,FLinside,FLprimeinside);
			if(ll<1){
				printf("--- ll=%d ----\n",ll);
				printf("FL=(%g,%g), GL=(%g,%g), FLprime=(%g,%g), GLprime=(%g,%g)\n",real(FL[ll]),imag(FL[ll]),
				real(GL[ll]),imag(GL[ll]),real(FLprime[ll]),imag(FLprime[ll]),real(GLprime[ll]),imag(GLprime[ll]));
				printf("inside FL-FLprime=(%g,%g), (%g,%g)\n",real(FLinside),imag(FLinside),
				real(FLprimeinside),imag(FLprimeinside));
			}
			HL=FL[ll]-ci*GL[ll];
			HLstar=conj(HL);
			HLprime=FLprime[ll]-ci*GLprime[ll];
			HLprimestar=conj(HLprime);
			numer=kR*FLinside*HLprime-qR*FLprimeinside*HL;
			denom=qR*FLprimeinside*HLstar-kR*FLinside*HLprimestar;
			B[iq][ll]=numer/denom;
			A[iq][ll]=0.5*(HL+B[iq][ll]*HLstar)/FLinside;
			if(ll<1){
				printf("ll=%d: checking (%g,%g)=?(%g,%g)\n",ll,real(A[iq][ll]*FLinside),imag(A[iq][ll]*FLinside),
				0.5*real(HL+B[iq][ll]*conj(HL)),0.5*imag(HL+B[iq][ll]*conj(HL)));
				printf("ll=%d: checking (%g,%g)=?(%g,%g)\n",ll,real(A[iq][ll]*FLprimeinside),imag(A[iq][ll]*FLprimeinside),
				0.5*real(HLprime+B[iq][ll]*conj(HLprime)),0.5*imag(HLprime+B[iq][ll]*conj(HLprime)));
				printf("checking |A|^2=%g\n",real(A[iq][ll]*conj(A[iq][ll])));
				printf("checking |B|^2=%g\n",real(B[iq][ll]*conj(B[iq][ll])));
				if(ll==0){
					printf("K=%g\n",K);
					printf("A[%d][L=%d]=(%g,%g)\n",iq,ll,real(A[iq][ll]),imag(A[iq][ll]));
					printf("B[%d][L=%d]=(%g,%g)\n",iq,ll,real(B[iq][ll]),imag(B[iq][ll]));
				}
			}
		}
	}
}

double CWaveFunction_optical::CalcPsiSquared(int iq,double r,double ctheta){
	complex<double> psi,Yl,x;
	vector<double> FL,GL,FLprime,GLprime;
	complex<double> FLinside,FLprimeinside,HL,HLstar;
	complex<double> hL,hLprime;
	double Amag,K=(0.5+iq)*delq,x0,psisquared,etak;
	int ll;
	if(int(FL.size())!=llmax[iq]+1){
		FL.resize(llmax[iq]+1);
		GL.resize(llmax[iq]+1);
		FLprime.resize(llmax[iq]+1);
		GLprime.resize(llmax[iq]+1);
	}
	
	etak=ALPHA*q1q2*mu/K;
	x0=K*r/HBARC;
	psi=planewave[iq]->planewave(r,ctheta);
	
	CoulWave::GetFGPrimeArray(llmax[iq],x0,etak,FL,GL,FLprime,GLprime);
	if(r<R){
		x=qinside[iq]*r/HBARC;
		for(ll=0;ll<=llmax[iq];ll++){
			Amag=fabs(real(A[iq][ll]*conj(A[iq][ll])));
			if(Amag>1.0E-20 && Amag<10000.0){
				GetF_Complex(iq,ll,qinside[iq]*r/HBARC,FLinside,FLprimeinside);
				Yl=SpherHarmonics::legendre(ll,ctheta);
				psi+=(2.0*ll+1.0)*pow(ci,ll)*Yl*(A[iq][ll]*FLinside/x0-FL[ll]/x0);
			}
		}
	}
	else{
		for(ll=0;ll<=llmax[iq];ll++){
			HL=FL[ll]-ci*GL[ll];
			HLstar=conj(HL);
			Yl=SpherHarmonics::legendre(ll,ctheta);
			psi+=0.5*(2.0*ll+1.0)*pow(ci,ll)*Yl*(B[iq][ll]-1.0)*HLstar/x0;
		}
	}
	psisquared=real(psi*conj(psi));	
	return psisquared;

}
