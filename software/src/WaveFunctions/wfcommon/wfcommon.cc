#ifndef __CWAVEFUNCTION_WFCOMMON_CC__
#define __CWAVEFUNCTION_WFCOMMON_CC__

#include "msu_coral/wavefunction.h"
#include "msu_commonutils/constants.h"
#include "msu_commonutils/randy.h"
#include "msu_commonutils/parametermap.h"
#include "msu_commonutils/sf.h"
#include "msu_commonutils/misc.h"
#include "msu_commonutils/log.h"

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <complex>
#include <cstring>

using namespace std;

Crandy *CWaveFunction::randy=NULL;

CWaveFunction::CWaveFunction(){
  generic=false;
  ci=complex<double>(0.0,1.0);
	randy=new Crandy(-12345);
  MPI=139.58;
  MKAON=493.677;
  MPROTON=938.271;
  MNEUTRON=939.565;
  MLAMBDA=1115.7;
	qarray=NULL;
}

CWaveFunction::~CWaveFunction(){
	int ichannel,iq,ell0;
	if(nchannels!=0){
		for(ichannel=0;ichannel<nchannels;ichannel++){
			delete Wepsilon[ichannel];
			delete delta[ichannel];
			delete ddeltadq[ichannel];
		}
		delete [] Wepsilon;
		delete [] delta;
		delete [] ddeltadq;
		delete [] qarray;
		delete [] eta;
		delete [] channelweight;
 
 
		for(iq=0;iq<nqmax;iq++)
			delete (planewave[iq]);
		delete [] planewave;
		for(ell0=0;ell0<=ellmax;ell0++){
			for(iq=0;iq<nqmax;iq++){
				if(partwave[ell0][iq]!=NULL) delete (partwave[ell0][iq]);
			}
			delete [] partwave[ichannel];
		}
		delete [] partwave;
	}
	
}

bool CWaveFunction::GetIDENTICAL(){
  return IDENTICAL;
}

int CWaveFunction::GetNQMAX(){
  return nqmax;
}

int CWaveFunction::GetNCHANNELS(){
  return nchannels;
}

double CWaveFunction::GetDELTA(int ichannel,int iq){
  return delta[ichannel][iq];
}

double CWaveFunction::GetDELQ(){
  return delq;
}

double CWaveFunction::GetQ(int iq){
  return qarray[iq];
}

void CWaveFunction::ParsInit(string parsfilename){
  FILE *qarrayfile;
  char qarrayfilename[120];
  string stemp;
  int iq;
  bool filetest=0;

  parameters.ReadParsFromFile(parsfilename);
  
  stemp=parameters.getS("QARRAYFILENAME","no_qarray_file");
  strcpy(qarrayfilename,stemp.c_str());
  delq=parameters.getD("DELQ",-999);
  nqmax=parameters.getI("NQMAX",-999);
  epsilon=parameters.getD("EPSILON",-999);
  COULOMB=parameters.getB("COULOMB",-999);
  STRONG=parameters.getB("STRONG",-999);
  IDENTICAL=parameters.getB("IDENTICAL",0);
	
  // If delq<0, read qarray from file (don't use this if for wf in kernels)
  if(delq<0) filetest=1;
  if(filetest==1){
    parameters.set("delq",-1.0);
    snprintf(message,CLog::CHARLENGTH,"will read qarray from %s\n",qarrayfilename);
	 CLog::Info(message);
    qarrayfile=fopen(qarrayfilename,"r");
    fscanf(qarrayfile,"%d",&nqmax);
    parameters.set("NQMAX",nqmax);
    qarray=new double[nqmax];
    for(iq=0;iq<nqmax;iq++){
      fscanf(qarrayfile,"%lf",&qarray[iq]);
    }
    fclose(qarrayfile);
  }
  else{
    delq=parameters.getD("DELQ",-999);
    qarray=new double[nqmax];
    for(iq=0;iq<nqmax;iq++){
      qarray[iq]=(iq+0.5)*delq;
    }
  }
}

void CWaveFunction::InitArrays(){
	int iq,l,ichannel;
	eta=new double[nqmax];
	planewave=new CPlaneWave*[nqmax];
	if(nchannels>0){
		partwave=new CPartWave **[ellmax+1];
		for(l=0;l<=ellmax;l++){
			partwave[l]=new CPartWave *[nqmax];
			for(iq=0;iq<nqmax;iq++)
				partwave[l][iq]=NULL;
		}

		delta=new double *[nchannels];
		ddeltadq=new double *[nchannels];
		Wepsilon=new double *[nchannels];
		for(ichannel=0;ichannel<nchannels;ichannel++){
			delta[ichannel]=new double[nqmax];
			ddeltadq[ichannel]=new double[nqmax];
			Wepsilon[ichannel]=new double[nqmax];
		}

		ell=new int[nchannels];
		channelweight=new double[nchannels];
	}
  
}

void CWaveFunction::InitWaves(){
  int iq,ichannel;
  double q,e1,e2,e;
	
  for(iq=0;iq<nqmax;iq++){
    q=qarray[iq];
    e1=sqrt(m1*m1+q*q);
    e2=sqrt(m2*m2+q*q);
    e=e1+e2;
    //eta[iq]=double(q1q2)*e1*e2/((e1+e2)*137.036*q); //old way (needs correction)
    eta[iq]=double(q1q2)*(pow(e,4)-pow(m1*m1-m2*m2,2))/(4.0*e*e*e*137.036*q);
		
    planewave[iq]=new CPlaneWave(eta[iq],q1q2,q);
  }
  for(ichannel=0;ichannel<nchannels;ichannel++){
    for(iq=0;iq<nqmax;iq++){
      q=qarray[iq];
      if(partwave[ell[ichannel]][iq]==NULL){
				partwave[ell[ichannel]][iq]
				=new CPartWave(eta[iq],q1q2,q,ell[ichannel],epsilon);
      }
    }
  }
}

double CWaveFunction::GetPsiSquared(double q,double r,double ctheta){
  int iq,iqlow,iqhigh;
  double wlow,whigh,interpolate,qscaled,rscaled;
	
	if(r<1000.0){
		if(generic==1 && q1q2!=0){
			qscaled=q*(muscale/mu)*q1q2scale/double(q1q2);
			rscaled=q*r/qscaled;
		}
		else{
			qscaled=q;
			rscaled=r;
		}
		
		if(delq<0){
			iq=0;
			while(qscaled>qarray[iq]+1.0E-5 && iq<nqmax){
				iq+=1;
			}
			iqlow=iq-1;
			iqhigh=iq;
		}
		else{
			iqhigh=lrint(qscaled/delq);
			if(iqhigh==0) iqhigh=1;
			iqlow=iqhigh-1;
		}
		if(iqhigh==nqmax && nqmax>1 && (qscaled-qarray[nqmax-1])<0.5*(qarray[nqmax-1]-qarray[nqmax-2])){
			iqhigh=nqmax-1;
			iqlow=iqhigh-1;
		}
		
		
		if(iqhigh==0)
			return CalcPsiSquared(0,rscaled,ctheta);
		else if(iqhigh>=nqmax && qscaled-qarray[nqmax-1]>1.0E-5)
			return 1.0;
		else{
			wlow=(qarray[iqhigh]-qscaled)/(qarray[iqhigh]-qarray[iqlow]);
			whigh=1.0-wlow;
			if(fabs(wlow)<1.0E-5)
				interpolate=CalcPsiSquared(iqhigh,rscaled,ctheta);
			else if(fabs(whigh)<1.0E-5) 
				interpolate=CalcPsiSquared(iqlow,rscaled,ctheta);
			else{
				interpolate=wlow*CalcPsiSquared(iqlow,rscaled,ctheta)
					+whigh*CalcPsiSquared(iqhigh,rscaled,ctheta);
			}
			if(rscaled>1.0 && qscaled>qarray[0] && qscaled<qarray[nqmax-1] && interpolate<-0.01){
				snprintf(message,CLog::CHARLENGTH,"interpolate=%g, qscaled=%g, qlow=%g, qhigh=%g, rscaled=%g, r=%g\n",
				interpolate,qscaled,qarray[iqlow],qarray[iqhigh],rscaled,r);
				CLog::Info(message);
				snprintf(message,CLog::CHARLENGTH,"wlow=%g, whigh=%g\n",wlow,whigh);
				CLog::Info(message);
				snprintf(message,CLog::CHARLENGTH,"iqlow=%d,  %g\n",iqlow,CalcPsiSquared(iqlow,rscaled,ctheta));
				CLog::Info(message);
				snprintf(message,CLog::CHARLENGTH,"iqhigh=%d, %g\n",iqhigh,CalcPsiSquared(iqhigh,rscaled,ctheta));
				CLog::Info(message);
				snprintf(message,CLog::CHARLENGTH,"nqmax=%d\n",nqmax);
				CLog::Info(message);
				for(iq=0;iq<nqmax;iq++){
					snprintf(message,CLog::CHARLENGTH,"qarray[%d]=%g\n",iq,qarray[iq]);
					CLog::Info(message);
				}
				CLog::Fatal("dying in CWaveFunction::GetPsiSquared\n");
			}
			return interpolate;
		}
	}
	else return 1.0;
}

double CWaveFunction::GetPsiSquared(double *pa,double *xa,double *pb,double *xb){
	double q,r,ctheta;
	getqrctheta(pa,xa,pb,xb,&q,&r,&ctheta);
	return GetPsiSquared(q,r,ctheta);
}

double CWaveFunction::CalcPsiSquared(int iq,double r,double ctheta){
	return 1.0+iq+r+ctheta;  // This is a dummy function to be overwritten by inherited class
}

void CWaveFunction::getqrctheta(double *p1,double *r1,double *p2,double *r2,double *q,double *r,double *ctheta){
	int alpha;
	const double g[4]={1.0,-1.0,-1.0,-1.0};
	double n[4],qvec[4],rvec[4],nnorm,ndotq,ndotr;
	
	nnorm=0.0;
	ndotq=0.0;
	ndotr=0.0;
	for(alpha=0;alpha<4;alpha++){
		n[alpha]=p1[alpha]+p2[alpha];
		qvec[alpha]=0.5*(p1[alpha]-p2[alpha]);
		rvec[alpha]=r1[alpha]-r2[alpha];
		nnorm+=g[alpha]*n[alpha]*n[alpha];
		ndotq+=n[alpha]*qvec[alpha]*g[alpha];
		ndotr+=n[alpha]*rvec[alpha]*g[alpha];
	}
	nnorm=sqrt(nnorm);
	ndotq=ndotq/nnorm;
	ndotr=ndotr/nnorm;
	
	*ctheta=0.0;
	*r=0.0;
	*q=0.0;
	for(alpha=0;alpha<4;alpha++){
		n[alpha]=n[alpha]/nnorm;
		rvec[alpha]=rvec[alpha]-ndotr*n[alpha];
		qvec[alpha]=qvec[alpha]-ndotq*n[alpha];
		*r-=g[alpha]*rvec[alpha]*rvec[alpha];
		*q-=g[alpha]*qvec[alpha]*qvec[alpha];
		*ctheta-=g[alpha]*rvec[alpha]*qvec[alpha];
	}
	if(*r<0.0 || *q<0.0 || fabs(*ctheta)/sqrt(*q**r)>1.0000001){
		snprintf(message,CLog::CHARLENGTH,"Disaster, r^2=%g, q^2=%g, ctheta=%g\n",*r,*q,*ctheta/sqrt(*r**q));
		CLog::Fatal(message);
	}
	*r=sqrt(*r);
	*q=sqrt(*q);
	if(*q**r>0.1)	*ctheta=*ctheta/(*r**q);
	else{
		*ctheta=(1.0-2.0*randy->ran());
	}
	if(*ctheta>1.0)*ctheta=1.0;
}

void CWaveFunction::getqrctheta(FourVector &p1,FourVector &r1,FourVector &p2,FourVector &r2,double &q,double &r,double &ctheta){
	int alpha;
	const double g[4]={1.0,-1.0,-1.0,-1.0};
	double n[4],qvec[4],rvec[4],nnorm,ndotq,ndotr;
	
	nnorm=0.0;
	ndotq=0.0;
	ndotr=0.0;
	for(alpha=0;alpha<4;alpha++){
		n[alpha]=p1[alpha]+p2[alpha];
		qvec[alpha]=0.5*(p1[alpha]-p2[alpha]);
		rvec[alpha]=r1[alpha]-r2[alpha];
		nnorm+=g[alpha]*n[alpha]*n[alpha];
		ndotq+=n[alpha]*qvec[alpha]*g[alpha];
		ndotr+=n[alpha]*rvec[alpha]*g[alpha];
	}
	nnorm=sqrt(nnorm);
	ndotq=ndotq/nnorm;
	ndotr=ndotr/nnorm;
	
	ctheta=0.0;
	r=0.0;
	q=0.0;
	for(alpha=0;alpha<4;alpha++){
		n[alpha]=n[alpha]/nnorm;
		rvec[alpha]=rvec[alpha]-ndotr*n[alpha];
		qvec[alpha]=qvec[alpha]-ndotq*n[alpha];
		r-=g[alpha]*rvec[alpha]*rvec[alpha];
		q-=g[alpha]*qvec[alpha]*qvec[alpha];
		ctheta-=g[alpha]*rvec[alpha]*qvec[alpha];
	}
	if(r<0.0 || q<0.0 || fabs(ctheta)/sqrt(q*r)>1.0000001){
		snprintf(message,CLog::CHARLENGTH,"Disaster, r^2=%g, q^2=%g, ctheta=%g\n",r,q,ctheta/sqrt(r*q));
		CLog::Fatal(message);
	}
	r=sqrt(r);
	q=sqrt(q);
	if(q*r>0.1)
		ctheta=ctheta/(r*q);
	else{
		ctheta=(1.0-2.0*randy->ran());
	}
	if(ctheta>1.0)
		ctheta=1.0;
}

void CWaveFunction::PrintPhaseShifts(){
	int iq,ichannel;
	double q,KE1,mu;
	printf("printing phaseshifts\n");
	snprintf(message,CLog::CHARLENGTH,"-------- PHASE SHIFTS --------\n");
	CLog::Info(message);
	snprintf(message,CLog::CHARLENGTH,"q(MeV/c)  KE_1(MeV)  ");
	CLog::Info(message);
	for(ichannel=0;ichannel<nchannels;ichannel++){
		snprintf(message,CLog::CHARLENGTH,"    l=%d   ",ell[ichannel]);
		CLog::Info(message);
	}
	snprintf(message,CLog::CHARLENGTH,"\n");
	CLog::Info(message);
	
	for(iq=0;iq<nqmax;iq++){
		q=GetQ(iq);
		mu=m1*m2/(m1+m2);
		KE1=0.5*m1*q*q/(mu*mu);
		snprintf(message,CLog::CHARLENGTH,"%7.2f %7.2f ",q,KE1);
		CLog::Info(message);
		for(ichannel=0;ichannel<nchannels;ichannel++){
			snprintf(message,CLog::CHARLENGTH,"% 10.3f",(180.0/PI)*delta[ichannel][iq]-180.0);
			CLog::Info(message);
		}
		CLog::Info("\n");
	}
}

void CWaveFunction::PrintCdelta(double Rx,double Ry,double Rz){
	double q,clocal;
	int ichannel,iq;
	snprintf(message,CLog::CHARLENGTH,"! Qinv  C(Q)_estimated ~ ddelta/dq\n");
	CLog::Info(message);
	for(iq=0;iq<nqmax;iq++){
		q=qarray[iq];
		clocal=1.0;
		for(ichannel=0;ichannel<nchannels;ichannel++){
			clocal+=(channelweight[ichannel]*(2.0*PI)*pow(HBARC,3)
				/(q*q*Rx*Ry*Rz*pow(4.0*PI,1.5)))
			*ddeltadq[ichannel][iq];
		}
		snprintf(message,CLog::CHARLENGTH,"%6.2f  %8.4f  %g\n",q,clocal,4.0*q*q*(clocal-1.0));
		CLog::Info(message); 
	}
	snprintf(message,CLog::CHARLENGTH,"_________________________________\n");
	CLog::Info(message);
}

double CWaveFunction::RelativisticCorrection(double r,int iq){
	if(q1q2==0) return 1.0;
	else{
		double q,E,E1,E2,dmudE;
		q=GetQ(iq);
		E1=sqrt(m1*m1+q*q);
		E2=sqrt(m2*m2+q*q);
		E=E1+E2;
		dmudE=0.25*(1.0-(3.0/pow(E,4))*pow(m1*m1-m2*m2,2));
		return 1.0-dmudE*(E/(E1*E2))*q1q2*HBARC*ALPHA/r;	
	}
}

void CWaveFunction::EffectiveRange(int ichannel,double scattlength,double Reff){
	int iq;
	double q,tandel;
	for(iq=0;iq<nqmax;iq++){
		q=qarray[iq];
		tandel=1.0/( (-HBARC/(q*scattlength))+0.5*q*Reff/HBARC);
		delta[ichannel][iq]=atan(tandel);
		ddeltadq[ichannel][iq]=(tandel*tandel/(1.0+tandel*tandel))
		*((-HBARC/(q*q*scattlength))-0.5*Reff/HBARC);
	}
}

double CWaveFunction::GetIW(int ellval,double epsval,double q,int q1q2val,double etaval,double deltaval){
	complex<double> psi0,psi,psi2,psiminus0,psiminus,psiminus2;
	complex<double> psiplus0,psiplus,psiplus2;
	complex<double> ddeta_psi,ddeta_psiminus,ddeta_psiplus;
	complex<double> I,e2idelta;
	double x,deleta,a2,root;
	
	if(q1q2val!=0){
		deleta=0.002*etaval;
		x=q*epsval/HBARC;
		e2idelta=Misc::ceiphi(2.0*deltaval);
		if(ellval==0){
			a2=1.0+etaval*etaval;
			root=sqrt(a2);
			
			psi0=CoulWave::CWoutgoing(ellval,x,etaval-0.5*deleta);
			psi0=psi0*e2idelta+conj(psi0);
			
			psiplus0=CoulWave::CWoutgoing(ellval+1,x,etaval-0.5*deleta);
			psiplus0=psiplus0*e2idelta+conj(psiplus0);
			
			psi=CoulWave::CWoutgoing(ellval,x,etaval);
			psi=psi*e2idelta+conj(psi);
			
			psiplus=CoulWave::CWoutgoing(ellval+1,x,etaval);
			psiplus=psiplus*e2idelta+conj(psiplus);
			
			psi2=CoulWave::CWoutgoing(ellval,x,etaval+0.5*deleta);
			psi2=psi2*e2idelta+conj(psi2);
			
			psiplus2=CoulWave::CWoutgoing(ellval+1,x,etaval+0.5*deleta);
			psiplus2=psiplus2*e2idelta+conj(psiplus2);
			
			ddeta_psi=(psi2-psi0)/deleta;
			ddeta_psiplus=(psiplus2-psiplus0)/deleta;
			
			I=(conj(psi)*psi+conj(psiplus)*psiplus)*x*a2;
			I-=conj(psi)*psiplus*(a2*(1.0+2.0*etaval*x)+etaval*etaval)/root;
			I+=(conj(psiplus)*ddeta_psi-conj(psi)*ddeta_psiplus)*etaval*root;
		}
		else{
			a2=double(ellval*ellval)+etaval*etaval;
			root=sqrt(a2);
			
			psi0=CoulWave::CWoutgoing(ellval,x,etaval-0.5*deleta);
			psi0=psi0*e2idelta+conj(psi0);
			
			psiminus0=CoulWave::CWoutgoing(ellval-1,x,etaval-0.5*deleta);
			psiminus0=psiminus0*e2idelta+conj(psiminus0);
			
			psi=CoulWave::CWoutgoing(ellval,x,etaval);
			psi=psi*e2idelta+conj(psi);
			
			psiminus=CoulWave::CWoutgoing(ellval-1,x,etaval);
			psiminus=psiminus*e2idelta+conj(psiminus);
			
			psi2=CoulWave::CWoutgoing(ellval,x,etaval+0.5*deleta);
			psi2=psi2*e2idelta+conj(psi2);
			
			psiminus2=CoulWave::CWoutgoing(ellval-1,x,etaval+0.5*deleta);
			psiminus2=psiminus2*e2idelta+conj(psiminus2);
			
			ddeta_psi=(psi2-psi0)/deleta;
			ddeta_psiminus=(psiminus2-psiminus0)/deleta;
			
			I=(conj(psi)*psi+conj(psiminus)*psiminus)*a2*x/double(ellval*ellval);
			I-=conj(psiminus)*psi
			*((2.0*ellval+1.0)*a2*double(ellval)+2.0*etaval*x*a2
				-etaval*etaval*double(ellval))/(double(ellval*ellval)*root);
				I+=(conj(ddeta_psiminus)*psi-conj(ddeta_psi)*psiminus)
				*etaval*root/double(ellval);
				
		}
		I=0.25*I;
	}
	else{
		x=q*epsval/HBARC;
		if(ellval==0){
			psi=sin(x+deltaval);
			psiplus=-cos(x+deltaval)+sin(x+deltaval)/x;
			psi=0.5*x*(Bessel::hstarn(0,x)+Misc::ceiphi(2.0*deltaval)*Bessel::hn(0,x));
			psiplus=0.5*x*(Bessel::hstarn(1,x)
				+Misc::ceiphi(2.0*deltaval)*Bessel::hn(1,x));
			psi*=Misc::ceiphi(-deltaval); psiplus*=Misc::ceiphi(-deltaval);
			
			I=(conj(psi)*psi+conj(psiplus)*psiplus)*x;
			I-=conj(psi)*psiplus;
		}
		else{
			psi=x*Bessel::hn(ellval,x);
			psi=(Misc::ceiphi(2.0*deltaval)*psi+conj(psi))*0.5;
			psiminus=x*Bessel::hn(ellval-1,x);
			psiminus=0.5*(Misc::ceiphi(2.0*deltaval)*psiminus+conj(psiminus));
			
			I=(conj(psi)*psi+conj(psiminus)*psiminus)*x;
			I-=(2.0*double(ellval)+1)*conj(psiminus)*psi;
			
		}
	}
	return -real(I)/q;
}

#endif
