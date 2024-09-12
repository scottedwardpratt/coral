#include "msu_coral/wavefunction.h"
#include "msu_commonutils/constants.h"
#include "msu_commonutils/misc.h"
#include "msu_commonutils/sf.h"
#include "msu_commonutils/gslmatrix.h"
using namespace std;
using namespace NMSUPratt;

void CWaveFunction::SquareWell_Init(){
	// To set up the wave functions and phase shifts
	CGSLMatrix_Complex *cmatrix;
	double q,beta;
	double F1b,G1b,F1bprime,G1bprime;
	double F2a,G2a,F2aprime,G2aprime;
	double F2b,G2b,F2bprime,G2bprime;
	double F3a,G3a,F3aprime,G3aprime;
	double F3b,G3b,F3bprime,G3bprime;
	double F,G,Fprime,Gprime;
	double F1,G1,F1prime,G1prime;
	complex<double> *Avec;
	int ichannel0,ichannelf;
	ichannel0=0; ichannelf=nchannels;
	
	complex<double> eta0,eta1,eta2,eta3;
	complex<double> x1b,x2a,x2b,x3a,x3b,x,q1,q2,q3;
	complex<double> **M,*Y;
	complex<double> x1,x2;
	double F2,G2,F2prime,G2prime,qsquared,r;
	int i,j,iq,ir,ichannel;
	mu=m1*m2/(m1+m2);
	
	//for(ichannel=0;ichannel<nchannels;ichannel++){
	for(ichannel=ichannel0;ichannel<ichannelf;ichannel++){
		//printf("nwells[%d]=%d\n",ichannel,nwells[ichannel]);
		Avec=new complex<double>[2*nwells[ichannel]];
		if(nwells[ichannel]==1){
			for(iq=0;iq<nqmax;iq++){
				q=GetQ(iq);
				eta0=q1q2*mu*ALPHA/q;
				
				qsquared=q*q-2.0*mu*V0[ichannel][0];
				if(qsquared>0) q1=sqrt(qsquared);
				else q1=ci*sqrt(abs(qsquared));
				x1=q1*a[ichannel][0]/HBARC;
				eta1=q1q2*mu*ALPHA/q1;
				CoulWave::GetFGprime_ComplexQ(ell[ichannel],x1,eta1,&F1,&G1,&F1prime,&G1prime);
				x2=q*a[ichannel][0]/HBARC;
				CoulWave::GetFGprime_ComplexQ(ell[ichannel],x2,eta0,&F2,&G2,&F2prime,&G2prime);     
				beta=(abs(q1)/q)*F1prime/F1;
				delta[ichannel][iq]=-atan2(beta*F2-F2prime,beta*G2-G2prime);
				A[ichannel][iq][0]=0.5*(exp(-2.0*ci*delta[ichannel][iq])
					*(F2+ci*G2)+(F2-ci*G2))/F1;
				A[ichannel][iq][1]=exp(-2.0*ci*delta[ichannel][iq]);
				
			}			
		}
		else if(nwells[ichannel]==2){
			cmatrix=new CGSLMatrix_Complex(4);
			Y=new complex<double>[4];
			M=new complex<double> *[4];
			for(i=0;i<4;i++)
				M[i]=new complex<double>[4];
			
			for(iq=0;iq<nqmax;iq++){		
				q=GetQ(iq);
				q1=sqrt(abs(q*q-2.0*mu*V0[ichannel][0]));
				if(q*q-2.0*mu*V0[ichannel][0]<0.0)
					q1=ci*q1;
				q2=sqrt(abs(q*q-2.0*mu*V0[ichannel][1]));
				if(q*q-2.0*mu*V0[ichannel][1]<0.0)
					q2=ci*q2;
				x1b=a[ichannel][0]*q1/HBARC;
				x2a=a[ichannel][0]*q2/HBARC;
				x2b=a[ichannel][1]*q2/HBARC;
				x=a[ichannel][1]*q/HBARC;
				eta1=q1q2*mu*ALPHA/q1;
				eta2=q1q2*mu*ALPHA/q2;
				eta0=q1q2*mu*ALPHA/q;
				
				CoulWave::GetFGprime_ComplexQ(ell[ichannel],x1b,eta1,&F1b,&G1b,&F1bprime,&G1bprime);
				CoulWave::GetFGprime_ComplexQ(ell[ichannel],x2a,eta2,&F2a,&G2a,&F2aprime,&G2aprime);
				CoulWave::GetFGprime_ComplexQ(ell[ichannel],x2b,eta2,&F2b,&G2b,&F2bprime,&G2bprime);
				CoulWave::GetFGprime_ComplexQ(ell[ichannel],x,eta0,&F,&G,&Fprime,&Gprime);
				
				for(i=0;i<4;i++){
					Y[i]=0.0;
					for(j=0;j<4;j++) 
						M[i][j]=0.0;
				}
				M[0][0]=F1b;              M[0][1]=-F2a;              M[0][2]=-G2a;
				M[1][0]=abs(q1)*F1bprime; M[1][1]=-abs(q2)*F2aprime; M[1][2]=-abs(q2)*G2aprime;
				M[2][1]=F2b;              M[2][2]=G2b;               M[2][3]=-0.5*(F+ci*G);
				M[3][1]=abs(q2)*F2bprime; M[3][2]=abs(q2)*G2bprime;  M[3][3]=-0.5*q*(Fprime+ci*Gprime);
				
				Y[2]=0.5*(F-ci*G); Y[3]=0.5*q*(Fprime-ci*Gprime);
				//cmatrix->SolveLinearEqs(Y,M,A[ichannel][iq]);
				cmatrix->SolveLinearEqs(Y,M,Avec);
				for(int ia=0;ia<2*nwells[ichannel];ia++)
					A[ichannel][iq][ia]=Avec[ia];
				delta[ichannel][iq]=-0.5*atan2(imag(A[ichannel][iq][3]),real(A[ichannel][iq][3]));
				if(delta[ichannel][iq]<0.0)
					delta[ichannel][iq]=delta[ichannel][iq]+PI;
			
			}
			delete(cmatrix);
			delete [] Y;
			for(i=0;i<4;i++) delete [] M[i];
			delete [] M;
		}
		else if(nwells[ichannel]==3){
			cmatrix=new CGSLMatrix_Complex(6);
			Y=new complex<double>[6];
			M=new complex<double> *[6];
			for(i=0;i<6;i++)
				M[i]=new complex<double>[6];
			for(iq=0;iq<nqmax;iq++){
				q=GetQ(iq);
				q1=sqrt(fabs(q*q-2.0*mu*V0[ichannel][0]));
				if(q*q-2.0*mu*V0[ichannel][0]<0.0)
					q1=ci*q1;
				q2=sqrt(fabs(q*q-2.0*mu*V0[ichannel][1]));
				if(q*q-2.0*mu*V0[ichannel][1]<0.0)
					q2=ci*q2;
				q3=sqrt(fabs(q*q-2.0*mu*V0[ichannel][2]));
				if(q*q-2.0*mu*V0[ichannel][2]<0.0)
					q3=ci*q3;
				x1b=a[ichannel][0]*q1/HBARC;
				x2a=a[ichannel][0]*q2/HBARC;
				x2b=a[ichannel][1]*q2/HBARC;
				x3a=a[ichannel][1]*q3/HBARC;
				x3b=a[ichannel][2]*q3/HBARC;
				x=a[ichannel][2]*q/HBARC;
				eta1=q1q2*mu*ALPHA/q1;
				eta2=q1q2*mu*ALPHA/q2;
				eta3=q1q2*mu*ALPHA/q3;
				eta0=q1q2*mu*ALPHA/q;
				CoulWave::GetFGprime_ComplexQ(ell[ichannel],x1b,eta1,&F1b,&G1b,&F1bprime,&G1bprime);
				CoulWave::GetFGprime_ComplexQ(ell[ichannel],x2a,eta2,&F2a,&G2a,&F2aprime,&G2aprime);
				CoulWave::GetFGprime_ComplexQ(ell[ichannel],x2b,eta2,&F2b,&G2b,&F2bprime,&G2bprime);
				CoulWave::GetFGprime_ComplexQ(ell[ichannel],x3a,eta3,&F3a,&G3a,&F3aprime,&G3aprime);
				CoulWave::GetFGprime_ComplexQ(ell[ichannel],x3b,eta3,&F3b,&G3b,&F3bprime,&G3bprime);
				CoulWave::GetFGprime_ComplexQ(ell[ichannel],x,eta0,&F,&G,&Fprime,&Gprime);
				for(i=0;i<6;i++){
					Y[i]=0.0;
					for(j=0;j<6;j++) M[i][j]=0.0;
				}
				M[0][0]=F1b;              M[0][1]=-F2a;              M[0][2]=-G2a;
				M[1][0]=abs(q1)*F1bprime; M[1][1]=-abs(q2)*F2aprime; M[1][2]=-abs(q2)*G2aprime;
				M[2][1]=F2b;              M[2][2]=G2b;               M[2][3]=-F3a;               M[2][4]=-G3a;
				M[3][1]=abs(q2)*F2bprime; M[3][2]=abs(q2)*G2bprime;  M[3][3]=-abs(q3)*F3aprime;  M[3][4]=-abs(q3)*G3aprime;
				M[4][3]=F3b;              M[4][4]=G3b;               M[4][5]=-0.5*(F+ci*G);
				M[5][3]=abs(q3)*F3bprime; M[5][4]=abs(q3)*G3bprime;  M[5][5]=-0.5*q*(Fprime+ci*Gprime);
				
				Y[4]=0.5*(F-ci*G); Y[5]=0.5*q*(Fprime-ci*Gprime);
				//cmatrix->SolveLinearEqs(Y,M,A[ichannel][iq]);	
				cmatrix->SolveLinearEqs(Y,M,Avec);
				for(int ia=0;ia<2*nwells[ichannel];ia++)
					A[ichannel][iq][ia]=Avec[ia];
				delta[ichannel][iq]=-0.5*atan2(imag(A[ichannel][iq][5]),real(A[ichannel][iq][5]));
				if(delta[ichannel][iq]<0.0)
					delta[ichannel][iq]=delta[ichannel][iq]+PI;				
			}
			delete(cmatrix);
			delete [] Y;
			for(i=0;i<6;i++) delete [] M[i];
			delete [] M;
			
		}
		else{
			snprintf(message,CLog::CHARLENGTH,"nwells[%d] not equal to 1, 2 or 3??? =%d\n",ichannel,nwells[ichannel]);
			CLog::Fatal(message);
		}
		delete [] Avec;
	}
	
	vector<complex<double>> DPhiTestArray(nchannels);
	for(iq=0;iq<nqmax;iq++){
		for(ichannel=ichannel0;ichannel<ichannelf;ichannel++)
			DelPhiArray[iq][0][ichannel]=0.0;
		for(ir=1;ir<=DelPhiArray_NRMAX;ir++){
			for(ichannel=0;ichannel<nchannels;ichannel++){
				DPhiTestArray[ichannel]=0.0;
			}
			r=ir*DelPhiArray_DELR;
			SquareWell_CalcDelPhi(iq,r,DPhiTestArray);
			for(ichannel=0;ichannel<6;ichannel++){
				DelPhiArray[iq][ir][ichannel]=DPhiTestArray[ichannel];
			}
		}
	}
	
}

void CWaveFunction::SquareWell_GetDelPhi(int iq,double r,complex<double> *DelPhi){
	double wa,wb;
	vector<complex<double>> DPhiVec;
	DPhiVec.resize(nchannels);
	int ichannel;
	int ir=int(floor(r/DelPhiArray_DELR));
	if(ir<DelPhiArray_NRMAX){
		wb=(r-ir*DelPhiArray_DELR)/DelPhiArray_DELR;
		wa=1.0-wb;
		for(ichannel=0;ichannel<nchannels;ichannel++){
			DelPhi[ichannel]=wa*DelPhiArray[iq][ir][ichannel]+wb*DelPhiArray[iq][ir+1][ichannel];
		}
	}
	else{
		SquareWell_CalcDelPhi(iq,r,DPhiVec);
		for(ichannel=0;ichannel<nchannels;ichannel++)
			DelPhi[ichannel]=DPhiVec[ichannel];
	}
}


void CWaveFunction::SquareWell_CalcDelPhi(int iq,double r,vector<complex<double>> &DelPhi){
	complex<double> q1;
	double qsquared,eta0;
	complex<double> x1, eta1;
	double F,G,Fprime,Gprime;
	double F0[5],G0[5],F0prime[5],G0prime[5]; // assuming L is never bigger than 4
	int lexist[5]={0};
	int ichannel,iwell,l;
	double q=GetQ(iq);
	
	mu=m1*m2/(m1+m2);
	eta0=q1q2*mu*ALPHA/q;
	for(ichannel=0;ichannel<nchannels;ichannel++){
		if(lexist[ell[ichannel]]==0){
			l=ell[ichannel];
			CoulWave::GetFGprime_ComplexQ(l,0.0*ci+q*r/HBARC,0.0*ci+eta0,&F0[l],&G0[l],&F0prime[l],&G0prime[l]);
			lexist[l]=1;
		}
	}
	
	complex<double> phi,phi0,AA;
	complex<double> cgs;

	for(ichannel=0;ichannel<nchannels;ichannel++){
		l=ell[ichannel];
		cgs=cgsqwell[iq][l];
		if(r>a[ichannel][nwells[ichannel]-1]){
			AA=A[ichannel][iq][2*nwells[ichannel]-1];
			phi0=0.5*(F0[l]+ci*G0[l]);
			phi=AA*phi0;
			DelPhi[ichannel]=phi-phi0;
		}
		else{
			iwell=0;
			while(r>a[ichannel][iwell])
				iwell+=1;
			qsquared=q*q-2.0*mu*V0[ichannel][iwell];
			if(qsquared > 0.0)
				q1=sqrt(qsquared);
			else
				q1=ci*sqrt(fabs(qsquared));
			
			eta1=q1q2*mu*ALPHA/q1;
			x1=q1*r/HBARC;
			CoulWave::GetFGprime_ComplexQ(l,0.0*ci+x1,0.0*ci+eta1,&F,&G,&Fprime,&Gprime);
			phi0=F0[l];
			if(iwell==0){
				phi=A[ichannel][iq][0]*F;
			}
			else{
				phi=A[ichannel][iq][2*iwell-1]*F +A[ichannel][iq][2*iwell]*G;
			}
			DelPhi[ichannel]=phi-phi0;
		}
		DelPhi[ichannel]*=cgs;
	}
}

/*
void CWaveFunction::SquareWell_CalcDelPhi(int iq,double r,complex<double> *DelPhi){
	complex<double> q1;
	double qsquared,eta0;
	complex<double> x1, eta1;
	double F,G,Fprime,Gprime;
	double F0[5],G0[5],F0prime[5],G0prime[5]; // assuming L is never bigger than 4
	int lexist[5]={0};
	int ichannel,iwell,l;
	double q=GetQ(iq);
	
	printf("nchannels=%d, r=%g\n",nchannels,r);
	for(ichannel=0;ichannel<nchannels;ichannel++){
		printf("Delphi[%d]=(%g,%g)\n",ichannel,real(DelPhi[ichannel]),imag(DelPhi[ichannel]));
		for(int ia=0;ia<2*nwells[ichannel];ia++){
			printf("A[%d][%d][%d]=(%g,%g)\n",ichannel,iq,ia,real(A[ichannel][iq][ia]),imag(A[ichannel][iq][ia]));
		}
	}
	
	mu=m1*m2/(m1+m2);
	eta0=q1q2*mu*ALPHA/q;
	for(ichannel=0;ichannel<nchannels;ichannel++){
		if(lexist[ell[ichannel]]==0){
			l=ell[ichannel];
			CoulWave::GetFGprime_ComplexQ(l,0.0*ci+q*r/HBARC,0.0*ci+eta0,&F0[l],&G0[l],&F0prime[l],&G0prime[l]);
			lexist[l]=1;
		}
	}
	
	complex<double> phi,phi0,AA;
	complex<double>cgs;

	
	for(ichannel=0;ichannel<nchannels;ichannel++){
		l=ell[ichannel];
		printf("---- beginning loop for ichannel=%d, iq=%d, l=%d\n",ichannel,iq,l);
		cgs=cgsqwell[iq][l];
		if(r>a[ichannel][nwells[ichannel]-1]){
			printf("check cc\n");
			AA=A[ichannel][iq][2*nwells[ichannel]-1];
			phi0=0.5*(F0[l]+ci*G0[l]);
			phi=AA*phi0;
			DelPhi[ichannel]=phi-phi0;
			printf("check ccc\n");
		}
		else{
			printf("check dd\n");
			iwell=0;
			while(r>a[ichannel][iwell])
				iwell+=1;
			if(iwell>=nwells[ichannel]){
				printf("disaster, iwell too big\n");
				exit(1);
			}
			printf("check, iwell=%d\n",iwell);
			qsquared=q*q-2.0*mu*V0[ichannel][iwell];
			if(qsquared > 0.0)
				q1=sqrt(qsquared);
			else
				q1=ci*sqrt(fabs(qsquared));
			
			printf("check ee\n");
			eta1=q1q2*mu*ALPHA/q1;
			x1=q1*r/HBARC;
			printf("ell[ichannel]=%d\n",l);
			CoulWave::GetFGprime_ComplexQ(l,x1,eta1,&F,&G,&Fprime, &Gprime);
			printf("check ff, ichannel=%d, iq=%d\n",ichannel,iq);
			phi0=F0[l];
			if(iwell==0){
				phi=A[ichannel][iq][0]*F;
			}
			else{
				printf("F=%g, G=%g, ichannel=%d, iq=%d, iwell=%d\n",
				F,G,ichannel,iq,iwell);
				phi=A[ichannel][iq][2*iwell-1]*F +A[ichannel][iq][2*iwell]*G;
				printf("check gg\n");
			}
			printf("DDD\n");
			//printf("DDD, phi=(%g,%g), phi0=(%g,%g), ichannel=%d, iwell=%d, nwells=%d\n",
			//real(phi),imag(phi),real(phi0),imag(phi0),ichannel,iwell,nwells[ichannel]);
			DelPhi[ichannel]=phi-phi0;
			printf("check ddd ----\n");
		}
		DelPhi[ichannel]*=cgs;
		printf("check ddddddddddd ----\n");
	}
	printf("leaving happily\n");
}
*/

void CWaveFunction::SquareWell_CalcDelPhi2(int iq,double r,
complex<double> *DelPhi,complex<double> *DelPhiPrime,double *DelPhi2){
	complex<double> q1;
	double qsquared,eta0;
	complex<double> x1, eta1;
	double F,G,Fprime,Gprime;
	double F0[5],G0[5],F0prime[5],G0prime[5]; // assuming L is never bigger than 4
	int lexist[5]={0};
	int ichannel,iwell,l;
	double q=GetQ(iq);
	
	mu=m1*m2/(m1+m2);
	eta0=q1q2*mu*ALPHA/q;
	for(ichannel=0;ichannel<nchannels;ichannel++){
		if(lexist[ell[ichannel]]==0){
			l=ell[ichannel];
			CoulWave::GetFGprime_ComplexQ(l,0.0*ci+q*r/HBARC,0.0*ci+eta0,&F0[l],&G0[l],&F0prime[l],&G0prime[l]);
			lexist[l]=1;
		}
	}
	
	
	complex<double> phi,phi0,AA,phiprime,phi0prime;
	complex<double>cgs;
	
	
	iwell=0;

	
	for(ichannel=0;ichannel<nchannels;ichannel++){

		cgs=cgsqwell[iq][ell[ichannel]];
		if(r>a[ichannel][nwells[ichannel]-1]){
			AA=A[ichannel][iq][2*nwells[ichannel]-1];
			phi0=0.5*(F0[ell[ichannel]]+ci*G0[ell[ichannel]]);
			phi0prime=0.5*q*(F0prime[ell[ichannel]]+ci*G0prime[ell[ichannel]]);
			phi=AA*phi0;
			phiprime=AA*phi0prime;
			DelPhi[ichannel]=phi-phi0;
			DelPhiPrime[ichannel]=phiprime-phi0prime;
			
			phi=AA*phi0+conj(phi0);
			phi0=phi0+conj(phi0);
			//if(r>100 && ichannel==3)
			//printf("%9.3f (%g,%g)\n",r,real(phi0),imag(phi0));
			
			DelPhi2[ichannel]=real(phi*conj(phi))-real(phi0*conj(phi0));
		}
		else{
			iwell=0;
			while(r>a[ichannel][iwell])
				iwell+=1;
			qsquared=q*q-2.0*mu*V0[ichannel][iwell];
			if(qsquared > 0.0)
				q1=sqrt(qsquared);
			else
				q1=ci*sqrt(fabs(qsquared));
			
			eta1=q1q2*mu*ALPHA/q1;
			x1=q1*r/HBARC;
			CoulWave::GetFGprime_ComplexQ(ell[ichannel],x1,eta1,&F,&G,&Fprime, &Gprime);
			phi0=F0[ell[ichannel]];
			phi0prime=q*F0prime[ell[ichannel]];
			if(iwell==0){
				phi=A[ichannel][iq][0]*F;
				phiprime=A[ichannel][iq][0]*abs(q1)*Fprime;
			}
			else{
				phi=A[ichannel][iq][2*iwell-1]*F +A[ichannel][iq][2*iwell]*G;
				phiprime=abs(q1)*A[ichannel][iq][2*iwell-1]*Fprime+abs(q1)*A[ichannel][iq][2*iwell]*Gprime;
			}
			DelPhi[ichannel]=phi-phi0;
			DelPhiPrime[ichannel]=phiprime-phi0prime;
			DelPhi2[ichannel]=real(phi*conj(phi))-real(phi0*conj(phi0));
		}
		DelPhi[ichannel]*=cgs;
		DelPhiPrime[ichannel]*=cgs;
		DelPhi2[ichannel]*=real(cgs*conj(cgs));
	}
}

void CWaveFunction::SquareWell_MakeArrays(){
	int ichannel,iq,ir,l;
	double q,eta0;
	int *lexist;
	V0=new double *[nchannels];
	a=new double *[nchannels];
	A=new complex<double> **[nchannels];
	lexist=new int[ellmax+1];
	mu=m1*m2/(m1+m2);
	for(l=0;l<=ellmax;l++)
		lexist[l]=0;
	
	for(ichannel=0;ichannel<nchannels;ichannel++){
		if(lexist[ell[ichannel]]==0)
			lexist[ell[ichannel]]=1;
		V0[ichannel]=new double[nwells[ichannel]];
		a[ichannel]=new double[nwells[ichannel]];
		for(int ia=0;ia<2*nwells[ichannel];ia++)
			V0[ichannel][ia]=a[ichannel][ia]=0.0;
		A[ichannel]=new complex<double> *[nqmax];
		for(iq=0;iq<nqmax;iq++){
			A[ichannel][iq]=new complex<double>[2*nwells[ichannel]];
			for(int ia=0;ia<2*nwells[ichannel];ia++)
				A[ichannel][iq][ia]=0.0;
		}
	}
  
	cgsqwell=new complex<double> *[nqmax];
	for(iq=0;iq<nqmax;iq++){
		cgsqwell[iq]=new complex<double>[ellmax+1];
		q=GetQ(iq);
		eta0=q1q2*mu*ALPHA/q;
		for(l=0;l<=ellmax;l++){
			if(lexist[l]==1){
				cgsqwell[iq][l]=CoulWave::cgamma(l+1.0+ci*eta0);
				cgsqwell[iq][l]=conj(cgsqwell[iq][l]/std::abs(cgsqwell[iq][l]));
			}
		}
	}
	delete [] lexist;
	
	DelPhiArray_NRMAX=400;
	DelPhiArray_DELR=0.1;
	DelPhiArray=new complex<double> **[nqmax];
	for(iq=0;iq<nqmax;iq++){
		DelPhiArray[iq]=new complex<double> *[DelPhiArray_NRMAX+1];
		for(ir=0;ir<=DelPhiArray_NRMAX;ir++){
			DelPhiArray[iq][ir]=new complex<double>[nchannels];
			for(ichannel=0;ichannel<nchannels;ichannel++)
				DelPhiArray[iq][ir][ichannel]=0.0;
		}
	}
	
}

void CWaveFunction::SquareWell_DeleteArrays(){
	int ichannel,iq,ir;
	for(ichannel=0;ichannel<nchannels;ichannel++){
		for(iq=0;iq<nqmax;iq++)	delete [] A[ichannel][iq];
		delete [] A[ichannel];
		delete [] a[ichannel];
		delete [] V0[ichannel];
	}
	delete [] A;
	delete [] a;
	delete [] V0;
	delete [] nwells;
	
	for(iq=0;iq<nqmax;iq++){
		delete [] cgsqwell[iq];
	}
	delete [] cgsqwell;
	
	for(iq=0;iq<nqmax;iq++){
		for(ir=0;ir<=DelPhiArray_NRMAX;ir++) delete [] DelPhiArray[iq][ir];
		delete [] DelPhiArray[iq];
	}
	delete [] DelPhiArray;
	
}
