#include "msu_coral/source2cf.h"
#include "msu_coral/wavefunction.h"
#include "msu_coral/kernel.h"
#include "msu_commonutils/arrays.h"
#include "msu_commonutils/randy.h"
#include "msu_commonutils/log.h"

using namespace std;
using namespace NMSUPratt;

void S2CF::s2c(int lx,int ly,int lz,CCHArray *S,CKernel *kernel,CCHArray *CF){
	char message[CLog::CHARLENGTH];
	int lmax=kernel->GetLMAX();
	if(lmax>CF->GetLMAX() || lmax>S->GetLMAX()){
		snprintf(message,CLog::CHARLENGTH,"FATAL: Kernel LMAX=%d larger than either S LMAX=%d or CF LMAX=%d\n",
			lmax,S->GetLMAX(),CF->GetLMAX());
			CLog::Fatal(message);
	}
	if( (CF->GetLMAX()!=lmax) ){
		CLog::Info("WARNING: Array parameters for kernel calculations don't match\n");
		CLog::Info("___ CORRELATION FUNCTION PARAMETERS ___\n");
		CF->PrintPars();
		CLog::Info("_____ SOURCE FUNCTION PARAMETERS _____\n");
		S->PrintPars();
		snprintf(message,CLog::CHARLENGTH,"For kernel, LMAX=%d\n",kernel->GetLMAX());
		CLog::Info(message);
	}
	if(kernel->GetIDENTICAL() && (!CF->GetXSYM() || !CF->GetYSYM() || !CF->GetZSYM())){
		snprintf(message,CLog::CHARLENGTH,"FATAL: kernel has no odd L components, but CF wants them\n");
		CLog::Info(message);
		CLog::Fatal("Make sure CF array has XSYM=YSYM=ZSYM='true'\n");
	}
	int iq,ir,nqmax,nrmax,L;
	double r,delr,q,delq,norm;
	bool match=0;
	delr=S->GetRADSTEP();
	nrmax=S->GetNRADIAL();
	delq=CF->GetRADSTEP();
	nqmax=CF->GetNRADIAL();
	if(fabs(delr-kernel->GetDELR())<1.0E-5 
		&& fabs(delq-kernel->GetDELQ())<1.0E-5
		&& nrmax==kernel->GetNRMAX()
		&& nqmax==kernel->GetNQMAX()) match=1;
	CF->ZeroArray(lx,ly,lz);
	L=lx+ly+lz;
	for(iq=0;iq<nqmax;iq++){
		q=(0.5+iq)*delq;
		for(ir=0;ir<nrmax;ir++){
			r=(0.5+ir)*delr;
			norm=4.0*PI*r*r*delr;
			if(match) CF->IncrementElement(lx,ly,lz,iq,kernel->GetValue(L,iq,ir)*norm*S->GetElement(lx,ly,lz,ir));
			else CF->IncrementElement(lx,ly,lz,iq,kernel->GetValue(L,q,r)*norm*S->GetElement(lx,ly,lz,r));
		}
	}
}

void S2CF::s2c(CCHArray *S,CKernel *kernel,CCHArray *CF){
	char message[CLog::CHARLENGTH];
	int lmax=kernel->GetLMAX();
	if(lmax>CF->GetLMAX() || lmax>S->GetLMAX()){
		snprintf(message,CLog::CHARLENGTH,"Kernel LMAX=%d larger than either S LMAX=%d or CF LMAX=%d\n",
			lmax,S->GetLMAX(),CF->GetLMAX());
		CLog::Fatal(message);
	}
	if( (CF->GetLMAX()!=lmax) ){
		snprintf(message,CLog::CHARLENGTH,"WARNING: Array parameters for kernel calculations don't match\n");
		CLog::Info(message);
		CLog::Info("___ CORRELATION FUNCTION PARAMETERS ___\n");
		CF->PrintPars();
		CLog::Info("_____ SOURCE FUNCTION PARAMETERS _____\n");
		S->PrintPars();
		snprintf(message,CLog::CHARLENGTH,"For kernel, LMAX=%d\n",kernel->GetLMAX());
		CLog::Info(message);
	}
	if(kernel->GetIDENTICAL() && (!CF->GetXSYM() || !CF->GetYSYM() || !CF->GetZSYM())){
		CLog::Info("FATAL: kernel has no odd L components, but CF wants them\n");
		CLog::Fatal("Make sure CF array has XSYM=YSYM=ZSYM='true'\n");;
	}
	int iq,ir,nqmax,nrmax,L,lx,ly,lz,dlx,dly,dlz;
	double r,delr,q,delq,norm;
	bool match=0;
	delr=S->GetRADSTEP();
	nrmax=S->GetNRADIAL();
	delq=CF->GetRADSTEP();
	nqmax=CF->GetNRADIAL();
	dlx=dly=dlz=1;
	if(CF->GetXSYM()) dlx=2;
	if(CF->GetYSYM()) dly=2;
	if(CF->GetZSYM()) dlz=2;
	if(fabs(delr-kernel->GetDELR())<1.0E-5 
		&& fabs(delq-kernel->GetDELQ())<1.0E-5
		&& nrmax==kernel->GetNRMAX()
		&& nqmax==kernel->GetNQMAX()) match=1;
	CF->ZeroArray();
	for(iq=0;iq<nqmax;iq++){
		q=(0.5+iq)*delq;
		for(lx=0;lx<2;lx+=dlx){
			for(ly=0;ly<=lmax-lx;ly+=dly){
				for(lz=0;lz<=lmax-lx-ly;lz+=dlz){
					L=lx+ly+lz;
					for(ir=0;ir<nrmax;ir++){
						r=(0.5+ir)*delr;
						norm=4.0*PI*r*r*delr;
						if(match) CF->IncrementElement(lx,ly,lz,iq,kernel->GetValue(L,iq,ir)
							*norm*S->GetElement(lx,ly,lz,ir));
						else CF->IncrementElement(lx,ly,lz,iq,kernel->GetValue(L,q,r)
							*norm*S->GetElement(lx,ly,lz,r));
					}
				}
			}
		}
	}
}

void S2CF::s2c(C3DArray *S,CKernelWF *kernel,C3DArray *CF){
	char message[CLog::CHARLENGTH];
	int ix,iy,iz,isx,isy,isz,jx,jy,jz,jsx,jsy,jsz;
	int nsx,nsy,nsz;
	long long int icalc,ncalc;
	int nxmax=S->GetNXMAX();
	int nymax=S->GetNYMAX();
	int nzmax=S->GetNZMAX();
	double delx=S->GetDELX();
	double dely=S->GetDELY();
	double delz=S->GetDELZ();
	int nqxmax=CF->GetNXMAX();
	int nqymax=CF->GetNYMAX();
	int nqzmax=CF->GetNZMAX();
	double delqx=CF->GetDELX();
	double delqy=CF->GetDELY();
	double delqz=CF->GetDELZ();
	double x,y,z,qx,qy,qz,norm,q,r,ctheta,wf2,svalue;
	bool IDENTICAL=kernel->GetIDENTICAL();
	bool XSYM=S->GetXSYM();
	bool YSYM=S->GetYSYM();
	bool ZSYM=S->GetZSYM();
	if(IDENTICAL&&!(XSYM&&YSYM&&ZSYM)){
		CLog::Info("S2CF::s2c, kernel says particles are identical, but symmetries are not all even\n");
		snprintf(message,CLog::CHARLENGTH,"XSYM=%d, YSYM=%d, ZSYM=%d\n",int(XSYM),int(YSYM),int(ZSYM));
		CLog::Fatal(message);
	}

	norm=delx*dely*delz;
	nsx=nsy=nsz=2;
	if(XSYM){
		nsx=1;
		norm*=2.0;
	}
	if(YSYM){
		nsy=1;
		norm*=2.0;
	}
	if(ZSYM){
		nsz=1;
		norm*=2.0;
	}
	CF->ZeroArray();
	ncalc=nxmax*nymax*nzmax*nsx*nsy*nsz;
	ncalc=ncalc/10;
	icalc=0;

	for(isx=0;isx<nsx;isx++){
		for(ix=0;ix<nxmax;ix++){
			x=(0.5+ix)*delx;
			if(isx==1) x=-x;

			for(isy=0;isy<nsy;isy++){
				for(iy=0;iy<nymax;iy++){
					y=(0.5+iy)*dely;
					if(isy==1) y=-y;

					for(isz=0;isz<nsz;isz++){
						for(iz=0;iz<nzmax;iz++){
							z=(0.5+iz)*delz;
							if(isz==1) z=-z;

							r=sqrt(x*x+y*y+z*z);
							svalue=S->GetElement(isx,ix,isy,iy,isz,iz);
							//svalue=exp(-x*x/(4.0*9.0)-y*y/(4.0*25)-z*z/(4.0*49));
							//svalue=svalue/(pow(4.0*PI*9.0*25.0*49.0);

							for(jsx=0;jsx<nsx;jsx++){
								for(jx=0;jx<nqxmax;jx++){
									qx=(0.5+jx)*delqx;
									if(jsx==1) qx=-qx;

									for(jsy=0;jsy<nsy;jsy++){
										for(jy=0;jy<nqymax;jy++){
											qy=(0.5+jy)*delqy;
											if(jsy==1) qy=-qy;

											for(jsz=0;jsz<nsz;jsz++){
												for(jz=0;jz<nqzmax;jz++){
													qz=(0.5+jz)*delqz;
													if(jsz==1) qz=-qz;
													q=sqrt(qx*qx+qy*qy+qz*qz);
													wf2=0.0;

													if(XSYM&YSYM&&ZSYM){
														ctheta=(qx*x+qy*y+qz*z)/(q*r);
														wf2+=0.25*kernel->GetPsiSquared(q,r,ctheta);
														if(!IDENTICAL) wf2+=0.25*kernel->GetPsiSquared(q,r,-ctheta);
														ctheta=(-qx*x+qy*y+qz*z)/(q*r);
														wf2+=0.25*kernel->GetPsiSquared(q,r,ctheta);
														if(!IDENTICAL) wf2+=0.25*kernel->GetPsiSquared(q,r,-ctheta);
														ctheta=(qx*x-qy*y+qz*z)/(q*r);
														wf2+=0.25*kernel->GetPsiSquared(q,r,ctheta);
														if(!IDENTICAL) wf2+=0.25*kernel->GetPsiSquared(q,r,-ctheta);																																																				
														ctheta=(qx*x+qy*y-qz*z)/(q*r);
														wf2+=0.25*kernel->GetPsiSquared(q,r,ctheta);
														if(!IDENTICAL) wf2+=0.25*kernel->GetPsiSquared(q,r,-ctheta);
														if(!IDENTICAL) wf2=0.5*wf2;
													}
													else if(!XSYM && YSYM && ZSYM){
														ctheta=(qx*x+qy*y+qz*z)/(q*r);
														wf2+=0.25*kernel->GetPsiSquared(q,r,ctheta);
														ctheta=(qx*x-qy*y+qz*z)/(q*r);
														wf2+=0.25*kernel->GetPsiSquared(q,r,ctheta);
														ctheta=(qx*x+qy*y-qz*z)/(q*r);
														wf2+=0.25*kernel->GetPsiSquared(q,r,ctheta);
														ctheta=(qx*x-qy*y-qz*z)/(q*r);
														wf2+=0.25*kernel->GetPsiSquared(q,r,ctheta);
													}
													else if(XSYM && !YSYM && ZSYM){
														ctheta=(qx*x+qy*y+qz*z)/(q*r);
														wf2+=0.25*kernel->GetPsiSquared(q,r,ctheta);
														ctheta=(-qx*x+qy*y+qz*z)/(q*r);
														wf2+=0.25*kernel->GetPsiSquared(q,r,ctheta);
														ctheta=(qx*x+qy*y-qz*z)/(q*r);
														wf2+=0.25*kernel->GetPsiSquared(q,r,ctheta);
														ctheta=(-qx*x+qy*y-qz*z)/(q*r);
														wf2+=0.25*kernel->GetPsiSquared(q,r,ctheta);
													}
													else if(XSYM && YSYM && !ZSYM){
														ctheta=(qx*x+qy*y+qz*z)/(q*r);
														wf2+=0.25*kernel->GetPsiSquared(q,r,ctheta);
														ctheta=(-qx*x+qy*y+qz*z)/(q*r);
														wf2+=0.25*kernel->GetPsiSquared(q,r,ctheta);
														ctheta=(qx*x-qy*y+qz*z)/(q*r);
														wf2+=0.25*kernel->GetPsiSquared(q,r,ctheta);
														ctheta=(-qx*x-qy*y+qz*z)/(q*r);
														wf2+=0.25*kernel->GetPsiSquared(q,r,ctheta);
													}
													else if(XSYM && !YSYM && !ZSYM){
														ctheta=(qx*x+qy*y+qz*z)/(q*r);
														wf2+=0.5*kernel->GetPsiSquared(q,r,ctheta);
														ctheta=(-qx*x+qy*y+qz*z)/(q*r);
														wf2+=0.5*kernel->GetPsiSquared(q,r,ctheta);
													}
													else if(!XSYM && YSYM && !ZSYM){
														ctheta=(qx*x+qy*y+qz*z)/(q*r);
														wf2+=0.5*kernel->GetPsiSquared(q,r,ctheta);
														ctheta=(qx*x-qy*y+qz*z)/(q*r);
														wf2+=0.5*kernel->GetPsiSquared(q,r,ctheta);
													}
													else if(!XSYM && !YSYM && ZSYM){
														ctheta=(qx*x+qy*y+qz*z)/(q*r);
														wf2+=0.5*kernel->GetPsiSquared(q,r,ctheta);
														ctheta=(qx*x+qy*y-qz*z)/(q*r);
														wf2+=0.5*kernel->GetPsiSquared(q,r,ctheta);
													}
													else{
														ctheta=(qx*x+qy*y+qz*z)/(q*r);
														wf2+=kernel->GetPsiSquared(q,r,ctheta);
													}
													CF->IncrementElement(jsx,jx,jsy,jy,jsz,jz,
														norm*(wf2-1.0)*svalue);
												}
											}
										}
									}
								}
							}
							icalc+=1;
							if(icalc%ncalc==0){
								snprintf(message,CLog::CHARLENGTH,"s2c, finished %g percent\n",10*double(icalc)/double(ncalc));
								CLog::Info(message);
							}
						}
					}
				}
			}
		}
	}
}

void S2CF::s2c(C3DArray *S,CWaveFunction *wf,C3DArray *CF){
	char message[CLog::CHARLENGTH];
	int ix,iy,iz,isx,isy,isz,jx,jy,jz,jsx,jsy,jsz;
	int nsx,nsy,nsz;
	int nxmax=S->GetNXMAX();
	int nymax=S->GetNYMAX();
	int nzmax=S->GetNZMAX();
	double delx=S->GetDELX();
	double dely=S->GetDELY();
	double delz=S->GetDELZ();
	int nqxmax=CF->GetNXMAX();
	int nqymax=CF->GetNYMAX();
	int nqzmax=CF->GetNZMAX();
	double delqx=CF->GetDELX();
	double delqy=CF->GetDELY();
	double delqz=CF->GetDELZ();
	double x,y,z,qx,qy,qz,norm,q,r,ctheta,wf2,svalue;
	bool IDENTICAL=wf->GetIDENTICAL();
	bool XSYM=S->GetXSYM();
	bool YSYM=S->GetYSYM();
	bool ZSYM=S->GetZSYM();
	if(IDENTICAL&&!(XSYM&&YSYM&&ZSYM)){
		CLog::Info("S2CF::s2c, wf says particles are identical, but symmetries are not all even\n");
		snprintf(message,CLog::CHARLENGTH,"XSYM=%d, YSYM=%d, ZSYM=%d\n",int(XSYM),int(YSYM),int(ZSYM));
		CLog::Fatal(message);
	}

	norm=delx*dely*delz;
	nsx=nsy=nsz=2;
	if(XSYM){
		nsx=1;
		norm*=2.0;
	}
	if(YSYM){
		nsy=1;
		norm*=2.0;
	}
	if(ZSYM){
		nsz=1;
		norm*=2.0;
	}
	CF->ZeroArray();

	for(isx=0;isx<nsx;isx++){
		for(ix=0;ix<nxmax;ix++){
			x=(0.5+ix)*delx;
			if(isx==1) x=-x;

			for(isy=0;isy<nsy;isy++){
				for(iy=0;iy<nymax;iy++){
					y=(0.5+iy)*dely;
					if(isy==1) y=-y;

					for(isz=0;isz<nsz;isz++){
						for(iz=0;iz<nzmax;iz++){
							z=(0.5+iz)*delz;
							if(isz==1) z=-z;

							r=sqrt(x*x+y*y+z*z);
							//svalue=S->GetElement(isx,ix,isy,iy,isz,iz);
							svalue=exp(-0.25*((x*x/9.0)+(y*y/25.0)+(z*z/49.0)));
							svalue=svalue/(3.0*5.0*7.0*pow(4.0*PI,1.5));

							for(jsx=0;jsx<nsx;jsx++){
								for(jx=0;jx<nqxmax;jx++){
									qx=(0.5+jx)*delqx;
									if(jsx==1) qx=-qx;

									for(jsy=0;jsy<nsy;jsy++){
										for(jy=0;jy<nqymax;jy++){
											qy=(0.5+jy)*delqy;
											if(jsy==1) qy=-qy;

											for(jsz=0;jsz<nsz;jsz++){
												for(jz=0;jz<nqzmax;jz++){
													qz=(0.5+jz)*delqz;
													if(jsz==1) qz=-qz;
													q=sqrt(qx*qx+qy*qy+qz*qz);
													wf2=0.0;

													if(XSYM&YSYM&&ZSYM){
														ctheta=(qx*x+qy*y+qz*z)/(q*r);
														wf2+=0.25*wf->GetPsiSquared(q,r,ctheta);
														if(!IDENTICAL) wf2+=0.25*wf->GetPsiSquared(q,r,-ctheta);
														ctheta=(-qx*x+qy*y+qz*z)/(q*r);
														wf2+=0.25*wf->GetPsiSquared(q,r,ctheta);
														if(!IDENTICAL) wf2+=0.25*wf->GetPsiSquared(q,r,-ctheta);
														ctheta=(qx*x-qy*y+qz*z)/(q*r);
														wf2+=0.25*wf->GetPsiSquared(q,r,ctheta);
														if(!IDENTICAL) wf2+=0.25*wf->GetPsiSquared(q,r,-ctheta);																																																				
														ctheta=(qx*x+qy*y-qz*z)/(q*r);
														wf2+=0.25*wf->GetPsiSquared(q,r,ctheta);
														if(!IDENTICAL) wf2+=0.25*wf->GetPsiSquared(q,r,-ctheta);
														if(!IDENTICAL) wf2=0.5*wf2;
													}
													else if(!XSYM && YSYM && ZSYM){
														ctheta=(qx*x+qy*y+qz*z)/(q*r);
														wf2+=0.25*wf->GetPsiSquared(q,r,ctheta);
														ctheta=(qx*x-qy*y+qz*z)/(q*r);
														wf2+=0.25*wf->GetPsiSquared(q,r,ctheta);
														ctheta=(qx*x+qy*y-qz*z)/(q*r);
														wf2+=0.25*wf->GetPsiSquared(q,r,ctheta);
														ctheta=(qx*x-qy*y-qz*z)/(q*r);
														wf2+=0.25*wf->GetPsiSquared(q,r,ctheta);
													}
													else if(XSYM && !YSYM && ZSYM){
														ctheta=(qx*x+qy*y+qz*z)/(q*r);
														wf2+=0.25*wf->GetPsiSquared(q,r,ctheta);
														ctheta=(-qx*x+qy*y+qz*z)/(q*r);
														wf2+=0.25*wf->GetPsiSquared(q,r,ctheta);
														ctheta=(qx*x+qy*y-qz*z)/(q*r);
														wf2+=0.25*wf->GetPsiSquared(q,r,ctheta);
														ctheta=(-qx*x+qy*y-qz*z)/(q*r);
														wf2+=0.25*wf->GetPsiSquared(q,r,ctheta);
													}
													else if(XSYM && YSYM && !ZSYM){
														ctheta=(qx*x+qy*y+qz*z)/(q*r);
														wf2+=0.25*wf->GetPsiSquared(q,r,ctheta);
														ctheta=(-qx*x+qy*y+qz*z)/(q*r);
														wf2+=0.25*wf->GetPsiSquared(q,r,ctheta);
														ctheta=(qx*x-qy*y+qz*z)/(q*r);
														wf2+=0.25*wf->GetPsiSquared(q,r,ctheta);
														ctheta=(-qx*x-qy*y+qz*z)/(q*r);
														wf2+=0.25*wf->GetPsiSquared(q,r,ctheta);
													}
													else if(XSYM && !YSYM && !ZSYM){
														ctheta=(qx*x+qy*y+qz*z)/(q*r);
														wf2+=0.5*wf->GetPsiSquared(q,r,ctheta);
														ctheta=(-qx*x+qy*y+qz*z)/(q*r);
														wf2+=0.5*wf->GetPsiSquared(q,r,ctheta);
													}
													else if(!XSYM && YSYM && !ZSYM){
														ctheta=(qx*x+qy*y+qz*z)/(q*r);
														wf2+=0.5*wf->GetPsiSquared(q,r,ctheta);
														ctheta=(qx*x-qy*y+qz*z)/(q*r);
														wf2+=0.5*wf->GetPsiSquared(q,r,ctheta);
													}
													else if(!XSYM && !YSYM && ZSYM){
														ctheta=(qx*x+qy*y+qz*z)/(q*r);
														wf2+=0.5*wf->GetPsiSquared(q,r,ctheta);
														ctheta=(qx*x+qy*y-qz*z)/(q*r);
														wf2+=0.5*wf->GetPsiSquared(q,r,ctheta);
													}
													else{
														ctheta=(qx*x+qy*y+qz*z)/(q*r);
														wf2+=wf->GetPsiSquared(q,r,ctheta);
													}
													CF->IncrementElement(jsx,jx,jsy,jy,jsz,jz,
														norm*(wf2-1.0)*svalue);
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
}

void S2CF::s2c(CMCList *lista,CMCList *listb,CKernelWF *kernel,C3DArray *cf){
	char message[CLog::CHARLENGTH];
	bool samelists;
	int ia,ib,na,nb,jx,jy,jz,jsx,jsy,jsz;
	int nsx,nsy,nsz;
	int nqxmax=cf->GetNXMAX();
	int nqymax=cf->GetNYMAX();
	int nqzmax=cf->GetNZMAX();
	int npairs,ipair,jpair;
	double delqx=cf->GetDELX();
	double delqy=cf->GetDELY();
	double delqz=cf->GetDELZ();
	double x,y,z,qx,qy,qz,norm,q,r,ctheta,wf2,*ra,*rb;
	double rmax=kernel->GetNRMAX()*kernel->GetDELR();
	if(lista==listb){
		npairs=lista->GetNMC()*(lista->GetNMC()-1)/2;
		samelists=true;
	}
	else{
		npairs=lista->GetNMC()*listb->GetNMC();
		samelists=false;
	}
	norm=lista->GetNorm()*listb->GetNorm()/double(npairs);
	if(samelists && !kernel->GetIDENTICAL()) norm=0.5*norm;

	nsx=nsy=nsz=2;
	if(cf->GetXSYM()){
		nsx=1;
	}
	if(cf->GetYSYM()){
		nsy=1;
	}
	if(cf->GetZSYM()){
		nsz=1;
	}
	cf->ZeroArray();

	na=lista->GetNMC();
	nb=lista->GetNMC();
	ipair=jpair=0;
	for(ia=0;ia<na;ia++){
		ra=lista->GetR(ia);
		if(samelists) nb=ia-1;
		for(ib=0;ib<nb;ib++){
			rb=listb->GetR(ib);
			x=ra[1]-rb[1];
			y=ra[2]-rb[2];
			z=ra[3]-rb[3];
			r=sqrt(x*x+y*y+z*z);
			if(r<rmax){
				for(jsx=0;jsx<nsx;jsx++){
					for(jx=0;jx<nqxmax;jx++){
						qx=(0.5+jx)*delqx;
						if(jsx==1) qx=-qx;
						for(jsy=0;jsy<nsy;jsy++){
							for(jy=0;jy<nqymax;jy++){
								qy=(0.5+jy)*delqy;
								if(jsy==1) qy=-qy;
								for(jsz=0;jsz<nsz;jsz++){
									for(jz=0;jz<nqzmax;jz++){
										qz=(0.5+jz)*delqz;
										if(jsz==1) qz=-qz;
										q=sqrt(qx*qx+qy*qy+qz*qz);
										ctheta=(qx*x+qy*y+qz*z)/(q*r);
										wf2=kernel->GetPsiSquared(q,r,ctheta);
										cf->IncrementElement(jsx,jx,jsy,jy,jsz,jz,norm*(wf2-1.0));
										if(samelists && !kernel->GetIDENTICAL()){
											wf2=kernel->GetPsiSquared(q,r,-ctheta);
											cf->IncrementElement(jsx,jx,jsy,jy,jsz,jz,norm*(wf2-1.0));
										}

									}
								}
							}
						}
					}
				}
				ipair+=1;
				if(ipair>=npairs/10){
					jpair+=1;
					snprintf(message,CLog::CHARLENGTH,"Finished %d percent\n",jpair*10);
					CLog::Info(message);
					ipair=0;
				}
			}
		}
	}
}

void S2CF::s2c(CMCList *lista,CMCList *listb,CKernel *kernel,C3DArray *cf){
	char message[CLog::CHARLENGTH];
	bool samelists;
	int ia,ib,na,nb,npairs,jx,jy,jz,jsx,jsy,jsz,ipair,jpair;
	int nsx,nsy,nsz;
	int nqxmax=cf->GetNXMAX();
	int nqymax=cf->GetNYMAX();
	int nqzmax=cf->GetNZMAX();
	double delqx=cf->GetDELX();
	double delqy=cf->GetDELY();
	double delqz=cf->GetDELZ();
	double x,y,z,qx,qy,qz,norm,q,r,ctheta,wf2,*ra,*rb;
	double rmax=kernel->GetDELR()*kernel->GetNRMAX();
	if(lista==listb){
		npairs=lista->GetNMC()*(lista->GetNMC()-1)/2;
		samelists=true;
	}
	else{
		npairs=lista->GetNMC()*listb->GetNMC();
		samelists=false;
	}
	norm=lista->GetNorm()*listb->GetNorm()/double(npairs);
	if(samelists && !kernel->GetIDENTICAL()) norm=0.5*norm;

	nsx=nsy=nsz=2;
	if(cf->GetXSYM()){
		nsx=1;
	}
	if(cf->GetYSYM()){
		nsy=1;
	}
	if(cf->GetZSYM()){
		nsz=1;
	}
	cf->ZeroArray();

	na=lista->GetNMC();
	nb=listb->GetNMC();
	ipair=jpair=0;
	for(ia=0;ia<na;ia++){
		ra=lista->GetR(ia);
		if(samelists) nb=ia-1;
		for(ib=0;ib<nb;ib++){
			rb=listb->GetR(ib);
			x=ra[1]-rb[1];
			y=ra[2]-rb[2];
			z=ra[3]-rb[3];
			r=sqrt(x*x+y*y+z*z);
			if(r<rmax){
				for(jsx=0;jsx<nsx;jsx++){
					for(jx=0;jx<nqxmax;jx++){
						qx=(0.5+jx)*delqx;
						if(jsx==1) qx=-qx;
						for(jsy=0;jsy<nsy;jsy++){
							for(jy=0;jy<nqymax;jy++){
								qy=(0.5+jy)*delqy;
								if(jsy==1) qy=-qy;
								for(jsz=0;jsz<nsz;jsz++){
									for(jz=0;jz<nqzmax;jz++){
										qz=(0.5+jz)*delqz;
										if(jsz==1) qz=-qz;
										q=sqrt(qx*qx+qy*qy+qz*qz);
										ctheta=(qx*x+qy*y+qz*z)/(q*r);
										wf2=kernel->GetPsiSquared(q,r,ctheta);
										cf->IncrementElement(jsx,jx,jsy,jy,jsz,jz,norm*(wf2-1.0));
										if(samelists && !kernel->GetIDENTICAL()){
											wf2=kernel->GetPsiSquared(q,r,-ctheta);
											cf->IncrementElement(jsx,jx,jsy,jy,jsz,jz,norm*(wf2-1.0));
										}

									}
								}
							}
						}
					}
				}
			}
			ipair+=1;
			if(ipair==npairs/10){
				jpair+=1;
				snprintf(message,CLog::CHARLENGTH,"Finished %d percent\n",jpair*10);
				CLog::Info(message);
				ipair=0;
			}
		}
	}
}

void S2CF::s2c(CMCList *lista,CMCList *listb,CKernelWF *kernel,C3DArray *cf,int NMC){
	char message[CLog::CHARLENGTH];
	Crandy randy(-time(NULL));
	bool samelists;
	int ia,ib,na,jx,jy,jz,jsx,jsy,jsz,imc=0,jmc=0;
	int nsx,nsy,nsz;
	int nqxmax=cf->GetNXMAX();
	int nqymax=cf->GetNYMAX();
	int nqzmax=cf->GetNZMAX();
	double delqx=cf->GetDELX();
	double delqy=cf->GetDELY();
	double delqz=cf->GetDELZ();
	double x,y,z,qx,qy,qz,norm,q,r,ctheta,wf2,*ra,*rb;
	double rmax=kernel->GetDELR()*kernel->GetNRMAX();
	if(lista==listb){
		samelists=true;
	}
	else{
		samelists=false;
	}
	norm=lista->GetNorm()*listb->GetNorm()/double(NMC);

	nsx=nsy=nsz=2;
	if(cf->GetXSYM()){
		nsx=1;
	}
	if(cf->GetYSYM()){
		nsy=1;
	}
	if(cf->GetZSYM()){
		nsz=1;
	}
	cf->ZeroArray();
	na=lista->GetNMC();

	for(jsx=0;jsx<nsx;jsx++){
		for(jx=0;jx<nqxmax;jx++){
			qx=(0.5+jx)*delqx;
			if(jsx==1) qx=-qx;
			for(jsy=0;jsy<nsy;jsy++){
				for(jy=0;jy<nqymax;jy++){
					qy=(0.5+jy)*delqy;
					if(jsy==1) qy=-qy;
					for(jsz=0;jsz<nsz;jsz++){
						for(jz=0;jz<nqzmax;jz++){
							qz=(0.5+jz)*delqz;
							if(jsz==1) qz=-qz;
							q=sqrt(qx*qx+qy*qy+qz*qz);
							imc=jmc=0;
							for(imc=0;imc<NMC;imc++){
								jmc+=1;
								ia=rint(floor(randy.ran()*na));
								do{
									ib=rint(floor(randy.ran()*na));
								} while(ia==ib && samelists);
								ra=lista->GetR(ia);
								rb=listb->GetR(ib);
								x=ra[1]-rb[1];
								y=ra[2]-rb[2];
								z=ra[3]-rb[3];
								r=sqrt(x*x+y*y+z*z);
								if(r<rmax){
									ctheta=(qx*x+qy*y+qz*z)/(q*r);
									wf2=kernel->GetPsiSquared(q,r,ctheta);
									cf->IncrementElement(jsx,jx,jsy,jy,jsz,jz,norm*(wf2-1.0));
									if(samelists && !kernel->GetIDENTICAL()){
										wf2=kernel->GetPsiSquared(q,r,-ctheta);
										cf->IncrementElement(jsx,jx,jsy,jy,jsz,jz,norm*(wf2-1.0));
									}

								}
							}
						}
					}
				}
			}
		}
		if(jmc==NMC/10){
			snprintf(message,CLog::CHARLENGTH,"finished %g percent\n",double(imc+1)*100/double(NMC));
			CLog::Info(message);
			jmc=0;
		}
	}
}

void S2CF::s2c(CMCList *lista,CMCList *listb,CWaveFunction *wf,C3DArray *cf){
	char message[CLog::CHARLENGTH];
	bool samelists;
	int ia,ib,na,nb,npairs,ipair,jpair,jx,jy,jz,jsx,jsy,jsz;
	int nsx,nsy,nsz;
	int nqxmax=cf->GetNXMAX();
	int nqymax=cf->GetNYMAX();
	int nqzmax=cf->GetNZMAX();
	double delqx=cf->GetDELX();
	double delqy=cf->GetDELY();
	double delqz=cf->GetDELZ();
	double x,y,z,qx,qy,qz,norm,q,r,ctheta,wf2,*ra,*rb;
	if(lista==listb){
		npairs=lista->GetNMC()*(lista->GetNMC()-1)/2;
		samelists=true;
	}
	else{
		npairs=lista->GetNMC()*listb->GetNMC();
		samelists=false;
	}
	norm=lista->GetNorm()*listb->GetNorm()/double(npairs);
	if(samelists && !wf->GetIDENTICAL()) norm=0.5*norm;

	nsx=nsy=nsz=2;
	if(cf->GetXSYM()){
		nsx=1;
	}
	if(cf->GetYSYM()){
		nsy=1;
	}
	if(cf->GetZSYM()){
		nsz=1;
	}
	cf->ZeroArray();

	na=lista->GetNMC();
	nb=listb->GetNMC();
	ipair=jpair=0;
	for(ia=0;ia<na;ia++){
		ra=lista->GetR(ia);
		if(samelists) nb=ia-1;
		for(ib=0;ib<nb;ib++){
			rb=listb->GetR(ib);
			x=ra[1]-rb[1];
			y=ra[2]-rb[2];
			z=ra[3]-rb[3];
			r=sqrt(x*x+y*y+z*z);
			for(jsx=0;jsx<nsx;jsx++){
				for(jx=0;jx<nqxmax;jx++){
					qx=(0.5+jx)*delqx;
					if(jsx==1) qx=-qx;
					for(jsy=0;jsy<nsy;jsy++){
						for(jy=0;jy<nqymax;jy++){
							qy=(0.5+jy)*delqy;
							if(jsy==1) qy=-qy;
							for(jsz=0;jsz<nsz;jsz++){
								for(jz=0;jz<nqzmax;jz++){
									qz=(0.5+jz)*delqz;
									if(jsz==1) qz=-qz;
									q=sqrt(qx*qx+qy*qy+qz*qz);
									ctheta=(qx*x+qy*y+qz*z)/(q*r);
									wf2=wf->GetPsiSquared(q,r,ctheta);
									cf->IncrementElement(jsx,jx,jsy,jy,jsz,jz,norm*(wf2-1.0));
									if(samelists && !wf->GetIDENTICAL()){
										wf2=wf->GetPsiSquared(q,r,-ctheta);
										cf->IncrementElement(jsx,jx,jsy,jy,jsz,jz,norm*(wf2-1.0));
									}

								}
							}
						}
					}
				}
			}
			ipair+=1;
			if(ipair>=npairs/10){
				jpair+=1;
				snprintf(message,CLog::CHARLENGTH,"Finished %d percent\n",jpair*10);
				CLog::Info(message);
				ipair=0;
			}
		}
	}
}

void S2CF::s2c_bosons(CMCList *list,C3DArray *cf){
	int i,n,jx,jy,jz;
	int nqxmax=cf->GetNXMAX();
	int nqymax=cf->GetNYMAX();
	int nqzmax=cf->GetNZMAX();
	double delqx=cf->GetDELX();
	double delqy=cf->GetDELY();
	double delqz=cf->GetDELZ();
	double c,r2,norm,arg,qx,qy,qz,*r;
	complex<double> alpha;
	n=list->GetNMC();
	norm=list->GetNorm()*list->GetNorm();

	cf->ZeroArray();
	for(jx=0;jx<nqxmax;jx++){
		qx=(0.5+jx)*delqx;
		for(jy=0;jy<nqymax;jy++){
			qy=(0.5+jy)*delqy;
			for(jz=0;jz<nqzmax;jz++){
				qz=(0.5+jz)*delqz;
				alpha=0.0;
				for(i=0;i<n;i++){
					r=list->GetR(i);
					r2=r[1]*r[1]+r[2]*r[2]+r[3]*r[3];
					arg=2.0*(r[1]*qx+r[2]*qy+r[3]*qz)/HBARC;
					if(arg>2.0*PI) arg=arg-2.0*PI*floor(arg/(2.0*PI));
					if(arg<0.0) arg=arg+2.0*PI*floor(fabs(arg/(2.0*PI)));
					if(r2<1.0E10)	alpha+=exp(ci*arg);
				}
				c=abs(alpha*conj(alpha))-double(n);
				c=c/(double(n)*double(n-1));
				cf->IncrementElement(0,jx,0,jy,0,jz,norm*c);
			}
		}
	}
}

void S2CF::s2c_gauss(CSourceCalc *sourcecalc,CKernelWF *kernel,C3DArray *cf){
	char message[CLog::CHARLENGTH];
	int imc;
	double root2=sqrt(2.0);
	double Rx=(sourcecalc->spars).getD("Rx",4.0);
	double Ry=(sourcecalc->spars).getD("Ry",4.0);
	double Rz=(sourcecalc->spars).getD("Rz",4.0);
	int NMC=(sourcecalc->spars).getI("NMC",500);
	double lambda=(sourcecalc->spars).getD("lambda",1.0);
	double Xoff=(sourcecalc->spars).getD("Xoff",0.0);
	double Yoff=(sourcecalc->spars).getD("Yoff",0.0);
	double Zoff=(sourcecalc->spars).getD("Zoff",0.0);
	//double Euler_Phi=spars.set("Euler_Phi",0.0);
	//double Euler_Theta=spars.set("Euler_Theta",0.0);
	//double Euler_Psi=spars.set("Euler_Psi",0.0);
	Crandy *randy=sourcecalc->randy;
	double x,y,z,r,q,ctheta,qx,qy,qz,wf2,xarray[3],gauss[2];
	int igauss=2,ix;

	bool IDENTICAL=kernel->GetIDENTICAL();
	bool XSYM=cf->GetXSYM();
	bool YSYM=cf->GetYSYM();
	bool ZSYM=cf->GetZSYM();
	if(IDENTICAL&&!(XSYM&&YSYM&&ZSYM)){
		CLog::Info("S2CF::s2c, kernel says particles are identical, but symmetries are not all even\n");
		snprintf(message,CLog::CHARLENGTH,"XSYM=%d, YSYM=%d, ZSYM=%d\n",int(XSYM),int(YSYM),int(ZSYM));
		CLog::Fatal(message);
	}
	int jsx,jsy,jsz,jx,jy,jz;
	int nsx=2,nsy=2,nsz=2;
	if(XSYM) nsx=1;
	if(YSYM) nsy=1;
	if(ZSYM) nsz=1;
	int nqxmax=cf->GetNXMAX();
	int nqymax=cf->GetNYMAX();
	int nqzmax=cf->GetNZMAX();
	double delqx=cf->GetDELX();
	double delqy=cf->GetDELY();
	double delqz=cf->GetDELZ();
	double norm=lambda/double(NMC);
	cf->ZeroArray();

	for(jsx=0;jsx<nsx;jsx++){
		for(jx=0;jx<nqxmax;jx++){
			qx=(0.5+jx)*delqx;
			if(jsx==1) qx=-qx;

			for(jsy=0;jsy<nsy;jsy++){
				for(jy=0;jy<nqymax;jy++){
					qy=(0.5+jy)*delqy;
					if(jsy==1) qy=-qy;

					for(jsz=0;jsz<nsz;jsz++){
						for(jz=0;jz<nqzmax;jz++){
							qz=(0.5+jz)*delqz;
							if(jsz==1) qz=-qz;
							q=sqrt(qx*qx+qy*qy+qz*qz);

							for(imc=0;imc<NMC;imc++){
								for(ix=0;ix<3;ix++){
									if(igauss==2){
										randy->ran_gauss2(gauss[0],gauss[1]);
										igauss=0;
									}
									xarray[ix]=gauss[igauss]; igauss++;
								}
								x=root2*Rx*xarray[0]+Xoff;
								y=root2*Ry*xarray[1]+Yoff;
								z=root2*Rz*xarray[2]+Zoff;
								r=sqrt(x*x+y*y+z*z);
								wf2=0.0;

								if(XSYM&YSYM&&ZSYM){
									ctheta=(qx*x+qy*y+qz*z)/(q*r);
									wf2+=0.25*kernel->GetPsiSquared(q,r,ctheta);
									if(!IDENTICAL) wf2+=0.25*kernel->GetPsiSquared(q,r,-ctheta);
									ctheta=(-qx*x+qy*y+qz*z)/(q*r);
									wf2+=0.25*kernel->GetPsiSquared(q,r,ctheta);
									if(!IDENTICAL) wf2+=0.25*kernel->GetPsiSquared(q,r,-ctheta);
									ctheta=(qx*x-qy*y+qz*z)/(q*r);
									wf2+=0.25*kernel->GetPsiSquared(q,r,ctheta);
									if(!IDENTICAL) wf2+=0.25*kernel->GetPsiSquared(q,r,-ctheta);																																																				
									ctheta=(qx*x+qy*y-qz*z)/(q*r);
									wf2+=0.25*kernel->GetPsiSquared(q,r,ctheta);
									if(!IDENTICAL) wf2+=0.25*kernel->GetPsiSquared(q,r,-ctheta);
									if(!IDENTICAL) wf2=0.5*wf2;
								}
								else if(!XSYM && YSYM && ZSYM){
									ctheta=(qx*x+qy*y+qz*z)/(q*r);
									wf2+=0.25*kernel->GetPsiSquared(q,r,ctheta);
									ctheta=(qx*x-qy*y+qz*z)/(q*r);
									wf2+=0.25*kernel->GetPsiSquared(q,r,ctheta);
									ctheta=(qx*x+qy*y-qz*z)/(q*r);
									wf2+=0.25*kernel->GetPsiSquared(q,r,ctheta);
									ctheta=(qx*x-qy*y-qz*z)/(q*r);
									wf2+=0.25*kernel->GetPsiSquared(q,r,ctheta);
								}
								else if(XSYM && !YSYM && ZSYM){
									ctheta=(qx*x+qy*y+qz*z)/(q*r);
									wf2+=0.25*kernel->GetPsiSquared(q,r,ctheta);
									ctheta=(-qx*x+qy*y+qz*z)/(q*r);
									wf2+=0.25*kernel->GetPsiSquared(q,r,ctheta);
									ctheta=(qx*x+qy*y-qz*z)/(q*r);
									wf2+=0.25*kernel->GetPsiSquared(q,r,ctheta);
									ctheta=(-qx*x+qy*y-qz*z)/(q*r);
									wf2+=0.25*kernel->GetPsiSquared(q,r,ctheta);
								}
								else if(XSYM && YSYM && !ZSYM){
									ctheta=(qx*x+qy*y+qz*z)/(q*r);
									wf2+=0.25*kernel->GetPsiSquared(q,r,ctheta);
									ctheta=(-qx*x+qy*y+qz*z)/(q*r);
									wf2+=0.25*kernel->GetPsiSquared(q,r,ctheta);
									ctheta=(qx*x-qy*y+qz*z)/(q*r);
									wf2+=0.25*kernel->GetPsiSquared(q,r,ctheta);
									ctheta=(-qx*x-qy*y+qz*z)/(q*r);
									wf2+=0.25*kernel->GetPsiSquared(q,r,ctheta);
								}
								else if(XSYM && !YSYM && !ZSYM){
									ctheta=(qx*x+qy*y+qz*z)/(q*r);
									wf2+=0.5*kernel->GetPsiSquared(q,r,ctheta);
									ctheta=(-qx*x+qy*y+qz*z)/(q*r);
									wf2+=0.5*kernel->GetPsiSquared(q,r,ctheta);
								}
								else if(!XSYM && YSYM && !ZSYM){
									ctheta=(qx*x+qy*y+qz*z)/(q*r);
									wf2+=0.5*kernel->GetPsiSquared(q,r,ctheta);
									ctheta=(qx*x-qy*y+qz*z)/(q*r);
									wf2+=0.5*kernel->GetPsiSquared(q,r,ctheta);
								}
								else if(!XSYM && !YSYM && ZSYM){
									ctheta=(qx*x+qy*y+qz*z)/(q*r);
									wf2+=0.5*kernel->GetPsiSquared(q,r,ctheta);
									ctheta=(qx*x+qy*y-qz*z)/(q*r);
									wf2+=0.5*kernel->GetPsiSquared(q,r,ctheta);
								}
								else{
									ctheta=(qx*x+qy*y+qz*z)/(q*r);
									wf2+=kernel->GetPsiSquared(q,r,ctheta);
								}
								cf->IncrementElement(jsx,jx,jsy,jy,jsz,jz,
									norm*(wf2-1.0));
							}
						}
					}
				}
			}
		}
	}
}

void S2CF::s2c_bowlersinyukov(CSourceCalc *sourcecalc,CKernel *kernel,C3DArray *cf3d){
	char message[CLog::CHARLENGTH];
	double Rx=(sourcecalc->spars).getD("Rx",4.0);
	double Ry=(sourcecalc->spars).getD("Ry",4.0);
	double Rz=(sourcecalc->spars).getD("Rz",4.0);
	double HBARC2=HBARC*HBARC;
	double lambda=(sourcecalc->spars).getD("lambda",1.0);
	//double Rbar=pow(Rx*Ry*Rz,1.0/3.0);
	double Rbar=5.0;
	double snorm=1.0/pow(4.0*PI*Rbar*Rbar,1.5);
	double nrmax=kernel->GetNRMAX();
	double nqmax=kernel->GetNQMAX();
	double delr=kernel->GetDELR();
	double delq=kernel->GetDELQ();
	int nqxmax=cf3d->GetNXMAX();
	int nqymax=cf3d->GetNYMAX();
	int nqzmax=cf3d->GetNZMAX();
	double delqx=cf3d->GetDELX();
	double delqy=cf3d->GetDELY();
	double delqz=cf3d->GetDELZ();
	
	double r,ss;
	CCHArray *sf=new CCHArray(0,nrmax,delr,true,true,true);
	double norm=0.0;
	for(int ir=0;ir<nrmax;ir++){
		r=(0.5+ir)*delr;
		ss=snorm*exp(-0.25*r*r/(Rbar*Rbar));
		sf->SetElement(0,0,0,ir,ss);
		norm+=ss*4.0*PI*r*r*delr;
	}
	snprintf(message,CLog::CHARLENGTH,"norm=%g\n",norm);
	CLog::Info(message);
	CCHArray *cf=new CCHArray(0,nqmax,delq,true,true,true);
	S2CF::s2c(sf,kernel,cf);
	//cf->Print(0,0,0);
	//Misc::Pause();

	int jx,jy,jz,iq,iq1,iq2;
	double qx,qy,qz,q,cfbar,w1,w2,q2,bowlsin;	
	for(jx=0;jx<nqxmax;jx++){
		qx=(0.5+jx)*delqx;
		for(jy=0;jy<nqymax;jy++){
			qy=(0.5+jy)*delqy;
			for(jz=0;jz<nqzmax;jz++){
				qz=(0.5+jz)*delqz;
				q=sqrt(qx*qx+qy*qy+qz*qz);
				iq=lrint(q/delq);
				iq1=iq-1;
				iq2=iq;
				if(iq1<0){
					iq1+=1;
					iq2+=1;
				}
				q2=(0.5+iq2)*delq;
				w1=(q2-q)/delq;
				w2=1.0-w1;
				cfbar=w1*cf->GetElement(0,0,0,iq1)+w2*cf->GetElement(0,0,0,iq2);
				bowlsin=lambda*(1.0+exp(-4.0*qx*qx*Rx*Rx/HBARC2-4.0*qy*qy*Ry*Ry/HBARC2-4.0* qz*qz*Rz*Rz/HBARC2))*(1.0+cfbar);
				bowlsin-=lambda;
				cf3d->SetElement(0,jx,0,jy,0,jz,bowlsin);
			}
		}
	}
	delete(sf);
	delete(cf);
}
