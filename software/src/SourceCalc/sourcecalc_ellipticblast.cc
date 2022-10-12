#include "msu_coral/sourcecalc.h"
#include "msu_commonutils/randy.h"

using namespace std;

CSourceCalc_EllipticBlast::CSourceCalc_EllipticBlast(){
  InitSPars();
  randy=new Crandy(-1234);
}

void CSourceCalc_EllipticBlast::InitSPars(){
  // DEFAULT VALUES
  spars.set("Rx",13);
  spars.set("Ry",13);
  spars.set("Tau",12);
  spars.set("BetaX",0.7);
  spars.set("BetaY",0.7);
  spars.set("T",110);
  spars.set("Pt",600);
  spars.set("Phi",0.0);
  spars.set("EtaG",2.0);  
  spars.set("Ma",938.28);
  spars.set("Mb",139.58);
  spars.set("Nsample",1000);
}

void CSourceCalc_EllipticBlast::SetSPars(double Rxset,double Ryset,double Tauset,
				double BetaXset,double BetaYset,
				double Tset,double Ptset,
				double Phiset,double EtaGset,
				double Maset,double Mbset){
  InitSPars();
  spars.set("Rx",Rxset);
  spars.set("Ry",Ryset);
  spars.set("Tau",Tauset);
  spars.set("BetaX",BetaXset);
  spars.set("BetaY",BetaYset);
  spars.set("T",Tset);
  spars.set("Pt",Ptset);
  spars.set("Phi",Phiset);
  spars.set("EtaG",EtaGset);  
  spars.set("Ma",Maset);  
  spars.set("Mb",Mbset);

}

void CSourceCalc_EllipticBlast::SetSPars(double Rset,double Tauset,
				    double Betaset,double Tset,double Ptset){

  InitSPars();
  spars.set("Rx",Rset);
  spars.set("Ry",Rset);
  spars.set("Tau",Tauset);
  spars.set("BetaX",Betaset);
  spars.set("BetaY",Betaset);
  spars.set("T",Tset);
  spars.set("Pt",Ptset);

}

void CSourceCalc_EllipticBlast::CalcS(CCHArray *A){
  double pa[4],pb[4];
  double cphi,sphi,ma,mb,phi,Pt,volume,delr;
  int nsample;
  double rcm[4],**ra,**rb;
  int ir,imc,ia,ib,nbmax,nrmax;
  double x,y,z,xbar,ybar,zbar,x2bar,y2bar,z2bar,ex,ey,ez,snorm,r;
  bool SameMass;
  const double PI=4.0*atan(1.0);
  delr=A->GetRADSTEP();
  nrmax=A->GetNRADIAL();

  nsample=spars.getI("Nsample",1000);
  ma=spars.getD("Ma",-999);
  mb=spars.getD("Mb",-999);
  SameMass=0;
  if(fabs(ma-mb)<1.0) SameMass=1;
  pa[3]=pb[3]=0.0;
  phi=spars.getD("Phi",-999);
  Pt=spars.getD("Pt",-999);
  sphi=sin(phi);
  cphi=cos(phi);
  pa[1]=cphi*Pt*ma/(ma+mb);
  pb[1]=cphi*Pt*mb/(ma+mb);
  pa[2]=sphi*Pt*ma/(ma+mb);
  pb[2]=sphi*Pt*mb/(ma+mb);
  pa[0]=sqrt(ma*ma+pa[1]*pa[1]+pa[2]*pa[2]+pa[3]*pa[3]);
  pb[0]=sqrt(mb*mb+pb[1]*pb[1]+pb[2]*pb[2]+pb[3]*pb[3]);

  ra=new double *[nsample];
  for(imc=0;imc<nsample;imc++) ra[imc]=new double[4];
  Get_r(pa,nsample,ra);
  if(SameMass){
    rb=ra;
  }
  else{
    rb=new double *[nsample];
    for(imc=0;imc<nsample;imc++) rb[imc]=new double[4];
    Get_r(pb,nsample,rb);
  }

  rcm[0]=0.0;
  xbar=ybar=zbar=x2bar=y2bar=z2bar=0.0;
  for(ia=0;ia<nsample;ia++){
    nbmax=nsample;
    if(SameMass) nbmax=ia-1;
    for(ib=0;ib<nbmax;ib++){
      rcm[1]=ra[ia][1]-rb[ib][1];
      rcm[2]=ra[ia][2]-rb[ib][2];
      rcm[3]=ra[ia][3]-rb[ib][3];

      r=sqrt(rcm[1]*rcm[1]+rcm[2]*rcm[2]+rcm[3]*rcm[3]);
      x=rcm[1];
      y=rcm[2];
      z=rcm[3];
      xbar+=x;
      x2bar+=x*x;
      ybar+=y;
      y2bar+=y*y;
      zbar+=z;
      z2bar+=z*z;

      ir=int(floor(r/delr));
      if(ir<nrmax){
	ex=x/r; ey=y/r; ez=z/r;
	A->IncrementAExpArrayFromE(ex,ey,ez,1.0,ir);
      }
    }
    if(10*(ia+1)%nsample==0){
      sprintf(message,"finished %g percent\n",100*double(ia+1)/double(nsample));
		CLog::Fatal(message);
	}
  }
  A->FillRemainderX();
  if(SameMass) snorm=2.0/double(nsample*(nsample-1));
  else snorm=1.0/double(nsample*nsample);
  xbar*=snorm;
  x2bar*=snorm;
  ybar*=snorm;
  y2bar*=snorm;
  zbar*=snorm;
  z2bar*=snorm;
  x2bar=x2bar-xbar*xbar;
  y2bar=y2bar-ybar*ybar;
  z2bar=z2bar-zbar*zbar;

  sprintf(message,"xbar=%g, ybar=%g, zbar=%g\n",xbar,ybar,zbar);
  CLog::Info(message);
  sprintf(message,"Effective Gaussian Radii: Rout=%g, Rside=%g, Rlong=%g\n",
  sqrt(0.5*x2bar),sqrt(0.5*y2bar),sqrt(0.5*z2bar));
  CLog::Info(message);

  for(ir=0;ir<nrmax;ir++){
    volume=(4.0*PI/3)*(pow((ir+1)*delr,3)
		       -pow(double(ir)*delr,3));
    A->ScaleArray(snorm/volume,ir);
  }

  for(imc=0;imc<nsample;imc++){
    delete ra[imc];
    if(!SameMass) delete rb[imc];
  }
  delete ra;
  delete rb;

}

void CSourceCalc_EllipticBlast::Get_r(double *p,int nsample,double **r){
  double eta,etaG,u[4],x,y,z,t,Rx,Ry,weight;
  double uxmax,uymax,betaxmax,betaymax,eprime,eu;
  double pt,gamma,gammav,rap,rout,rlong,rside,sinhy,coshy,T,tau;
  int imc;
  double m=sqrt(p[0]*p[0]-p[1]*p[1]-p[2]*p[2]-p[3]*p[3]);
  etaG=spars.getD("EtaG",-999);
  Rx=spars.getD("Rx",-999);
  Ry=spars.getD("Ry",-999);
  betaxmax=spars.getD("BetaX",-999);
  betaymax=spars.getD("BetaY",-999);
  T=spars.getD("T",-999);
  tau=spars.getD("tau",-999);
  uxmax=betaxmax/sqrt(1.0-betaxmax*betaxmax);
  uymax=betaymax/sqrt(1.0-betaymax*betaymax);
  rap=atanh(p[3]/p[0]);
  pt=sqrt(p[1]*p[1]+p[2]*p[2]);
  gammav=pt/m;
  gamma=sqrt(1.0+gammav*gammav);

  for(imc=0;imc<nsample;imc++){
	  do{
		  eta=etaG*randy->ran_gauss();
		  TRY_AGAIN:
		  x=(1.0-2.0*randy->ran());
		  y=(1.0-2.0*randy->ran());
		  if(x*x+y*y>1.0) goto TRY_AGAIN;
		  u[1]=uxmax*x;
		  u[2]=uymax*y;
		  x=x*Rx;
		  y=y*Ry;
		  u[3]=sinh(eta);
		  u[0]=sqrt(1.0+u[1]*u[1]+u[2]*u[2]+u[3]*u[3]);
      
		  eprime=p[0]*cosh(eta)-p[3]*u[3];
		  eu=u[0]*p[0]-u[1]*p[1]-u[2]*p[2]-u[3]*p[3];

		  weight=(eprime/p[0])*exp(-(eu-m)/T);
		  if(weight>1.0){
			  CLog::Fatal("DISASTER! weight > 1.0\n");
		  }
	  } while(weight<randy->ran());

	  z=tau*sinh(eta);
	  t=tau*cosh(eta);

	  rout=(p[1]*x+p[2]*y)/pt;
	  rside=(p[1]*y-p[2]*x)/pt;
	  sinhy=sinh(rap);
	  coshy=cosh(rap);
	  rlong=coshy*z-sinhy*t;
	  t=coshy*t-sinhy*z;
	  r[imc][0]=gamma*t-gammav*rout;
	  r[imc][1]=gamma*rout-gammav*t;
	  r[imc][2]=rside;
	  r[imc][3]=rlong;

  }
}
