#include "msu_coral/sourcecalc.h"
#include "msu_commonutils/randy.h"
#include "msu_commonutils/constants.h"

using namespace std;

CSourceCalc_Blast::CSourceCalc_Blast(){
	InitSPars();
	randy=new Crandy(-1234);
}

void CSourceCalc_Blast::InitSPars(){
	// DEFAULT VALUES
	spars.set("lambda",0.5);
	spars.set("R",13.0);
	spars.set("Tau",12.0);
	spars.set("DelTau",3.0);
	spars.set("Beta",0.7);
	spars.set("T",110.0);
	spars.set("Pt",600.0);
	spars.set("Phi",0.0);
	spars.set("EtaG",1.5);  
	spars.set("Ma",938.28);
	spars.set("Mb",139.58);
	spars.set("Nsample",1000);
}

void CSourceCalc_Blast::SetSPars(double lambdaset,double Rset,double Tauset,double DelTauset,double Betaset,double Tset,double Ptset,double EtaGset,double Maset,double Mbset){
	InitSPars();
	spars.set("lambda",lambdaset);
	spars.set("R",Rset);
	spars.set("Tau",Tauset);
	spars.set("Beta",Betaset);
	spars.set("DelTau",DelTauset);
	spars.set("T",Tset);
	spars.set("Pt",Ptset);
	spars.set("EtaG",EtaGset);  
	spars.set("Ma",Maset);  
	spars.set("Mb",Mbset);
}

void CSourceCalc_Blast::SetSPars(double lambdaset,double Rset,double Tauset,double DelTauset){
	InitSPars();
	spars.set("lambda",lambdaset);
	spars.set("R",Rset);
	spars.set("Tau",Tauset);
	spars.set("DelTau",DelTauset);
}


void CSourceCalc_Blast::GetMCList(double *p,CMCList *mclist){
	double eta,etaG,u[4],x,y,z,t,tt,R,weight;
	double umax,betamax,eprime,eu,eumin,T,tau,tau0,deltau;
	double pt,gamma,gammav;
	int imc,nsample;
	double m=sqrt(p[0]*p[0]-p[1]*p[1]);
	etaG=spars.getD("EtaG",-999);
	R=spars.getD("R",-999);
	T=spars.getD("T",-999);
	tau0=spars.getD("Tau",-999);
	betamax=spars.getD("Beta",-999);
	deltau=spars.getD("DelTau",-999);
	nsample=spars.getI("Nsample",-999);
	umax=betamax/sqrt(1.0-betamax*betamax);
	pt=fabs(p[1]);
	gammav=pt/m;
	gamma=sqrt(1.0+gammav*gammav);
	if(gammav<umax) eumin=m;
	else eumin=sqrt(1.0+umax*umax)*p[0]-umax*p[1];

	for(imc=0;imc<nsample;imc++){
		do{
			eta=etaG*randy->ran_gauss();
			TRY_AGAIN:
			x=(1.0-2.0*randy->ran());
			y=(1.0-2.0*randy->ran());
			if(x*x+y*y>1.0) goto TRY_AGAIN;
			u[1]=umax*x;
			u[2]=umax*y;
			x=x*R;
			y=y*R;
			u[3]=sinh(eta);
			u[0]=sqrt(1.0+u[1]*u[1]+u[3]*u[3]);

			eprime=p[0]*cosh(eta);
			eu=u[0]*p[0]-u[1]*p[1];

			weight=(eprime/p[0])*exp(-(eu-eumin)/T);
			if(weight>1.0){
				printf("DISASTER! weight=%g which is > 1.0, eu=%g, eumin=%g\n",weight,eu,eumin);
				exit(1);
			}
		} while(weight<randy->ran());

		tau=GetTau(tau0,deltau);
		z=tau*sinh(eta);
		t=tau*cosh(eta);

		tt=gamma*t-gammav*x;
		x=gamma*x-gammav*t;
		t=tt;
		mclist->SetR(imc,t,x,y,z);
	}
}

double CSourceCalc_Blast::GetTau(double tau0,double deltau){
	double tau;
	do{
		tau=tau0+deltau*randy->ran_gauss();
	}while(tau<0.0);
	return tau;
}

void CSourceCalc_Blast::CalcS(CCHArray *A){
	int nsample;
	double rcm[4],*ra,*rb;
	double delr,ma,mb,volume;
	int ir,ia,ib,nbmax,nrmax;
	double x,y,z,ex,ey,ez,snorm,r,lambda;
	bool SameMass;
	CMCList *lista,*listb;
	delr=A->GetRADSTEP();
	nrmax=A->GetNRADIAL();

	ma=spars.getD("Ma",-999);
	mb=spars.getD("Mb",-999);
	nsample=spars.getI("Nsample",1000);
	lambda=spars.getD("lambda",-999);
	SameMass=0;
	if(fabs(ma-mb)<1.0) SameMass=1;

	lista=new CMCList(nsample);
	if(!SameMass) listb=new CMCList(nsample);
	else listb=lista;
	CalcS(lista,listb);

	rcm[0]=0.0;
	for(ia=0;ia<nsample;ia++){
		nbmax=nsample;
		if(SameMass) nbmax=ia-1;
		for(ib=0;ib<nbmax;ib++){
			ra=lista->GetR(ia);
			rb=listb->GetR(ib);

			rcm[1]=ra[1]-rb[1];
			rcm[2]=ra[2]-rb[2];
			rcm[3]=ra[3]-rb[3];
			r=sqrt(rcm[1]*rcm[1]+rcm[2]*rcm[2]+rcm[3]*rcm[3]);
			x=rcm[1];
			y=rcm[2];
			z=rcm[3];
			ir=int(floor(r/delr));
			if(ir<nrmax){
				ex=x/r; ey=y/r; ez=z/r;
				A->IncrementAExpArrayFromE(ex,ey,ez,lambda,ir);
			}
		}
	}
	A->FillRemainderX();
	if(SameMass)
		snorm=2.0/double(nsample*(nsample-1));
	else
		snorm=1.0/double(nsample*nsample);

	for(ir=0;ir<nrmax;ir++){
		volume=(4.0*PI/3)*(pow((ir+1)*delr,3)
			-pow(double(ir)*delr,3));
		A->ScaleArray(snorm/volume,ir);
	}

	delete lista;
	if(!SameMass) delete listb;

}

void CSourceCalc_Blast::CalcS(CMCList *lista,CMCList *listb){
	double pa[4],pb[4];
	double ma,mb,Pt,lambda;
	int nsample;

	lambda=spars.getD("lambda",1.0);
	lista->SetNorm(sqrt(lambda));
	listb->SetNorm(sqrt(lambda));
	nsample=spars.getI("Nsample",1000);
	ma=spars.getD("Ma",-999);
	mb=spars.getD("Mb",-999);
	pa[3]=pb[3]=0.0;
	Pt=spars.getD("Pt",-999);
	pa[1]=Pt*ma/(ma+mb);
	pb[1]=Pt*mb/(ma+mb);
	pa[2]=pa[3]=0.0;
	pb[2]=pb[3]=0.0;
	pa[0]=sqrt(ma*ma+pa[1]*pa[1]);
	pb[0]=sqrt(mb*mb+pb[1]*pb[1]);

	if(lista->GetNMC()!=nsample) lista->Resize(nsample);
	GetMCList(pa,lista);

	if(lista!=listb){
		if(listb->GetNMC()!=nsample) listb->Resize(nsample);
		GetMCList(pb,listb);
	}

}
