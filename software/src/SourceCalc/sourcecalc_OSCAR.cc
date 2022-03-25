#ifndef __INCLUDE_SOURCECALC_OSCAR__
#define __INCLUDE_SOURCECALC_OSCAR__
#include "sourcecalc.h"
#include "msu_commonutils/randy.h"
#include "msu_commonutils/constants.h"

using namespace std;

CSourceCalc_OSCAR::CSourceCalc_OSCAR(){
	InitSPars();	
	spars.PrintPars();
	randy=new Crandy(1234);
}

CSourceCalc_OSCAR::CSourceCalc_OSCAR(string sparsfilename){
	InitSPars();
	spars.ReadParsFromFile(sparsfilename);
	spars.PrintPars();
	randy=new Crandy(1234);
}

void CSourceCalc_OSCAR::InitSPars(){
	// DEFAULT VALUES
	spars.set("PT",600.0);
	spars.set("DELPT",20.0);
	spars.set("PHIMIN_DEG",0.0);
	spars.set("PHIMAX_DEG",360.0);
	spars.set("YMIN",-1.0);
	spars.set("YMAX",1.0);
	spars.set("NPARTSMAX",20000);
	spars.set("OSCARfilename",string("UNDEFINED"));
	spars.set("NEVENTSMAX",10000);
	spars.set("ETA_GAUSS",1.2);
}

void CSourceCalc_OSCAR::SetSPars(double PT_set,double DELPT_set,double PHIMIN_DEG_set,double PHIMAX_DEG_set,double YMIN_set,double YMAX_set){
	spars.set("PT",PT_set);
	spars.set("DELPT",DELPT_set);
	spars.set("PHIMIN_DEG",PHIMIN_DEG_set);
	spars.set("PHIMAX_DEG",PHIMAX_DEG_set);
	spars.set("YMIN",YMIN_set);
	spars.set("YMAX",YMAX_set);
}

void CSourceCalc_OSCAR::SetIDs(int *idlista,int nida,int *idlistb,int nidb){
	idlist_a=idlista;
	nid_a=nida;
	nid_b=nidb;
	idlist_b=idlistb;
}

void CSourceCalc_OSCAR::CalcS(CMCList *&lista,CMCList *&listb){
	double **ra,**rb;
	int ia,ib,na=0,nb=0;
	bool AEQUALB=spars.getB("AEQUALB",true);
	int NPARTSMAX=spars.getI("NPARTSMAX",50000);

	ra=new double *[NPARTSMAX];
	for(ia=0;ia<NPARTSMAX;ia++) ra[ia]=new double[4];
	if(AEQUALB){
		rb=ra;
	}
	else{
		rb=new double *[NPARTSMAX];
		for(ib=0;ib<NPARTSMAX;ib++) rb[ib]=new double[4];
	}
	ReadR(ra,na,rb,nb);

	if(lista==NULL) lista=new CMCList(na);
	else if(lista->GetNMC()!=na) lista->Resize(na);
	if(!AEQUALB){
		if(listb==NULL) listb=new CMCList(nb);
		else if(listb->GetNMC()!=nb) listb->Resize(nb);
	}

	for(int i=0;i<na;i++) lista->SetR(i,ra[i]);
	if(lista!=listb){
		for(int i=0;i<nb;i++) listb->SetR(i,rb[i]);
	}

	for(ia=0;ia<NPARTSMAX;ia++){
		delete [] ra[ia];
	}
	delete [] ra;
	if(!AEQUALB){
		for(ib=0;ib<NPARTSMAX;ib++){
			delete [] rb[ib];
		}
		delete [] rb;
	}
	printf("______ FINISHED CREATING MCLISTS ___________\n");

}

void CSourceCalc_OSCAR::ReadR(double **ra,int &na,double **rb,int &nb){
	FILE *oscarfile;
	double r[4],p[4];
	double MA=spars.getD("MA",139.57);
	double MB=spars.getD("MB",139.57);
	double mass,rdummy1,rdummy2;
	int ident,idummy,i,alpha;
	int ndummy_header=3,ndummy_betweenevents=0;
	int npart,npartmax,ipart,ievent;
	int NEVENTSMAX=spars.getI("NEVENTSMAX",100);
	string OSCARfilename=spars.getS("OSCARfilename","UNDEFINED");
	bool AEQUALB=spars.getB("AEQUALB",false);
	char dummy[160];
	na=nb=0;
	printf("Opening %s\n",OSCARfilename.c_str());
	oscarfile=fopen(OSCARfilename.c_str(),"r");
	// Read Header Info
	for(i=0;i<ndummy_header;i++){
		fgets(dummy,120,oscarfile);
	}

	ievent=0;
	do{
		ievent+=1;
		fscanf(oscarfile,"%d %d %lf %lf",&idummy,&npartmax,&rdummy1,&rdummy2);
		fgets(dummy,120,oscarfile);
		for(npart=0;npart<npartmax;npart++){
			fscanf(oscarfile,"%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf",
				&ipart,&ident,&p[1],&p[2],&p[3],&p[0],&mass,&r[1],&r[2],&r[3],&r[0]);
			for(alpha=0;alpha<4;alpha++)
				p[alpha]*=1000.0;
			if(IDMatch(ident,idlist_a,nid_a))
				Check(p,r,MA,ra,na);
			if(!AEQUALB) if(IDMatch(ident,idlist_b,nid_b)) Check(p,r,MB,rb,nb);
		}
		for(i=0;i<ndummy_betweenevents;i++){
			if(!feof(oscarfile)) fgets(dummy,120,oscarfile);
		}
	} while(ievent<NEVENTSMAX && !feof(oscarfile));
	fclose(oscarfile);
	if(AEQUALB) nb=na;
	printf("OSCAR file read: %d events, na=%d, nb=%d\n",ievent,na,nb);
}

bool CSourceCalc_OSCAR::Check(double *p,double *r,double m,double **ra,int &n){
	double PT=spars.getD("PT",600.0);
	double DELPT=spars.getD("DELPT",20.0);
	double YMIN=spars.getD("YMIN",-1.0);
	double YMAX=spars.getD("YMAX",1.0);
	double PHIMIN=2.0*PI*spars.getD("PHIMIN_DEG",0.0)/360.0;
	double PHIMAX=2.0*PI*spars.getD("PHIMAX_DEG",0.0)/360.0;
	double ETA_GAUSS=spars.getD("ETA_GAUSS",1.2);
	double MA=spars.getD("MA",139.57);
	double MB=spars.getD("MB",139.57);
	double pt,phi,eta;
	double gammav=PT/(MA+MB);
	double gamma=sqrt(1.0+gammav*gammav);
	double pttarget=PT*m/(MA+MB);
	double rout,rlong,rside,sinhy,coshy,tau,vperp,y;
	const double TAUCOMPARE=12.0;
	int NPARTSMAX=spars.getI("NPARTSMAX",50000);
	bool success=false;
	
	pt=sqrt(p[1]*p[1]+p[2]*p[2]);
	if(fabs(pt-pttarget)<DELPT){
		if(p[1]!=p[1] || p[2]!=p[2] || p[3]!=p[3]){
			printf("bad particle has nan, p=(%g,%g,%g)\n",p[1],p[2],p[3]);
			return false;
		}
		phi=atan2(p[1],p[2]);
		if( ((PHIMIN<PHIMAX)&&(phi>PHIMIN && phi<PHIMAX)) 
		|| ((PHIMIN>PHIMAX)&&(phi>PHIMIN||phi<PHIMAX)) ){
			y=atanh(p[3]/p[0]);
			if(y>YMIN && y<YMAX){
				p[0]=sqrt(pt*pt+p[3]*p[3]+m*m);
				rout=(p[1]*r[1]+p[2]*r[2])/pt;
				rside=(p[1]*r[2]-p[2]*r[1])/pt;
				sinhy=sinh(y);
				coshy=cosh(y);
				rlong=coshy*r[3]-sinhy*r[0];
				tau=coshy*r[0]-sinhy*r[3];
				eta=asinh(rlong/tau);
				if(randy->ran()<exp(-0.5*eta*eta/(ETA_GAUSS*ETA_GAUSS))){
				//printf("%g %g %g %g\n",tau,rout,rside,rlong);
					vperp=pt/sqrt(m*m+pt*pt);
					rout=rout-vperp*(tau-TAUCOMPARE);
					tau=TAUCOMPARE;
					rout=gamma*rout-gammav*tau;
					ra[n][0]=0.0;
					ra[n][1]=rout;
					ra[n][2]=rside;
					ra[n][3]=rlong;
					n+=1;
					success=true;
					if(n==NPARTSMAX){
						printf("TOO MANY PARTICLES FIT CRITERIA, increase parameter NPARTSMAX=%d\n",NPARTSMAX);
						exit(1);
					}
				}
			}
		}
	}
	return success;
}

bool CSourceCalc_OSCAR::IDMatch(int ident,int *idlist,int nid){
	int i=0;
	bool answer=false;
	while(answer==false && i<nid){
		if(ident==idlist[i]) answer=true;
		i+=1;
	}
	if(i==nid && answer==false){
		printf("IDMATCH failed in CSourceCalc_OSCAR::IDMatch !!!\n");
		exit(1);
	}
	return answer;
}

#endif
