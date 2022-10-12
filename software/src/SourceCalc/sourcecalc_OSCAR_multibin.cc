#include "msu_coral/sourcecalc.h"
#include "msu_commonutils/randy.h"
//#include "part.h"
#include "msu_commonutils/constants.h"

using namespace std;

CSourceCalc_OSCAR_MultiBin::CSourceCalc_OSCAR_MultiBin(){
	InitSPars();
	randy=new Crandy(1234);
	B3D_BINARY_FORMAT=false;
}

CSourceCalc_OSCAR_MultiBin::CSourceCalc_OSCAR_MultiBin(string sparsfilename){
	InitSPars();
	spars.ReadParsFromFile(sparsfilename);
	NPHIBINS=spars.getI("NPHIBINS",6);
	NPTBINS=spars.getI("NPTBINS",24);
	DELPT=spars.getD("DELPT",25.0);
	PTMIN=spars.getD("PTMIN",150.0);
	OSCARfilename=spars.getS("OSCARfilename","OSCARfilename_undefined");
	PTMAX=PTMIN+NPTBINS*DELPT;
	DELPHI=90.0/double(NPHIBINS);
	spars.PrintPars();
	randy=new Crandy(1234);
	B3D_BINARY_FORMAT=false;
}

void CSourceCalc_OSCAR_MultiBin::InitSPars(){
	// DEFAULT VALUES
	spars.set("DELPT",25.0);
	spars.set("NPTBINS",24);
	spars.set("NPHIBINS",6);
	spars.set("PTMIN",150.0);
	spars.set("PHIMIN_DEG",0.0);
	spars.set("PHIMAX_DEG",360.0);
	spars.set("YMIN",-1.0);
	spars.set("YMAX",1.0);
	spars.set("NPARTSMAX",10000);
	spars.set("NEVENTSMAX",10000);
	spars.set("ETA_GAUSS",1.2);
	NPHIBINS=spars.getI("NPHIBINS",6);
	NPTBINS=spars.getI("NPTBINS",24);
	DELPT=spars.getD("DELPT",25.0);
	PTMIN=spars.getD("PTMIN",150.0);
	OSCARfilename=spars.getS("OSCARfilename","OSCARfilename_undefined");
	B3D_BINARY_FORMAT=spars.getB("B3D_BINARY_FORMAT",false);
	PTMAX=PTMIN+NPTBINS*DELPT;
	DELPHI=90.0/double(NPHIBINS);
}

void CSourceCalc_OSCAR_MultiBin::SetSPars(double PT_set,double DELPT_set,double PHIMIN_DEG_set,double PHIMAX_DEG_set,double YMIN_set,double YMAX_set){
	spars.set("PT",PT_set);
	spars.set("DELPT",DELPT_set);
	spars.set("PHIMIN_DEG",PHIMIN_DEG_set);
	spars.set("PHIMAX_DEG",PHIMAX_DEG_set);
	spars.set("YMIN",YMIN_set);
	spars.set("YMAX",YMAX_set);
}

void CSourceCalc_OSCAR_MultiBin::SetIDs(int *idlista,int nida,int *idlistb,int nidb){
	idlist_a=idlista;
	nid_a=nida;
	nid_b=nidb;
	idlist_b=idlistb;
}

void CSourceCalc_OSCAR_MultiBin::CalcS(CMCPRList ***&lista,CMCPRList ***&listb){
	double ****ra,****rb,****pa,****pb;
	int ia,ib,**na,**nb,ipt,iphi;
	bool AEQUALB=spars.getB("AEQUALB",true);
	int NPARTSMAX=spars.getI("NPARTSMAX",10000);
	//spars.PrintPars();
	
	if(lista!=NULL){
		for(ipt=0;ipt<NPTBINS;ipt++){
			for(iphi=0;iphi<NPHIBINS;iphi++){
				if(lista[ipt][iphi]!=NULL){
					delete lista[ipt][iphi];
					lista[ipt][iphi]=NULL;
				}
				delete [] lista[ipt];
			}
			delete [] lista;
		}
	}
	lista=NULL;
	if(!AEQUALB && listb!=NULL){
		for(ipt=0;ipt<NPTBINS;ipt++){
			for(iphi=0;iphi<NPHIBINS;iphi++){
				if(listb[ipt][iphi]!=NULL){
					delete listb[ipt][iphi];
					listb[ipt][iphi]=NULL;
				}
				delete [] listb[ipt];
			}
			delete [] listb;
		}
		listb=NULL;
	}
	
	lista=new CMCPRList**[NPTBINS];
	for(ipt=0;ipt<NPTBINS;ipt++){
		lista[ipt]=new CMCPRList*[NPHIBINS];
		for(iphi=0;iphi<NPHIBINS;iphi++) lista[ipt][iphi]=NULL;
	}
	if(!AEQUALB){
		listb=new CMCPRList**[NPTBINS];
		for(ipt=0;ipt<NPTBINS;ipt++){
			listb[ipt]=new CMCPRList*[NPHIBINS];
			for(iphi=0;iphi<NPHIBINS;iphi++) listb[ipt][iphi]=NULL;
		}
	}
	
	ra=new double ***[NPTBINS];
	pa=new double ***[NPTBINS];
	na=new int *[NPTBINS];
	for(ipt=0;ipt<NPTBINS;ipt++){
		ra[ipt]=new double **[NPHIBINS];
		pa[ipt]=new double **[NPHIBINS];
		na[ipt]=new int[NPHIBINS];
		for(iphi=0;iphi<NPHIBINS;iphi++){
			ra[ipt][iphi]=new double *[NPARTSMAX];
			pa[ipt][iphi]=new double*[NPARTSMAX];
			for(ia=0;ia<NPARTSMAX;ia++){
				ra[ipt][iphi][ia]=new double[4];
				pa[ipt][iphi][ia]=new double[4];
			}
		}
	}
	
	
	if(AEQUALB){
		rb=ra;
		pb=pa;
		nb=na;
	}
	else{
		rb=new double ***[NPTBINS];
		pb=new double ***[NPTBINS];
		nb=new int *[NPTBINS];
		for(ipt=0;ipt<NPTBINS;ipt++){
			rb[ipt]=new double **[NPHIBINS];
			pb[ipt]=new double **[NPHIBINS];
			nb[ipt]=new int[NPHIBINS];
			for(iphi=0;iphi<NPHIBINS;iphi++){
				rb[ipt][iphi]=new double *[NPARTSMAX];
				pb[ipt][iphi]=new double *[NPARTSMAX];
				for(ia=0;ia<NPARTSMAX;ia++){
					rb[ipt][iphi][ia]=new double[4];
					pb[ipt][iphi][ia]=new double[4];
				}
			}
		}
	}
	
	ReadPR(pa,ra,na,pb,rb,nb);

	for(ipt=0;ipt<NPTBINS;ipt++){
		for(iphi=0;iphi<NPHIBINS;iphi++){	
			lista[ipt][iphi]=new CMCPRList(na[ipt][iphi]);
			if(!AEQUALB){
				listb[ipt][iphi]=new CMCPRList(nb[iphi][iphi]);
			}
			for(ia=0;ia<na[ipt][iphi];ia++)
				lista[ipt][iphi]->SetPR(ia,pa[ipt][iphi][ia],ra[ipt][iphi][ia]);
			if(lista!=listb){
				for(ib=0;ib<nb[ipt][iphi];ib++)
					listb[ipt][iphi]->SetPR(ib,pb[ipt][iphi][ib],rb[ipt][iphi][ib]);
			}
		}
	}
	
	for(ipt=0;ipt<NPTBINS;ipt++){
		for(iphi=0;iphi<NPHIBINS;iphi++){
			for(ia=0;ia<NPARTSMAX;ia++){
				delete [] ra[ipt][iphi][ia];
				delete [] pa[ipt][iphi][ia];
			}
			delete [] ra[ipt][iphi];
			delete [] pa[ipt][iphi];
		}
		delete [] ra[ipt];
		delete [] pa[ipt];
		delete [] na[ipt];
	}
	delete [] ra;
	delete [] pa;
	delete [] na;
	if(!AEQUALB){
		for(ipt=0;ipt<NPTBINS;ipt++){
			for(iphi=0;iphi<NPHIBINS;iphi++){
				for(ia=0;ia<NPARTSMAX;ia++){
					delete [] rb[ipt][iphi][ia];
					delete [] pb[ipt][iphi][ia];
				}
				delete [] rb[ipt][iphi];
				delete [] pb[ipt][iphi];
			}
			delete [] rb[ipt];
			delete [] pb[ipt];
			delete [] nb[ipt];
		}
	}
}

/*
long long int CSourceCalc_OSCAR_MultiBin::ReadPR(double ****pa,double ****ra,int **&na,double ****pb,double ****rb,int **&nb){
	CB3DBinaryPartInfo bpart;
	long long int nread=0;
	int ipt,iphi,i;
	double r[4],p[4],phi,pt;
	double MA=spars.getD("MA",139.57);
	double MB=spars.getD("MB",139.57);
	double mass,mt,tau,rdummy1,rdummy2,rapidity;
	int ID,idummy;
	double weight;
	int ipart,nparts,nevents;
	int NEVENTSMAX=spars.getI("NEVENTSMAX",100);
	bool AEQUALB=spars.getB("AEQUALB",false);
	int NPARTSMAX=spars.getI("NPARTSMAX",10000);
	FILE *oscarfile;

	if(B3D_BINARY_FORMAT){
		oscarfile=fopen(OSCARfilename.c_str(),"rb");
	}
	else{
		oscarfile=fopen(OSCARfilename.c_str(),"r");
	}
	char dummy[160];
	for(ipt=0;ipt<NPTBINS;ipt++){
		for(iphi=0;iphi<NPHIBINS;iphi++){
			na[ipt][iphi]=0;
			if(!AEQUALB) nb[ipt][iphi]=0;
		}
	}

	if(!B3D_BINARY_FORMAT){
		for(i=0;i<3;i++)
			fgets(dummy,160,oscarfile);
	}
	nevents=1;
	while(nevents<=NEVENTSMAX && !feof(oscarfile)){
		if(B3D_BINARY_FORMAT){
			fread(&idummy,sizeof(int),1,oscarfile);
			fread(&nparts,sizeof(int),1,oscarfile);
		}
		else{
			fscanf(oscarfile,"%d %d %lf %lf",&idummy,&nparts,&rdummy1,&rdummy2);
			fgets(dummy,160,oscarfile);
		}
		if(!feof(oscarfile)){
			for(ipart=0;ipart<nparts;ipart++){
				if(B3D_BINARY_FORMAT){
					fread(&bpart,sizeof(bpart),1,oscarfile);
					nread+=1;
					ID=bpart.ID;
					p[1]=bpart.px;
					p[2]=bpart.py;
					p[3]=0.0;
					rapidity=bpart.rapidity;
					r[1]=bpart.x;
					r[2]=bpart.y;
					tau=bpart.tau;
					r[0]=tau*cosh(bpart.eta-rapidity);
					r[3]=tau*sinh(bpart.eta-rapidity);
					weight=bpart.weight;
				}
				else{
					fscanf(oscarfile,"%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
					&idummy,&ID,&p[1],&p[2],&p[3],&p[0],&mass,&r[1],&r[2],&r[3],&r[0],&weight);
				}
				pt=sqrt(p[1]*p[1]+p[2]*p[2]);
				if(pt>PTMIN && pt<PTMAX){			
					phi=fabs(atan(p[2]/p[1])*180.0/PI);
					iphi=int(lrint(floor(phi/DELPHI)));
					ipt=int(lrint(floor((pt-PTMIN)/DELPT)));
					if(ipt>=NPTBINS || ipt<0){
						sprintf(message,"ipt=%d is out of range\n",ipt);
						CLog::Fatal(message);
					}
					if(iphi<0 || iphi>=NPHIBINS){
						sprintf(message"iphi=%d is out of range, phi=%g\n",iphi,phi);
						bpart.Print();
						CLog::Fatal(message);
					}
					if(IDMatch(ID,idlist_a,nid_a) && na[ipt][iphi]<NPARTSMAX){
						mass=MA;
						mt=sqrt(mass*mass+p[1]*p[1]+p[2]*p[2]);
						p[0]=sqrt(mt*mt+p[3]*p[3]);	
						Check(p,r,MA,pa[ipt][iphi],ra[ipt][iphi],na[ipt][iphi]);
					}
					if(!AEQUALB){
						if(IDMatch(ID,idlist_b,nid_b) && nb[ipt][iphi]<NPARTSMAX){
							mass=MB;
							mt=sqrt(mass*mass+p[1]*p[1]+p[2]*p[2]);
							p[0]=sqrt(mt*mt+p[3]*p[3]);	
							Check(p,r,MB,pb[ipt][iphi],rb[ipt][iphi],nb[ipt][iphi]);
						}
					}
				}
			}
			nevents+=1;
		}
	}
	fclose(oscarfile);
	if(AEQUALB) nb=na;
	return nread;
}*/

bool CSourceCalc_OSCAR_MultiBin::Check(double *p,double *r,double m,double **pa,double **ra,int &n){
	double YMIN=spars.getD("YMIN",-1.0);
	double YMAX=spars.getD("YMAX",1.0);
	double PHIMIN=2.0*PI*spars.getD("PHIMIN_DEG",0.0)/360.0;
	double PHIMAX=2.0*PI*spars.getD("PHIMAX_DEG",0.0)/360.0;
	double ETA_GAUSS=spars.getD("ETA_GAUSS",1.2);
	double MA=spars.getD("MA",139.57);
	double phi,eta;
	double rout,rlong,rside,sinhy,coshy,tau,vperp,y;
	const double TAUCOMPARE=12.0;
	int NPARTSMAX=spars.getI("NPARTSMAX",20000);
	bool XREFLECTIONSYMMETRY=spars.getB("XREFLECTIONSYMMETRY",false);
	bool YREFLECTIONSYMMETRY=spars.getB("YREFLECTIONSYMMETRY",false);
	bool success=false;
	double pt=sqrt(p[1]*p[1]+p[2]*p[2]);
	double gammav=pt/MA;
	double gamma=sqrt(1.0+gammav*gammav);
	if(n<NPARTSMAX){
		if(p[1]!=p[1] || p[2]!=p[2] || p[3]!=p[3]){
			sprintf(message,"bad particle has nan, p=(%g,%g,%g)\n",p[1],p[2],p[3]);
			CLog::Info(message);
			return false;
		}
		if(XREFLECTIONSYMMETRY){
			if(p[1]<0.0){
				p[1]=-p[1];
				r[1]=-r[1];
			}
		}
		if(YREFLECTIONSYMMETRY){
			if(p[2]<0.0){
				p[2]=-p[2];
				r[2]=-r[2];
			}
		}
		phi=atan2(p[1],p[2]);
		if(phi>PHIMIN && phi<PHIMAX){
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
					if(n==NPARTSMAX){
						sprintf(message,"TOO MANY PARTICLES FIT CRITERIA, increase parameter NPARTSMAX=%d if you want more\n",NPARTSMAX);
						CLog::Fatal(message);
					}
					vperp=pt/sqrt(m*m+pt*pt);
					rout=rout-vperp*(tau-TAUCOMPARE);
					tau=TAUCOMPARE;
					rout=gamma*rout-gammav*tau;
					ra[n][0]=0.0;
					ra[n][1]=rout;
					ra[n][2]=rside;
					ra[n][3]=rlong;
					pa[n][0]=p[0];
					pa[n][1]=p[1];
					pa[n][2]=p[2];
					pa[n][3]=p[3];
					n+=1;
					success=true;
				}
			}
		}
	}
	return success;
}

bool CSourceCalc_OSCAR_MultiBin::IDMatch(int ident,int *idlist,int nid){
	int i=0;
	bool answer=false;
	while(answer==false && i<nid){
		if(ident==idlist[i]) answer=true;
		i+=1;
	}
	return answer;
}
