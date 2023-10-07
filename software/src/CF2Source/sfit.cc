#include  <ctime>
#include "msu_commonutils/arrays.h"
#include "msu_coral/sourcecalc.h"
#include "msu_coral/sfit.h"
#include "msu_coral/cfcalc.h"
#include "msu_coral/minimization.h"
#include "msu_commonutils/misc.h"
#include "msu_commonutils/randy.h"
#include "msu_commonutils/parametermap.h"
#include "msu_coral/source2cf.h"
#include "msu_commonutils/gslmatrix.h"

using namespace std;
using namespace NMSUPratt;

int  CCF2SFit::nmaxpars=10;
bool  CCF2SFit::MCsourceflag=0;

void  CParInfo::Set(string nameset, double  xset, double  errorset, double  xminset, double  xmaxset){
	strcpy(name,nameset.c_str());
	currentx=xset;
	error=errorset;
	xmin=xminset;
	xmax=xmaxset;
}

CParInfo::CParInfo(){
	name= new   char [20];
}

CParInfo::~CParInfo(){
	delete  [] name;
}

void  CCF2SFit::SetCalcFlag( int  calcflagset){
	calcflag=calcflagset;
}

void  CCF2SFit::SetMCSourceFlag( bool  MCsourceflagset){
	MCsourceflag=MCsourceflagset;
}

void  CCF2SFit::SetPar(string parstring, double  xset){
	int  i;
	char  parname[20];
	strcpy(parname,parstring.c_str());
	i=0;
	while (Misc::comparestrings(par[i]->name,parname)==0 && i<nmaxpars){
		i+=1;
	}
	if (i<nmaxpars) par[i]->currentx=xset;
	else{
		snprintf(message,CLog::CHARLENGTH,"Can not set %s, parameter with that name does not exist\n",parname);
		CLog::Info(message);
	}
}

double CCF2SFit::GetPar(string parstring){
	int  i;
	double xvalue;
	char  parname[20];
	strcpy(parname,parstring.c_str());
	i=0;
	while (Misc::comparestrings(par[i]->name,parname)==0 && i<nmaxpars){
		i+=1;
	}
	if(i<nmaxpars) xvalue=par[i]->bestx;
	else{
		snprintf(message,CLog::CHARLENGTH,"Can not set %s, parameter with that name does not exist\n",parname);
		CLog::Info(message);
		xvalue=0.0;
	}
	return xvalue;
}

void  CCF2SFit::SetPar(string parstring, double  xset, double  errorset,double  xminset, double  xmaxset){
	int  i;
	char  parname[20];
	strcpy(parname,parstring.c_str());
	i=0;
	while (Misc::comparestrings(par[i]->name,parname)==0 && i<nmaxpars){
		i+=1;
	}
	if (i<nmaxpars) par[i]->Set(parstring,xset,errorset,xminset,xmaxset);
	else {
		snprintf(message,CLog::CHARLENGTH,"Can not set %s, parameter with that name does not exist\n",
		parname);
		CLog::Info(message); 
	}
}

void  CCF2SFit::AddPar(string parstring, double  xset, double  errorset,double  xminset, double  xmaxset){
	if (npars>nfreepars) SwitchPars(nfreepars,npars);
	nfreepars+=1;
	npars+=1;
	if (nfreepars<nmaxpars){
		par[nfreepars-1]->Set(parstring,xset,errorset,xminset,xmaxset);
		(sourcecalc->spars).set(parstring,xset);
	}
	else{
		snprintf(message,CLog::CHARLENGTH,"Too Many Parameters! Increase static int CCF2SFit::nmaxpars\n");
		CLog::Fatal(message);
	}
	ErrorMatrix[nfreepars-1][nfreepars-1]=errorset*errorset;
	StepMatrix[nfreepars-1][nfreepars-1]=errorset;
}

void  CCF2SFit::UseBestPars(){
	int  i;
	currentchisquared=bestchisquared;
	for (i=0;i<nfreepars;i++) par[i]->currentx=par[i]->bestx;
}

void  CCF2SFit::SetL( int  lxset, int  lyset, int  lzset){
	lx=lxset;
	ly=lyset;
	lz=lzset;
}

void  CCF2SFit::SwitchPars( int  i, int  j){
	int  k;
	CParInfo *partemp;
	partemp=par[j];
	par[j]=par[i];
	par[i]=partemp;

	for (k=0;k<nmaxpars;k++){
		SwitchValues(&ErrorMatrix[i][k],&ErrorMatrix[j][k]);
		SwitchValues(&ErrorMatrix[k][i],&ErrorMatrix[k][j]);
		SwitchValues(&StepMatrix[i][k],&StepMatrix[j][k]);
		SwitchValues(&StepMatrix[k][i],&StepMatrix[k][j]);
	}
}

void  CCF2SFit::SwitchValues( double  *a, double  *b){
	double  dummy;
	dummy=*a;
	*a=*b;
	*b=dummy;
}

void  CCF2SFit::FreePar(string parstring){
	int  i;
	char  parname[20];
	strcpy(parname,parstring.c_str());
	i=nfreepars;
	while (Misc::comparestrings(par[i]->name,parname)==0 && i<npars){
		i+=1;
	}
	if (i==nmaxpars){
		snprintf(message,CLog::CHARLENGTH,"Par %s already free or does not exist\n",
		parstring.c_str());
		CLog::Info(message);
	}
	else{
		if (i!=nfreepars) SwitchPars(i,nfreepars);
		par[nfreepars]->fixed=false;
		snprintf(message,CLog::CHARLENGTH,"Freeing par[%d], name=%s\n",nfreepars,par[nfreepars]->name);
		CLog::Info(message);
		nfreepars+=1;
	}
}

void  CCF2SFit::FixPar(string parstring){
	int  i;
	char  parname[20];
	strcpy(parname,parstring.c_str());
	i=0;
	while (Misc::comparestrings(par[i]->name,parname)==0 && i<nfreepars){
		i+=1;
	}
	if (i==nfreepars){
		snprintf(message,CLog::CHARLENGTH,"Par %s already fixed or does not exist\n",
		parstring.c_str());
		CLog::Info(message);
	}
	else{
		if (i!=nfreepars-1) SwitchPars(i,nfreepars-1);
		nfreepars-=1;
		par[nfreepars]->fixed=true;
		snprintf(message,CLog::CHARLENGTH,"Fixing par[%d], name=%s\n",nfreepars,par[nfreepars]->name);
		CLog::Info(message);
	}
}

void  CCF2SFit::PrintPars(){
	int  i;
	snprintf(message,CLog::CHARLENGTH,"ipar        name     value       error       min         max     fixed\n");
	CLog::Info(message);
	for (i=0;i<npars;i++){
		snprintf(message,CLog::CHARLENGTH,"%2d : %12s %11.4e %11.4e %11.4e %11.4e  %d\n",
		i,par[i]->name,par[i]->currentx,par[i]->error,
		par[i]->xmin,par[i]->xmax,int(par[i]->fixed));
		CLog::Info(message);
	}
}

void  CCF2SFit::Init(){
	int  i,j;
	ResetChiSquared();
	nfreepars=npars=0;
	ncalls=0;
	randy= new  Crandy(-1234);
	par= new  CParInfo *[nmaxpars];
	ErrorMatrix= new   double  *[nmaxpars];
	StepMatrix= new   double  *[nmaxpars];
	for (i=0;i<nmaxpars;i++){
		par[i]= new  CParInfo();
		par[i]->bestx=0.0;
		par[i]->xmin=-1.0E10; par[i]->xmax=1.0E10;
		ErrorMatrix[i]= new   double [nmaxpars];
		StepMatrix[i]= new   double [nmaxpars];
		par[i]->fixed=false;
		for (j=0;j<nmaxpars;j++){
			ErrorMatrix[i][j]=0.0;
			StepMatrix[i][j]=0.0;
		}
	}
	snprintf(message,CLog::CHARLENGTH,"SFit Initialized\n");
	CLog::Info(message);
}

void CCF2SFit::ResetChiSquared(){
	currentchisquared=1.0E20;
	bestchisquared=1.0E20;
}

void  CCF2SFit::InitErrorMatrix(){
	int  i,j;
	for (i=0;i<nfreepars;i++){
		for (j=0;j<nfreepars;j++){
			if (i!=j){
				ErrorMatrix[i][j]=StepMatrix[i][j]=0.0;
			}
			else{
				StepMatrix[i][j]=par[i]->error;
				ErrorMatrix[i][j]=par[i]->error*par[i]->error;
			}
		}
	}
}

void  CCF2SFit::UpdateStepMatrix(){
	int  i,j;
	double  *SMeigenval;
	CGSLMatrix_Real *matrixcalc;
	matrixcalc= new  CGSLMatrix_Real(nfreepars);
	SMeigenval= new   double [nfreepars];
	matrixcalc->EigenFind(ErrorMatrix,StepMatrix,SMeigenval);
	for (i=0;i<nfreepars;i++){
		if (SMeigenval[i]<-1.0E-10){
			snprintf(message,CLog::CHARLENGTH,"FATAL: In UpdateStepSize, negative eigenvalue, =%g\n",
				SMeigenval[i]);
			CLog::Fatal(message);
		}
	}
	for (i=0;i<nfreepars;i++){
		par[i]->error=0.0;
		for (j=0;j<nfreepars;j++){
			StepMatrix[i][j]=StepMatrix[i][j]
				*sqrt(fabs(SMeigenval[j]));
			par[i]->error+=StepMatrix[i][j]*StepMatrix[i][j];
		}
		par[i]->error=sqrt(par[i]->error);
	}
	delete  [] SMeigenval;
	delete (matrixcalc);
}

void  CCF2SFit::PrintErrorMatrix(){
	int  ia,ib;
	snprintf(message,CLog::CHARLENGTH,"________ < delX_i delX_j > ________\n");
	CLog::Info(message);
	for (ia=0;ia<nfreepars;ia++){
		for (ib=0;ib<nfreepars;ib++){
			snprintf(message,CLog::CHARLENGTH,"%9.2e ",ErrorMatrix[ia][ib]);
			CLog::Info(message);
		}
		snprintf(message,CLog::CHARLENGTH,"\n");
		CLog::Info(message);
	}
	snprintf(message,CLog::CHARLENGTH,"___________________________________\n");
	CLog::Info(message);
}

void  CCF2SFit::PrintStepMatrix(){
	int  ia,ib;
	snprintf(message,CLog::CHARLENGTH,"________ < StepMatrix_ij > ________\n");
	CLog::Info(message);
	for (ia=0;ia<nfreepars;ia++){
		for (ib=0;ib<nfreepars;ib++){
			snprintf(message,CLog::CHARLENGTH,"%9.2e ",StepMatrix[ia][ib]);
			CLog::Info(message);
		}
		snprintf(message,CLog::CHARLENGTH,"\n");
		CLog::Info(message);
	}
	snprintf(message,CLog::CHARLENGTH,"___________________________________\n");
	CLog::Info(message);
}

CCF2SFit::CCF2SFit(){
	sourceCH=NULL;
	source3D=NULL;
	lista=NULL;
	listb=NULL;
	kernel=NULL;
	kernelwf=NULL;
	wf=NULL;
	cexp3D=NULL;
	cerror3D=NULL;
	ctheory3D=NULL;
	cexpCH=NULL;
	cerrorCH=NULL;
	ctheoryCH=NULL;
	Init();
}

CCF2SFit::CCF2SFit(CCHArray *sourceCHset,C3DArray *source3Dset,CMCList *listaset,CMCList *listbset,CKernel *kernelset,CKernelWF *kernelwfset,CWaveFunction *wfset,C3DArray *cexp3Dset,C3DArray *cerror3Dset,C3DArray *ctheory3Dset,CCHArray *cexpCHset,CCHArray *cerrorCHset,CCHArray *ctheoryCHset){
	sourceCH=sourceCHset;
	source3D=source3Dset;
	lista=listaset;
	listb=listbset;
	kernel=kernelset;
	kernelwf=kernelwfset;
	wf=wfset;
	cexp3D=cexp3Dset;
	cerror3D=cerror3Dset;
	ctheory3D=ctheory3Dset;
	cexpCH=cexpCHset;
	cerrorCH=cerrorCHset;
	ctheoryCH=ctheoryCHset;
	Init();
}

CCF2SFit::~CCF2SFit(){
	int  i;
	delete  [] par;
	for (i=0;i<nmaxpars;i++){
		delete  [] ErrorMatrix[i];
		delete  [] StepMatrix[i];
	}
	delete  [] ErrorMatrix;
	delete  [] StepMatrix;
}

double  CCF2SFit::GetChiSquared( double  *xx){
	double  chisquared=1.0E40;
	int  i;
	double  *x;
	ncalls+=1;
	x= new   double [nfreepars];
	for (i=0;i<nfreepars;i++) x[i]=xx[i];

	for (i=0;i<nfreepars;i++){
		(sourcecalc->spars).set(par[i]->name,x[i]);
	}
	if (calcflag==1){
		sourcecalc->CalcS(lx,ly,lz,sourceCH);
		S2CF::s2c(sourceCH,kernel,ctheoryCH);
		chisquared=CFCalc::GetChiSquared(lx,ly,lz,cexpCH,cerrorCH,ctheoryCH);
	}
	else   if (calcflag==2){
		sourcecalc->CalcS(sourceCH);
		S2CF::s2c(sourceCH,kernel,ctheoryCH);
		ctheoryCH->FillRemainderX();
		ArrayCalc::Calc3DArrayFromAExpArray(ctheoryCH,ctheory3D);
		chisquared=CFCalc::GetChiSquared(cexp3D,cerror3D,ctheory3D);

	}
	else   if (calcflag==3){
		sourcecalc->CalcS(lista,listb);
		S2CF::s2c(lista,listb,kernelwf,ctheory3D);
		chisquared=CFCalc::GetChiSquared(cexp3D,cerror3D,ctheory3D);
	}
	else   if (calcflag==4){
		sourcecalc->CalcS(lista,listb);
		S2CF::s2c(lista,listb,wf,ctheory3D);
		chisquared=CFCalc::GetChiSquared(cexp3D,cerror3D,ctheory3D);
	}

	else   if (calcflag==5){
		time_t start,end;
		time(&start);

		sourcecalc->CalcS(sourceCH);
		S2CF::s2c(sourceCH,kernel,ctheoryCH);
		ctheoryCH->FillRemainderX();
		ArrayCalc::Calc3DArrayFromAExpArray(ctheoryCH,ctheory3D);
		chisquared=CFCalc::GetChiSquared(cexp3D,cerror3D,ctheory3D);

		time(&end);
		double  dif = difftime(end, start);
		snprintf(message,CLog::CHARLENGTH,"%g seconds ",dif);
		CLog::Info(message);

	}
	else   if (calcflag==6){
		sourcecalc->CalcS(source3D);
		S2CF::s2c(source3D,kernelwf,ctheory3D);
		chisquared=CFCalc::GetChiSquared(cexp3D,cerror3D,ctheory3D);
	}
	else if(calcflag==7){
		sourcecalc->CalcS(sourceCH);
		S2CF::s2c(sourceCH,kernel,ctheoryCH);
		ctheoryCH->FillRemainderX();
		chisquared=CFCalc::GetChiSquared(cexpCH,cerrorCH,ctheoryCH);
	}
	else if(calcflag==8){
		sourcecalc->GaussCFCalc(ctheory3D);
		chisquared=CFCalc::GetChiSquared(cexp3D,cerror3D,ctheory3D);	
	}
	else if(calcflag==9){
		sourcecalc->randy->reset(-1234);
		S2CF::s2c_gauss(sourcecalc,kernelwf,ctheory3D);
		chisquared=CFCalc::GetChiSquared(cexp3D,cerror3D,ctheory3D);
	}
	else if(calcflag==10){
		S2CF::s2c_bowlersinyukov(sourcecalc,kernel,ctheory3D);
		chisquared=CFCalc::GetChiSquared(cexp3D,cerror3D,ctheory3D);
	}
	if (chisquared<bestchisquared){
		bestchisquared=chisquared;
		for (i=0;i<nfreepars;i++)
			par[i]->bestx=x[i];
	}
	if(chisquared!=chisquared){
		snprintf(message,CLog::CHARLENGTH,"GetChiSquared Failed: chisquared=%g\n",chisquared);
		CLog::Fatal(message);
	}
	return  chisquared;
	delete  [] x;
}

void  CCF2SFit::Newton(int  maxtries){
	bool  success,screwy;
	int  i,j,k,itry,nrescale,nsuccess;
	double  **curvature,*slope,**eigenvec;
	double  *x,*xnew,*xa,*xb,*xc,*xd,*delx,*dx,*eigenval;
	double  chi2,chi2a,chi2b,chi2c,chi2d,newchi2,scheck;
	CGSLMatrix_Real *matrixcalc;
	matrixcalc= new  CGSLMatrix_Real(nfreepars);
	slope= new   double [nfreepars];
	x= new   double [nfreepars];
	xnew= new   double [nfreepars];
	xa= new   double [nfreepars];
	xb= new   double [nfreepars];
	xc= new   double [nfreepars];
	xd= new   double [nfreepars];
	delx= new   double [nfreepars];
	dx= new   double [nfreepars];
	eigenval=new double[nfreepars];
	curvature= new   double *[nfreepars];
	eigenvec=new double *[nfreepars];
	for (i=0;i<nfreepars;i++){
		curvature[i]= new   double [nfreepars];
		eigenvec[i]=new double[nfreepars];
	}

	for (i=0;i<nfreepars;i++){
		dx[i]=par[i]->error;
		x[i]=par[i]->currentx;
	}

	nsuccess=itry=0;
	do{
		itry+=1;
		chi2=GetChiSquared(x);
		NEWSTART:
		for (i=1;i<nfreepars;i++){
			for (j=0;j<i;j++){
				for (k=0;k<nfreepars;k++) xa[k]=xb[k]=xc[k]=xd[k]=x[k];
				xa[i]+=dx[i]; xa[j]+=dx[j];
				xb[i]+=dx[i]; xb[j]-=dx[j];
				xc[i]-=dx[i]; xc[j]+=dx[j];
				xd[i]-=dx[i]; xd[j]-=dx[j];
				chi2a=GetChiSquared(xa);
				chi2b=GetChiSquared(xb);
				chi2c=GetChiSquared(xc);
				chi2d=GetChiSquared(xd);
				curvature[i][j]=(chi2a+chi2d-chi2b-chi2c)/(4.0*dx[i]*dx[j]);
				curvature[j][i]=curvature[i][j];
				snprintf(message,CLog::CHARLENGTH,"___ itry=%d, Curvature[%d][%d]=%g ___\n",
				itry,i,j,curvature[i][j]);
				CLog::Info(message);
			}
		}
		for (i=0;i<nfreepars;i++){
			for (k=0;k<nfreepars;k++) xa[k]=xb[k]=x[k];
			xa[i]+=dx[i];
			xb[i]-=dx[i];
			chi2a=GetChiSquared(xa);
			chi2b=GetChiSquared(xb);
			curvature[i][i]=(chi2a-2*chi2+chi2b)/(dx[i]*dx[i]);
			snprintf(message,CLog::CHARLENGTH,"___ itry=%d, Curvature[%d][%d]=%g ___\n",itry,i,i,curvature[i][i]);
			CLog::Info(message);
			slope[i]=(chi2a-chi2b)/(2.0*dx[i]);
		}
		matrixcalc->SolveLinearEqs(slope,curvature,delx);
		matrixcalc->EigenFind(curvature,eigenvec,eigenval);
		screwy=0;
		for(i=0;i<nfreepars;i++){
			if(eigenval[i]<0.0){
				snprintf(message,CLog::CHARLENGTH,"NEWTON's METHOD: screwy eigenval[i]=%g, is less than zero\n",eigenval[i]);
				CLog::Info(message);
				screwy=1;
			}
		}
		if(screwy==1){
			if(fabs(chi2-bestchisquared)<1.0E-20) exit(1);
			for(i=0;i<nfreepars;i++) x[i]=par[i]->bestx;
			chi2=bestchisquared;
			goto NEWSTART;
		}
		nrescale=0;
		TRY_RESCALED:
		snprintf(message,CLog::CHARLENGTH,"delx[] = ");
		CLog::Info(message);
		for (i=0;i<nfreepars;i++){
			snprintf(message,CLog::CHARLENGTH,"%g ",delx[i]);
			CLog::Info(message);
			xnew[i]=x[i]-delx[i];
		}
		snprintf(message,CLog::CHARLENGTH,"\n");
		CLog::Info(message);

		for (i=0;i<nfreepars;i++){
			if (par[i]->fixed==false){
				if (xnew[i]<par[i]->xmin || xnew[i]>par[i]->xmax){
					for (j=0;j<nfreepars;j++) delx[j]=0.5*delx[j];
					snprintf(message,CLog::CHARLENGTH,"Stepped outside min/max, will rescale delx[]\n");
					CLog::Info(message);
					if (nrescale>4) exit(1);
					nrescale+=1;
					goto  TRY_RESCALED;
				}
			}
		}

		newchi2=GetChiSquared(xnew);
		if (newchi2>chi2){
			for (j=0;j<nfreepars;j++){
				delx[j]=0.5*delx[j];
			}
			snprintf(message,CLog::CHARLENGTH,"new chi^2 bigger than previous, will rescale delx[]\n");
			CLog::Info(message);
			goto  TRY_RESCALED;
		}

		success=1;
		scheck=0.0;
		for (i=0;i<nfreepars;i++){
			for (j=0;j<nfreepars;j++){
				scheck+=curvature[i][j]*delx[i]*delx[j];
			}
		}
		if (scheck<1.0) success=1;

		if (success){
			snprintf(message,CLog::CHARLENGTH,"SUCCESS!!!!!!!!!!!\n");
			CLog::Info(message);
			nsuccess+=1;
			CalcErrorMatrixFromCurvature(curvature);
			for (i=0;i<nfreepars;i++) dx[i]=par[i]->error;
		}

		chi2=newchi2;
		for (i=0;i<nfreepars;i++) x[i]=xnew[i];
	}while (itry<maxtries && nsuccess<2);

	currentchisquared=chi2;
	for (i=0;i<nfreepars;i++) par[i]->currentx=x[i];

	delete  [] dx;
	delete  [] x;
	delete  [] xnew;
	delete  [] delx;
	delete  [] xa;
	delete  [] xb;
	delete  [] xc;
	delete  [] xd;
	delete  [] slope;
	delete [] eigenval;
	for (i=0;i<nfreepars;i++){
		delete [] curvature[i];
		delete [] eigenvec[i];
	}
	delete  [] curvature;
	delete [] eigenvec;
	delete (matrixcalc);
}

void  CCF2SFit::CalcErrorMatrixFromCurvature( double  **C){
	int  i,j,k;
	CGSLMatrix_Real *matrixcalc;
	double  **U;
	double  *EigenVal;
	matrixcalc= new  CGSLMatrix_Real(nfreepars);
	EigenVal= new   double [nfreepars];
	U= new   double  *[nfreepars];
	for (i=0;i<nfreepars;i++) U[i]= new   double [nfreepars];
	matrixcalc->EigenFind(C,U,EigenVal);

	for (i=0;i<nfreepars;i++){
		for (j=0;j<=i;j++){
			if (i<nfreepars && j<nfreepars){
				ErrorMatrix[i][j]=0.0;
				for (k=0;k<nfreepars;k++){
					ErrorMatrix[i][j]+=U[i][k]*U[j][k]/EigenVal[k];
				}
			}
			else  ErrorMatrix[i][j]=0.0;
			if (i!=j) ErrorMatrix[j][i]=ErrorMatrix[i][j];
		}
		par[i]->error=sqrt(fabs(ErrorMatrix[i][i]));
	}

	for (i=0;i<nfreepars;i++)  delete  [] U[i];
	delete  [] U;
	delete  [] EigenVal;
	delete (matrixcalc);
}

// Li changes (the argument maxcall is NOT used.)
void  CCF2SFit::ConjugateGradient(int  maxcalls){
	int i;
	double  * x =  new   double [nfreepars];
	for (i=0;i<nfreepars;i++)
	{
		x[i] = par[i]->currentx;
	}
	SetDimension(nfreepars);
	double  fmin;
	int  iter;
	if (conjugate_gradient(x, iter, fmin) ==  true )
	{
		snprintf(message,CLog::CHARLENGTH,"SUCCESS in Conjugate gradient method after %d iterations\n",iter);
		CLog::Info(message);
		currentchisquared=fmin;
		for (i=0;i<nfreepars;i++)
		{
			par[i]->currentx = x[i];
		}
	}
	if(iter>maxcalls){
		snprintf(message,CLog::CHARLENGTH,"In CCF2SFit::ConjugateGradient, iter=%d, maxcalls=%d\n",iter,maxcalls);
		CLog::Info(message);
	}

	snprintf(message,CLog::CHARLENGTH,"Best chi^2=%g, Best x[] = ",bestchisquared);
	CLog::Info(message);
	for (i=0;i<nfreepars;i++){
		snprintf(message,CLog::CHARLENGTH,"%g ",par[i]->bestx);
		CLog::Info(message);
	}
	snprintf(message,CLog::CHARLENGTH,"\n");
	CLog::Info(message);

	delete []  x;
}

// Li changes
double  CCF2SFit::fn( double  * x){
	return  GetChiSquared(x);
}

// Li changes
bool  CCF2SFit::dfn( double  * x){
	double  * x_minus =  new   double [nfreepars];
	double  * x_plus =  new   double [nfreepars];
	double  delta_x;

	for ( int  i=0;i<nfreepars;i++)
	{
		x_minus[i] = x[i];
		x_plus[i] = x[i];
	}
	for ( int  i=0;i<nfreepars;i++)
	{
	//delta_x = fabs(x[i] * 1.0e-6) > 1.0e-8 ? fabs(x[i] * 1.0e-6) : 1.0e-8;
		delta_x = fabs(x[i] * 1.0e-3) > 1.0e-4 ? fabs(x[i] * 1.0e-3) : 1.0e-4;
		x_minus[i] = x[i] - delta_x;
		x_plus[i] = x[i] + delta_x;
		vec_dx[i] = 0.5 * (GetChiSquared(x_plus) - GetChiSquared(x_minus)) / delta_x;
		x_minus[i] = x[i];
		x_plus[i] = x[i];
	}
	for ( int  i=0;i<nfreepars;i++)
	{
		snprintf(message,CLog::CHARLENGTH,"dx[%d] = %g, ", i, vec_dx[i]);
		CLog::Info(message);
	}
	snprintf(message,CLog::CHARLENGTH,"\n");
	CLog::Info(message);
	delete []  x_minus;
	delete []  x_plus;
	return  true ;
}

void  CCF2SFit::Metropolis(int  maxcalls){
	int  icall,i,j,Nsuccess=0;
	bool  success;
	double  stepscale,step,chisquared;
	double  *x,*xran,*xbar;
	x= new   double [nfreepars];
	xran= new   double [nfreepars];
	xbar= new   double [nfreepars];
	for (i=0;i<nfreepars;i++){
		x[i]=par[i]->currentx;
		xbar[i]=0.0;
		for (j=0;j<nfreepars;j++) ErrorMatrix[i][j]=0.0;
	}
	stepscale=1.0/sqrt( double (nfreepars));
	for (icall=0;icall<maxcalls;icall++){

		GETNEWXRAN:
		if (MCsourceflag){
			for (i=0;i<nfreepars;i++) x[i]=par[i]->currentx;
			currentchisquared=GetChiSquared(x);
		}
		for (j=0;j<nfreepars;j++) 
			xran[j]=randy->ran_gauss();
		for (i=0;i<nfreepars;i++){
			if (par[i]->fixed==false){
				x[i]=par[i]->currentx;
				for (j=0;j<nfreepars;j++){
					step=StepMatrix[i][j]*xran[j]*stepscale;
					x[i]+=step;
				}
			}
		}
		for (i=0;i<nfreepars;i++)
			if (x[i]<par[i]->xmin || x[i]>par[i]->xmax)  goto  GETNEWXRAN;

		chisquared=GetChiSquared(x);
		success=0;
		if (chisquared<currentchisquared){
			success=1;
		}
		else   if (randy->ran()<exp(-0.5*(chisquared-currentchisquared))){
			success=1;
		}

		if (success==1){
			Nsuccess+=1;
			currentchisquared=chisquared;
			for (i=0;i<nfreepars;i++){
				par[i]->currentx=x[i];
			}
			snprintf(message,CLog::CHARLENGTH,"SUCCESS, Nsuccess=%d\n",Nsuccess);
			CLog::Info(message);
		}

		for (i=0;i<nfreepars;i++){
			xbar[i]+=x[i];
			for (j=0;j<nfreepars;j++) ErrorMatrix[i][j]+=x[i]*x[j];
		}
	}

	for (i=0;i<nfreepars;i++){
		par[i]->xbar=xbar[i]/ double (maxcalls);
		for (j=0;j<nfreepars;j++) ErrorMatrix[i][j]=ErrorMatrix[i][j]/ double (maxcalls);
	}
	for (i=0;i<nfreepars;i++)
		for (j=0;j<nfreepars;j++) ErrorMatrix[i][j]=(ErrorMatrix[i][j]-par[i]->xbar*par[j]->xbar)
		* double (Nsuccess)/ double (Nsuccess-1);

	snprintf(message,CLog::CHARLENGTH,"Best chi^2=%g, Best x[] = ",bestchisquared);
	CLog::Info(message);
	for (i=0;i<nfreepars;i++){
		snprintf(message,CLog::CHARLENGTH,"%g ",par[i]->bestx);
		CLog::Info(message);
	}
	snprintf(message,CLog::CHARLENGTH,"\n");
	CLog::Info(message);

	delete  [] x;
	delete  [] xran;
	delete  [] xbar;
}

void  CCF2SFit::SteepestDescent( int  maxtries){
	int  i,j,itry,nfailure,nsuccess;
	double  chisquared,newchisquared,chi2a,chi2b,qstep;
	double  d2chi2dq2,dchi2dq,dq,qhatnorm;
	double  *delx,*delq,*x,*xa,*xb,*qxratio,*xnew,*qhat;
	delx= new   double [nfreepars];
	x= new   double [nfreepars];
	xa= new   double [nfreepars];
	xb= new   double [nfreepars];
	xnew= new   double [nfreepars];
	delq= new   double [nfreepars];
	qxratio= new   double [nfreepars];
	qhat=new double[nfreepars];
	for (i=0;i<nfreepars;i++){
		qxratio[i]=par[i]->error;
		x[i]=par[i]->currentx;
	}
	chisquared=GetChiSquared(x);

	for(itry=0;itry<maxtries;itry++){
		// First find qhat by investigating Grad_q chi2
		dq=0.01;
		qhatnorm=0.0;
		for (i=0;i<nfreepars;i++){
			delx[i]=qxratio[i]*dq;
			for (j=0;j<nfreepars;j++){
				xa[j]=x[j];
				xb[j]=x[j];
			}
			xb[i]=x[i]+0.5*delx[i];
			xa[i]=x[i]-0.5*delx[i];
			chi2a=GetChiSquared(xa);
			chi2b=GetChiSquared(xb);
			qhat[i]=-(chi2b-chi2a)*qxratio[i]/delx[i];
			qhatnorm+=qhat[i]*qhat[i];
		}
		for(i=0;i<nfreepars;i++) qhat[i]=qhat[i]/sqrt(qhatnorm);

		// Now do Newton's method along qhat
		nsuccess=0;
		while(nsuccess<3){
			for (i=0;i<nfreepars;i++){
				delq[i]=qhat[i]*dq;
				delx[i]=delq[i]*qxratio[i];
				xa[i]=x[i]-0.5*delx[i];
				xb[i]=x[i]+0.5*delx[i];
			}
			chi2a=GetChiSquared(xa);
			chi2b=GetChiSquared(xb);
			dchi2dq=(chi2b-chi2a)/dq;
			d2chi2dq2=4.0*(chi2b+chi2a-2.0*chisquared)/(dq*dq);
			qstep=-dchi2dq/d2chi2dq2;
			snprintf(message,CLog::CHARLENGTH,"||||||||| qstep=%g, ",qstep);
			CLog::Info(message);
			if(d2chi2dq2<0.0){
				snprintf(message,CLog::CHARLENGTH,"upside down curvature\n");
				CLog::Info(message);
			}
			if (fabs(qstep)>1.0) qstep=qstep/fabs(qstep);
			if(qstep*dchi2dq>0) qstep=-0.5*qstep/fabs(qstep);
			snprintf(message,CLog::CHARLENGTH,"-> %g |||||||||\n",qstep);
			CLog::Info(message);

			TRYNEWQSTEP:
			nfailure=0;
			for (i=0;i<nfreepars;i++){
				delq[i]=qstep*qhat[i];
				xnew[i]=x[i]+delq[i]*qxratio[i];
			}
			newchisquared=GetChiSquared(xnew);
			if(newchisquared>chisquared){
				qstep=0.5*qstep;
				nfailure+=1;
				if (nfailure>5){
					snprintf(message,CLog::CHARLENGTH,"STEEPEST DESCENT FAILURE\n");
					CLog::Info(message);
				}
				goto TRYNEWQSTEP;
			}
			snprintf(message,CLog::CHARLENGTH,"++++++++++++ qstep=%g ++++++++++++++\n",qstep);
			CLog::Info(message);
			chisquared=newchisquared;
			for(i=0;i<nfreepars;i++) x[i]=xnew[i];
			if(fabs(qstep)<2.0*dq){
				dq=dq/3.0;
				nsuccess+=1;
			}

		}
		snprintf(message,CLog::CHARLENGTH,"_______________________ finished itry=%d _______________________________\n",itry);
		CLog::Info(message);
	}

	chisquared=GetChiSquared(x);
	currentchisquared=chisquared;
	for (i=0;i<nfreepars;i++){
		if (par[i]->fixed==false) par[i]->currentx=x[i];
	}

	delete [] xa;
	delete [] xb;
	delete [] x;
	delete [] xnew;
	delete [] delx;
	delete [] delq;
	delete [] qxratio;
	delete [] qhat;
}
