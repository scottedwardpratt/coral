#include <sstream>
#include <iomanip>

#include "msu_coral/kernel.h"
#include "msu_commonutils/sf.h"
#include "msu_commonutils/parametermap.h"
#include "msu_commonutils/misc.h"
#include "msu_coral/wavefunction.h"

using namespace std;
using namespace NMSUPratt;

string CKernel::getKernelFilename( string datadir, int ell, double q )
{
	stringstream filename;
	filename << datadir 
		<< "/ell" << ell 
			<< "_q" << setfill('0') << fixed << setw(7) << setprecision(2) << q << ".tmp";
	return filename.str();
}

/**
*  \brief Constructor for CKernel
*  \param kparsfilename name of the file that has the parameters to set up the kernel
*/
CKernel::CKernel(string kparsfilename){
	int iq,ir,ell,dell;
	ParsInit(kparsfilename);
	// dell is the difference in angular kernelum
	dell=1;
	if(IDENTICAL) dell=2;

	kernel=new double **[ellmax+1];
	P=new double [ellmax+1];
	for(ell=0;ell<=ellmax;ell+=1){
		kernel[ell]=NULL;
		if(ell%dell==0){
			kernel[ell]=new double *[nqmax];
			for(iq=0;iq<nqmax;iq++){
				kernel[ell][iq]=new double[nrmax];
				for(ir=0;ir<nrmax;ir++) kernel[ell][iq][ir]=0.0;
			}
		}
	}
}

/**
*  \brief Destructor for CKernel
*/
CKernel::~CKernel(){
	int ell,iq,dell=1;
	if(IDENTICAL) dell=2;
	delete [] P;
	for(ell=0;ell<=ellmax;ell+=dell){
		for(iq=0;iq<nqmax;iq++) delete [] kernel[ell][iq];
		if(kernel[ell]!=NULL) delete [] kernel[ell];
	}
	delete [] kernel;
	CLog::Info("kernel object deleted\n");

}

//! Getter method for ellmax
//! \return CKernel::ellmax
int CKernel::GetLMAX(){
	return ellmax;
}

//! Getter method for delr
//! \return CKernel::delr
double CKernel::GetDELR(){
	return delr;
}
double CKernel::GetQ(int iq){
	return (iq+0.5)*delq;
}

//! Getter method for delq
//! \return CKernel::delq
double CKernel::GetDELQ(){
	return delq;
}

//! Getter method for nrmax
//! \return CKernel::nrmax
int CKernel::GetNRMAX(){
	return nrmax;
}

//! Getter method for nqmax
//! \return CKernel::nqmax
int CKernel::GetNQMAX(){
	return nqmax;
}

//! Getter method for IDENTICAL
//! \return CKernel::IDENTICAL
bool CKernel::GetIDENTICAL(){
	return IDENTICAL;
}

/**
*  Performs the linear interpolation to get the value of the kernel at a particular point
*  \param  ell The Legendre order of the requested kernel
*  \param  q   The value of relative momentum at which you want to evaluate the kernel (in MeV/c)
*  \param  r   The value of separation at which you want to evaluate the kernel (in fm)
*  \return Value of the kernel at the requested point
*/
double CKernel::GetValue(int ell,double q,double r){
	int iqlow,iqhigh,irlow,irhigh,iq,ir;
	double wrlow,wrhigh,wqlow,wqhigh,answer;
	double qlow,qhigh,rlow,rhigh,qi,ri;

	// Perform simple linear extrapolation
	if((IDENTICAL && pow(-1.0,ell)<0) || q>delq*nqmax || r>delr*nrmax) answer=0.0;
	else{
		iq=lrint(floor(q/delq));
		qi=(0.5+iq)*delq;
		if(q<qi){
			iqlow=iq-1;
			if(iqlow<0) iqlow=0;
			iqhigh=iq;
			qhigh=(0.5+iq)*delq;
			wqhigh=fabs(delq-fabs(qhigh-q))/delq;
			if(wqhigh>1.0) wqhigh=1.0;
			wqlow=1.0-wqhigh;
		}
		else{
			iqlow=iq;
			iqhigh=iq+1;
			if(iqhigh>=nqmax) iqhigh=iqlow;
			qlow=(iq+0.5)*delq;
			wqlow=fabs(delq-fabs(q-qlow))/delq;
			if(wqlow>1.0) wqlow=1.0;
			wqhigh=1.0-wqlow;
		}
	
		ir=int(floor(r/delr));
		ri=(0.5+ir)*delr;
		if(r<ri){
			irlow=ir-1;
			if(irlow<0) irlow=0;
			irhigh=ir;
			rhigh=(0.5+ir)*delr;
			wrhigh=fabs(delr-fabs(rhigh-r))/delr;
			//wrhigh=1.0;
			if(wrhigh>1.0) wrhigh=1.0;
			wrlow=1.0-wrhigh;
		}
		else{
			irlow=ir;
			irhigh=ir+1;
			if(irhigh>=nrmax) irhigh=irlow;
			rlow=(ir+0.5)*delr;
			wrlow=fabs(delr-fabs(r-rlow))/delr;
			if(wrlow>1.0) wrlow=1.0;
			wrhigh=1.0-wrlow;
		}
	
		if(fabs(wrlow+wrhigh-1.0)>1.0E-8 || fabs(wqlow+wqhigh-1.0)>1.0E-6){
			CLog::Fatal("weights do not add up in GetValue(ell,q,r) in CKernel::GetValue\n");
		}
	
		answer=wqlow*wrlow*GetValue(ell,iqlow,irlow)
			+wqlow*wrhigh*GetValue(ell,iqlow,irhigh)
				+wqhigh*wrlow*GetValue(ell,iqhigh,irlow)
					+wqhigh*wrhigh*GetValue(ell,iqhigh,irhigh);
	
		if(ell==0 && r>1.0 && answer<-1.0){
			CLog::Info("__________________________________________________\n");
			snprintf(message,CLog::CHARLENGTH,"answer below -1, =%g, q=%g, r=%g, iq=%d, ir=%d, nrmax=%d\n",
			answer,q,r,iq,ir,nrmax);
			CLog::Info(message);
			snprintf(message,CLog::CHARLENGTH,"wqlow=%g, wqhigh=%g, wrlow=%g, wrhigh=%g\n",wqlow,wqhigh,wrlow,wrhigh);
			CLog::Info(message);
			CLog::Info("__________________________________________________\n");
		}
	}
	return answer;
}

/**
* \brief  The value of a particular point in the interpolation table
* \param  ell The Legendre order of the table in question
* \param  iq  The q grid index
* \param  ir  The r grid index
* \return The value of the interpolation table at the specified index on the grid
*/
double CKernel::GetValue(int ell,int iq,int ir){
	if(iq>=nqmax||ir>=nrmax || iq<0 || ir<0) return 0.0;
	else{
		return kernel[ell][iq][ir];
	}
}

/**
* \brief Read in the parameters from a CparameterMap
* \param parameters The CparameterMap containing the parameters (see below)
* \return True==Success, False==Failure
*
* The parameters in the CparameterMap that this code uses are as follows:
*   - \c nqmax     The number of points in the q direction in the kernel interpolation table
*   - \c nrmax     The number of points in the r direction in the kernel interpolation table
*   - \c kellmax   The maximum order of the Legendre polynomial expansion of the kernel to compute
*   - \c delq      The step size in the q direction in the kernel interpolation table
*   - \c delr      The step size in the r direction in the kernel interpolation table
*   - \c IDENTICAL Bool to denote whether this kernel is for identical pairs
*/
bool CKernel::Read(CparameterMap& parameters){
	nqmax=parameters.getI("NQMAX",200);
	nrmax=parameters.getI("NRMAX",500);
	ellmax=parameters.getI("KLMAX",4);
	delq=parameters.getD("DELQ",0.5);
	delr=parameters.getD("DELR",0.1);
	IDENTICAL=parameters.getB("IDENTICAL",true);
	return true;
}


/**
* \brief Writes the parameters to a CparameterMap
* \param parameters The CparameterMap containing the parameters (see below)
* \return True==Success, False==Failure
*
* The parameters in the CparameterMap that this code uses are as follows:
*   - \c nqmax     The number of points in the q direction in the kernel interpolation table
*   - \c nrmax     The number of points in the r direction in the kernel interpolation table
*   - \c kellmax   The maximum order of the Legendre polynomial expansion of the kernel to compute
*   - \c delq      The step size in the q direction in the kernel interpolation table
*   - \c delr      The step size in the r direction in the kernel interpolation table
*   - \c IDENTICAL Bool to denote whether this kernel is for identical pairs
*/
bool CKernel::Write( CparameterMap& parameters){
	parameters.set("NQMAX",nqmax);
	parameters.set("NRMAX",nrmax);
	parameters.set("KLMAX",ellmax);
	parameters.set("DELQ",delq);
	parameters.set("DELR",delr);
	parameters.set("IDENTICAL",IDENTICAL);
	return true;
}

void CKernel::ParsInit(string kparsfilename){
	CparameterMap parameters;
	if ( kparsfilename!=string("") ) {
		parameters.ReadParsFromFile(kparsfilename);
	}
	Read(parameters);

	CLog::Info("Parameters for kernel: \n");
	snprintf(message,CLog::CHARLENGTH,"delq=%g,nqmax=%d, delr=%g, nrmax=%d,KLMAX=%d\n",delq,nqmax,delr,nrmax,ellmax);
	CLog::Info(message);
	snprintf(message,CLog::CHARLENGTH,"IDENTICAL set to %d\n",int(IDENTICAL));
	CLog::Info(message);
	CLog::Info("__________________________________________\n");
}

/**
* \brief Read in the interpolation table from a directory
* \param datadir The directory containing the interpolation table.  
*
* The data in datadir is stored in one file per q.  Each file contains 
* an interpolation table in r for that value of q.
*/
void CKernel::ReadData(string datadir){
	double q,delr0;//,r;
	int iq,ir,ell,nrmax0,dell;
	FILE *infile;
	dell=1;
	string filename;
	if(IDENTICAL) dell=2;
	// Check if all the needed files are there
	for(ell=0;ell<=ellmax;ell+=dell){
		bool this_l_OK = true;
		for(iq=0;iq<nqmax;iq++){
			q=(0.5+iq)*delq;
			filename = getKernelFilename( datadir, ell, q );
			if (! Misc::file_exists( filename ) ) this_l_OK = false;
		}
		if ( !this_l_OK )
		{
			cout << "Missing files for l = "<< ell<<", truncating kernel at lmax = "<<ell-dell<< endl;
			ellmax = ell-dell;
			break;
		}
	} 
	// Load the files
	for(ell=0;ell<=ellmax;ell+=dell){
		for(iq=0;iq<nqmax;iq++){
			q=(0.5+iq)*delq;
			filename = getKernelFilename( datadir, ell, q );
			infile=fopen(filename.c_str(),"r");
			if(infile==NULL){
				cout <<"Opening kernel data file "<<filename<<" failed!"<< endl;
				exit(1);
			}
			fscanf(infile,"%d %lf",&nrmax0,&delr0);
			if(nrmax0<nrmax || fabs(delr-delr0)>1.0E-8){
				CLog::Fatal("Inconsistent values for nrmax or delr in data files! in CKernel::ReadData\n");
			}
			for(ir=0;ir<nrmax;ir++){
				fscanf(infile,"%lf",&kernel[ell][iq][ir]);
			}
			fclose(infile);
		}
	}
}

/** 
* \brief Writes the interpolation table to bunch of files in a directory
* \param datadir The directory containing the interpolation table.  
*
* The data in datadir is stored in one file per q.  Each file contains 
* an interpolation table in r for that value of q.
*/
void CKernel::WriteData(string datadir){
	double q;
	int ell,iq,ir,dell;
	char filename[100];
	FILE *outfile;
	char shellcommand[120];
	snprintf(shellcommand,120,"mkdir -p %s",datadir.c_str());

	system(shellcommand);

	dell=1;
	if(IDENTICAL) dell=2;
	for(ell=0;ell<=ellmax;ell+=dell){
		for(iq=0;iq<nqmax;iq++){
			q=(0.5+iq)*delq;
			snprintf(filename,100,"%s/ell%d_q%07.2f.tmp",datadir.c_str(),ell,q);
			snprintf(message,CLog::CHARLENGTH,"For q=%g, Will write to %s\n",q,filename);
			CLog::Info(message);
			outfile=fopen(filename,"w");
			fprintf(outfile,"%d %g\n",nrmax,delr);
			for(ir=0;ir<nrmax;ir++){
				//r=(0.5+ir)*delr;
				fprintf(outfile,"%17.10e\n",kernel[ell][iq][ir]);
			}
			fclose(outfile);
		}
	}
}


/** 
* \brief Simple routine to print out the entire interpolation table to stdout for debugging
*/
void CKernel::Print(){
	double q,r;
	int ell,iq,ir,dell;//,nrmax0;

	dell=1;
	if(IDENTICAL)
		dell=2;
	snprintf(message,CLog::CHARLENGTH,"ellmax=%d, nqmax=%d, nrmax=%d\n",ellmax,nqmax,nrmax);
	CLog::Info(message);
	for(iq=0;iq<nqmax;iq++){
		q=(0.5+iq)*delq;
		snprintf(message,CLog::CHARLENGTH,"______________ q=%g MeV/c____________________\n",q);
		CLog::Info(message);
		for(ir=0;ir<nrmax;ir++){
			r=(0.5+ir)*delr;
			snprintf(message,CLog::CHARLENGTH,"%6.2f ",r);
			CLog::Info(message);
			for(ell=0;ell<=ellmax;ell+=dell){
				snprintf(message,CLog::CHARLENGTH,"%10.3e ",kernel[ell][iq][ir]);
				CLog::Info(message);
			}
			CLog::Info("\n");
		}
	}
}

/**
* \brief Builds the interpolation table, assuming classical Coulomb correlations.  
* \param ma       Mass of particle a, in MeV
* \param mb       Mass of particle b, in MeV
* \param zazb     Charge of particle a times charge or particle b.  Both in units of \f$e\f$.
*
*  THIS DOCUMENTATION NEEDS TO BE FLUSHED OUT!!!!!!
*/
void CKernel::Calc_ClassCoul(double ma,double mb,int zazb){
	double q,r,ctheta;//,delctheta;
	int ell,iq,ir,nu=512;
	double x,ea,eb,u,umax,delu;//,cweight;

	for(ell=0;ell<=ellmax;ell++){
		for(iq=0;iq<nqmax;iq++){
			for(ir=0;ir<nrmax;ir++) kernel[ell][iq][ir]=0.0;
		}
	}

	for(iq=0;iq<nqmax;iq++){
		q=(0.5+iq)*delq;
		for(ir=0;ir<nrmax;ir++){
			r=(0.5+ir)*delr;
			ea=sqrt(q*q+ma*ma);
			eb=sqrt(q*q+mb*mb);
			x=2.0*ea*eb*zazb*197.327/(q*q*r*137.036*(ea+eb));
		
			if(x<1.0){
				umax=2.0*sqrt(1.0-x);
				delu=umax/double(nu);
				for(u=0.5*delu;u<umax;u+=delu){
					ctheta=-1.0+x+sqrt(u*u+x*x);
					for(ell=0;ell<=ellmax;ell++)
						kernel[ell][iq][ir]+=0.5*delu*SpherHarmonics::legendre(ell,ctheta);
				}
				snprintf(message,CLog::CHARLENGTH,"q=%g, r=%g, kernel0=%g =? %g\n",
				q,r,kernel[0][iq][ir],sqrt(1.0-x));
				CLog::Info(message);
			}
			kernel[0][iq][ir]-=1.0;
		}
	}
}

/**
* \brief Builds the interpolation table, assuming HBT only.  
* 
*  For a pure HBT kernel, the pair relative wavefunction is simply 
*  \f$ \Psi = \frac{1}{\sqrt{2}}\left(e^{iq\cdot r}+e^{-iq\cdot r}\right) \f$
*  Where \f$q = \frac{1}{2}(q_1-q_2)\f$ is the relative four-momentum of the pair and \f$r\f$ is 
*  the space-time separation of the pair at freeze-out.  It is straight-forward to show 
*  that the kernel is
*  \f[
*      K_\ell(r,q) = (-1)^{\ell/2} j_\ell(2qr/\hbar c)
*  \f]
*  This is for even \f$\ell\f$ only.  For odd ones, this is zero.  Here \f$q=|\vec{q}|\f$ 
*  and \f$r=|\vec{r}|\f$, both in the pair rest frame.
*
*/
void CKernel::Calc_PureHBT(){
	int ir,iq,ell,sign;
	double qr,r,q;

	for(iq=0;iq<nqmax;iq++){
		q=(0.5+iq)*delq;
		for(ir=0;ir<nrmax;ir++){
			r=(0.5+ir)*delr;
			qr=q*r/197.3269602;
			sign=-1;
			for(ell=0;ell<=ellmax;ell+=2){
				sign=-sign;
				kernel[ell][iq][ir]=sign*Bessel::jn(ell,2.0*qr);
			}
			kernel[0][iq][ir]-=1.0;
		}
	}
}

/**
* \brief Builds the interpolation table, using the actual CWaveFunction in CKernel::wf
* \param wf pointer to the wavefunction to use to build the kernel
* 
* The table is built up on the regular grid defined by CKernel::nqmax, CKernel::delq,
* CKernel::nrmax, CKernel::delr, by a simple Simpson rule integration over the angle between
* \f$\vec{q}\f$ and \f$\vec{r}\f$:
* \f[
*     K_\ell(q,r) = \frac{1}{2}\int_{-1}^{1} d\mu \left[|\Psi(q,r,\mu)|^2 - 1\right] P_\ell(\mu)
* \f]
* Where \f$\mu = cos(\vec{q}\cdot\vec{r}/qr)\f$.
*/
void CKernel::Calc(CWaveFunction *wf){
	double q,r,ctheta,wf2;
	int ell,dell,iq,ir,ictheta,nctheta=720;
	if(IDENTICAL!=wf->GetIDENTICAL()){
		CLog::Info("Creating Kernel with different symmetry (value of IDENTICAL) than wave function\n");
		CLog::Fatal("You probably need to edit parameters file\n");
	}

	dell=1;
	if(IDENTICAL) dell=2;
	for(iq=0;iq<nqmax;iq++){
		q=(0.5+iq)*delq;
		for(ir=0;ir<nrmax;ir++){
			r=(0.5+ir)*delr;
			for(ell=0;ell<=ellmax;ell+=dell)
				kernel[ell][iq][ir]=0.0;
			for(ictheta=0;ictheta<nctheta;ictheta++){
				ctheta=-1.0+2.0*(0.5+ictheta)/double(nctheta);
				wf2=wf->GetPsiSquared(q,r,ctheta);
				for(ell=0;ell<=ellmax;ell+=dell){
					kernel[ell][iq][ir]+=(wf2-1.0)*SpherHarmonics::legendre(ell,ctheta);
				}
			}
			for(ell=0;ell<=ellmax;ell+=dell)
				kernel[ell][iq][ir]=kernel[ell][iq][ir]/double(nctheta);
		}
	}
}

double CKernel::CalcPsiSquared(int iq,int ir,double ctheta){
	int ell,dell=1,nsmall=0;
	double answer=1.0,dela;
	CalcP(ctheta);
	if(IDENTICAL) dell=2;
	ell=0;
	if(ir<nrmax && iq<nqmax){
		while(ell<=ellmax && nsmall<5){
			dela=kernel[ell][iq][ir]*(2.0*ell+1)*P[ell];
			//*SpherHarmonics::legendre(ell,ctheta);
			answer+=dela;
			if(fabs(dela)<1.0E-4) nsmall+=1;
			else(nsmall=0);
			ell+=dell;
		}
	}
	return answer;
}

double CKernel::GetPsiSquared(int iq,int ir,double ctheta){
	double answer=1.0;
	if(ir<nrmax && iq<nqmax){
		answer=CalcPsiSquared(iq,ir,ctheta);
	}
	return answer;
}

double CKernel::CalcPsiSquared(int iq,double r,double ctheta){
	double wa,wb;
	int ira,irb;
	ira=lrint(floor((r/delr)-0.5));
	if(ira<0) ira=0;
	irb=ira+1;
	wb=(r-(ira+0.5)*delr)/delr;
	wa=((irb+0.5)*delr-r)/delr;
	return wa*CalcPsiSquared(iq,ira,ctheta)+wb*CalcPsiSquared(iq,irb,ctheta);
}


double CKernel::GetPsiSquared(int iq,double r,double ctheta){
	double answer=1.0;
	if(iq<nqmax){
		CalcP(ctheta);
		answer=CalcPsiSquared(iq,r,ctheta);
	}
	return answer;
}



double CKernel::GetPsiSquared(double q,double r,double ctheta){
	double wa,wb,answer=1.0;
	int iqa,iqb;
	iqa=lrint(floor((q/delq)-0.5));
	if(iqa<0) iqa=0;
	iqb=iqa+1;
	if(iqb>=nqmax){
		iqb=nqmax-1;
		iqa=iqb-1;
	}
	wb=(q-(iqa+0.5)*delq)/delq;
	wa=((iqb+0.5)*delq-q)/delq;
	answer=wa*CalcPsiSquared(iqa,r,ctheta)+wb*CalcPsiSquared(iqb,r,ctheta);

	return answer;
}

void CKernel::CalcP(double ctheta){
	int ell;
	P[0]=1.0;
	P[1]=ctheta;
	for(ell=1;ell<ellmax;ell++){
		P[ell+1]=((2*ell+1)*ctheta*P[ell]-ell*P[ell-1])/double(ell+1);
	}
}    

// Routines Below for CKernelWF

CKernelWF::CKernelWF(string kparsfilename){
	int iq,ir,ictheta;
	ParsInit(kparsfilename);
	// dell is the difference in angular kernelum
	delctheta=2.0/double(nctheta);
	if(IDENTICAL) delctheta=delctheta*0.5;

	wfarray=new double **[nqmax];
	for(iq=0;iq<nqmax;iq+=1){
		wfarray[iq]=new double *[nrmax];
		for(ir=0;ir<nrmax;ir++){
			wfarray[iq][ir]=new double[nctheta+1];
			for(ictheta=0;ictheta<=nctheta;ictheta++) wfarray[iq][ir][ictheta]=0.0;
		}
	}
}

void CKernelWF::ParsInit(string kparsfilename){
	CparameterMap parameters;
	parameters.ReadParsFromFile(kparsfilename);

	nqmax=parameters.getI("NQMAX",25);
	nrmax=parameters.getI("NRMAX",120);
	nctheta=parameters.getI("NCTHETA",120);
	delq=parameters.getD("DELQ",4.0);
	delr=parameters.getD("DELR",0.5);
	IDENTICAL=parameters.getB("IDENTICAL",0);

	snprintf(message,CLog::CHARLENGTH,"reading from %s\n",kparsfilename.c_str());
	CLog::Info(message);
	snprintf(message,CLog::CHARLENGTH,"  _________ PARAMETERS FOR KERNELWF ________\n");
	CLog::Info(message);
	snprintf(message,CLog::CHARLENGTH,"delq set to %g\n",delq);
	CLog::Info(message);
	snprintf(message,CLog::CHARLENGTH,"nqmax set to %d\n",nqmax);
	CLog::Info(message);
	snprintf(message,CLog::CHARLENGTH,"delr set to %g\n",delr);
	CLog::Info(message);
	snprintf(message,CLog::CHARLENGTH,"nrmax set to %d\n",nrmax);
	CLog::Info(message);
	snprintf(message,CLog::CHARLENGTH,"nctheta set to %d\n",nctheta);
	CLog::Info(message);
	snprintf(message,CLog::CHARLENGTH,"IDENTICAL set to %d\n",IDENTICAL);
	CLog::Info(message);
	CLog::Info("__________________________________________\n");
}


CKernelWF::~CKernelWF(){
	int iq,ir;
	for(iq=0;iq<nqmax;iq++){
		for(ir=0;ir<nrmax;ir++) delete [] wfarray[iq][ir];
		delete [] wfarray[iq];
	}
	delete [] wfarray;
	CLog::Info("KernelWF object deleted\n");

}

void CKernelWF::Calc(CWaveFunction *wf){
	int iq,ir,ictheta;
	double q,r,ctheta,ps2;
	for(iq=0;iq<nqmax;iq++){
		q=(iq+0.5)*delq;
		for(ir=0;ir<nrmax;ir++){
			r=(ir+0.5)*delr;
			for(ictheta=0;ictheta<=nctheta;ictheta++){
				ctheta=1.0-ictheta*delctheta;
				if(ctheta>1.0) ctheta=1.0;
				ps2=wf->GetPsiSquared(q,r,ctheta);
				if(ps2<0.0 && r>1.0){
					snprintf(message,CLog::CHARLENGTH,"screwy, psi^2=%g, r=%g, q=%g, ctheta=%g\n",ps2,r,q,ctheta);
					CLog::Fatal(message);
				}
				wfarray[iq][ir][ictheta]=ps2-1.0;
			}
		}
	}
}

double CKernelWF::GetPsiSquared(int iq,int ir,int ictheta){
	if(iq<nqmax && ir<nrmax && ictheta<=nctheta)
		return 1.0+wfarray[iq][ir][ictheta];
	else return 1.0;  
}

double CKernelWF::GetPsiSquared(int iq,int ir,double ctheta){
	double wa,wb;
	int ictheta;
	if(iq<nqmax && ir<nrmax){
		if(IDENTICAL) ctheta=fabs(ctheta);
		ictheta=lrint(floor((1.0-ctheta)/delctheta));
		wb=((1.0-ctheta)-ictheta*delctheta)/delctheta;
		wa=1.0-wb;
		return 1.0+wa*wfarray[iq][ir][ictheta]+wb*wfarray[iq][ir][ictheta+1];
	}
	else return 1.0;
}

double CKernelWF::GetPsiSquared(int iq,double r,double ctheta){
	int ir;
	double wa,wb;
	if(iq<nqmax){
		ir=lrint(floor((r-0.5*delr)/delr));
		if(ir<0) ir=0;
		wb=(r-(ir+0.5)*delr)/delr;
		wa=1.0-wb;
		return wa*GetPsiSquared(iq,ir,ctheta)+wb*GetPsiSquared(iq,ir+1,ctheta);
	}
	else return 1.0;
}

double CKernelWF::GetPsiSquared(double q,double r,double ctheta){
	int iq;
	double wa,wb;
	iq=lrint(floor((q-0.5*delq)/delq));
	if(iq<0) iq=0;
	wb=(q-(iq+0.5)*delq)/delq;
	wa=1.0-wb;
	return wa*GetPsiSquared(iq,r,ctheta)+wb*GetPsiSquared(iq+1,r,ctheta);
}

double CKernelWF::GetDELR(){
	return delr;
}

double CKernelWF::GetDELQ(){
	return delq;
}

double CKernelWF::GetDELCTHETA(){
	return delctheta;
}

int CKernelWF::GetNRMAX(){
	return nrmax;
}

int CKernelWF::GetNQMAX(){
	return nqmax;
}

int CKernelWF::GetNCTHETA(){
	return nctheta;
}

bool CKernelWF::GetIDENTICAL(){
	return IDENTICAL;
}

void CKernelWF::WriteData(string datadir){
	double q;//,r;
	int iq,ir,ictheta;
	double meanwf2;
	char filename[100];
	FILE *outfile;
	char shellcommand[120];
	snprintf(shellcommand,120,"mkdir -p %s",datadir.c_str());
	system(shellcommand);
	for(iq=0;iq<nqmax;iq++){
		meanwf2=0.0;
		q=(0.5+iq)*delq;
		snprintf(filename,100,"%s/q%07.2f.tmp",datadir.c_str(),q);
		outfile=fopen(filename,"w");
		fprintf(outfile,"%d %d\n",nrmax,nctheta);
		for(ir=0;ir<nrmax;ir++){
			for(ictheta=0;ictheta<=nctheta;ictheta++){
				meanwf2+=wfarray[iq][ir][ictheta];
				fprintf(outfile,"%17.10e ",wfarray[iq][ir][ictheta]);
			}
			fprintf(outfile,"\n");
		}
		meanwf2=meanwf2/double(nctheta*nrmax);
		fclose(outfile);
	} 
}

void CKernelWF::ReadData(string datadir){
	double q;
	double meanwf2;
	int iq,ir,ictheta,nrmaxread,ncthetaread;
	char filename[100];
	FILE *infile;

	for(iq=0;iq<nqmax;iq++){
		meanwf2=0.0;
		q=(0.5+iq)*delq;
		snprintf(filename,100,"%s/q%07.2f.tmp",datadir.c_str(),q);
		infile=fopen(filename,"r");
		fscanf(infile,"%d %d\n",&nrmaxread,&ncthetaread);
		if(nrmaxread!=nrmax || ncthetaread!=nctheta){
			CLog::Info("CKernelWF : trying to read file with wrong dimensions\n");
			snprintf(message,CLog::CHARLENGTH,"nrmaxread=%d, nrmax=%d, ncthetaread=%d, nctheta=%d\n",
			nrmaxread,nrmax,ncthetaread,nctheta);
			CLog::Fatal(message);
		}
		for(ir=0;ir<nrmax;ir++){
			for(ictheta=0;ictheta<=nctheta;ictheta++){
				fscanf(infile,"%lf ",&wfarray[iq][ir][ictheta]);
				meanwf2+=wfarray[iq][ir][ictheta];
			}
		}
		fclose(infile);
		meanwf2=meanwf2/double(nctheta*nrmax);
	} 
}



