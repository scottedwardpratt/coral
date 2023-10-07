#ifndef __INCLUDE_KERNEL_H
#define __INCLUDE_KERNEL_H

#include "msu_commonutils/commondefs.h"
#include "msu_commonutils/sf.h"
#include "msu_commonutils/constants.h"

using namespace std;

namespace NMSUPratt{

	class CKernel{
	public:
		CKernel( string kparsfilename="" );
		virtual ~CKernel();
		virtual double GetValue( int ell, double q, double r );
		double GetValue( int ell, int iq, int ir );
		bool Read(CparameterMap& parameters );
		bool Write(CparameterMap& parameters);
		void ReadData(string datadirname );
		void WriteData(string datadirname );
		void Print();
		int GetLMAX();
		double GetDELR();
		double GetDELQ();
		int GetNQMAX();
		int GetNRMAX();
		double GetQ(int iq);
		virtual void Calc( CWaveFunction *wf );
		void Calc_ClassCoul( double ma, double mb, int zazb );
		void Calc_PureHBT();
		bool GetIDENTICAL();
		double GetPsiSquared( int iq, int ir, double ctheta );
		double GetPsiSquared( int iq, double r, double ctheta );
		double GetPsiSquared( double q, double r, double ctheta);
		string getKernelFilename( string datadir, int ell, double q );
		char message[300];
	
	private:
		bool IDENTICAL;
		int ellmax;
		int nrmax,nqmax;
		double delr,delq;
		double ***kernel;
		double *P;
		void ParsInit( string kparsfilename );
		double CalcPsiSquared( int iq, int ir, double ctheta );
		double CalcPsiSquared( int iq, double r, double ctheta );
		void CalcP( double ctheta );
	};

	class CKernelExactHBT: public CKernel {
	public:
  
		CKernelExactHBT( string kparsfilename="" ): CKernel( kparsfilename ){}
		double GetValue( int ell, double q, double r ){
			return pow(-1.0,ell)*Bessel::jn(ell,2.0*q*r/HBARC);
		}
		char message[300];
	};


	class CKernelWF{
	public:

		double GetPsiSquared(int iq,int ir,int ictheta);
		double GetPsiSquared(int iq,int ir,double ctheta);
		double GetPsiSquared(int iq,double r,double ctheta);
		double GetPsiSquared(double q,double r,double ctheta);
		void Calc(CWaveFunction *wf);
		void ReadData( string datadirname );
		void WriteData( string datadirname );
		double GetDELR();
		double GetDELQ();
		double GetDELCTHETA();
		int GetNQMAX();
		int GetNRMAX();
		int GetNCTHETA();
		bool GetIDENTICAL();
		CKernelWF( string kparsfilename );
		~CKernelWF();
		char message[300];

	private:

		bool IDENTICAL;
		int nctheta,nrmax,nqmax;
		double delr,delq,delctheta;
		double ***wfarray;
		void ParsInit( string kparsfilename );
  
	};

}

#endif
