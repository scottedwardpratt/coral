#ifndef  __INCLUDE_SFIT_H__
#define  __INCLUDE_SFIT_H__

#include  "msu_commonutils/commondefs.h"
#include "sourcecalc.h"
#include  "msu_coral/minimization.h"
#include "msu_commonutils/misc.h"
#include "msu_commonutils/randy.h"

using namespace std;
namespace NMSUPratt{

	class  CParInfo{
		public :
		bool  fixed;  // 'true' if parameter is not allowed to vary
		// searches are confined to  xmin < x < xmax
		double  xmin,xmax,error,xbar;
		double  bestx,currentx;  // bestx is the value that vave smallest chi^2
		char  *name;
		void  Set(string sname, double  xset, double  error,
		double  xminset, double  xmaxset);
		CParInfo();
		~CParInfo();
	};

	// Li changes
	class  CCF2SFit :  public  CMinimization{
		public :
		/*
		calcflag = 1 if CF is a CCHArray object of specific lx,ly,lz
		and CSourceCalc object makes CCHArray objects
		calcflag = 2 if CF is a C3DArray object 
		and CSourceCalc object makes CCHArray objects
		calcflag =3,4 are used if source functions are calculated through
		intermediate MC lists, but no such implementations yet exist
		calcflag=7, used for comparing 2 CCHArray CF functions
		*/  
		void  SetCalcFlag( int  calcflagset);
		void  SetMCSourceFlag( bool  MCsourceflagset);
	
		void  SetPar(string parstring, double  value);
		void  SetPar(string parstring, double  value, double  error,double  xmin, double  xmax);
		double GetPar(string parstring);
		void  AddPar(string parstring, double  value, double  error,double  xmin, double  xmax);
		void  PrintPars();
		void  PrintErrorMatrix();
		void  PrintStepMatrix();
		void  FixPar(string parname);
		void  FreePar(string parname);
		void  UseBestPars();
		void  SetL( int  lxset, int  lyset, int  lzset);
	
		// Li changes
		void  ConjugateGradient(int  maxcalls);
		double  fn( double  * x);
		bool  dfn( double  * x);
	
		void  Metropolis( int  maxcalls);
		void  SteepestDescent( int  maxtries);
		void  Newton( int  maxtries);
		void  UpdateStepMatrix();
		void  InitErrorMatrix();
	
		// This calculates source functions
		CSourceCalc *sourcecalc;
	
		// ______________ The objects below depend on calc_flag _____________
		// These are used to store information about the source
		CCHArray *sourceCH;
		C3DArray *source3D;
		CMCList *lista;
		CMCList *listb;
	
		// Wavefunctions or kernels
		CKernel *kernel;
		CKernelWF *kernelwf;
		CWaveFunction *wf;
	
		// These are arrays for storing correlation functions and errors
		C3DArray *cexp3D;
		C3DArray *cerror3D;
		C3DArray *ctheory3D;
		CCHArray *cexpCH;
		CCHArray *cerrorCH;
		CCHArray *ctheoryCH;
		// _____________________________________________________________________
	
	
		void ResetChiSquared();
		CCF2SFit();
		CCF2SFit(CCHArray *sourceCHset,C3DArray *source3Dset,
		CMCList *listaset,CMCList *listbset,
		CKernel *kernelset,CKernelWF *kernelwfset,
		CWaveFunction *wfset,
		C3DArray *cexp3Dset,C3DArray *cerror3Dset,
		C3DArray *ctheory3Dset,CCHArray *cexpCHset,
		CCHArray *cerrorCHset,CCHArray *ctheoryCHset);
		~CCF2SFit();
		char message[300];
	
		protected :
		int  calcflag;
		static   bool  MCsourceflag;  // If the source has a Monte Carlo nature, 
		// i.e., it fluctuates for a given parameter set, set this to true
	
		CParInfo **par;
		static   int  nmaxpars;
		int  nfreepars,npars;
		double  **ErrorMatrix;
		double  currentchisquared,bestchisquared;
	
		int  lx,ly,lz;
	
	
		Crandy *randy;
		void  SwitchPars( int  ipara, int  iparb);
		void  SwitchValues( double  *a, double  *b);
	
		int  ncalls;
		double  **StepMatrix;
		void  Init();
		double  GetChiSquared( double  *x);
		void  CalcErrorMatrixFromCurvature( double  **C);
	};

	class  CCF2SFit_Blast :  public  CCF2SFit{
		public :
		CCF2SFit_Blast(CSourceCalc *scset,C3DArray *cexpset,
		C3DArray *cerrorset,C3DArray *ctheory3Dset,
		CCHArray *ctheoryset,CCHArray *sourceset,
		CKernel *kernelset);
		char message[300];
	};

	class  CCF2SFit_GX1D :  public  CCF2SFit{
		public :
		CCF2SFit_GX1D(CSourceCalc *scset,CCHArray *cexpset,
		CCHArray *cerrorset,CCHArray *ctheoryset,
		CCHArray *sourceset,CKernel *kernelset);
		char message[300];
	};

	class  CCF2SFit_3DGaussian :  public  CCF2SFit{
		public :
		CCF2SFit_3DGaussian(CSourceCalc *scset,C3DArray *cexpset,
		C3DArray *cerrorset,C3DArray *ctheory3Dset,
		CCHArray *ctheoryset,CCHArray *sourceset,
		CKernel *kernelset);
		char message[300];
	};

}

#endif
