#ifndef __CWAVEFUNCTION_WF_H__
#define __CWAVEFUNCTION_WF_H__

#include "msu_commonutils/commondefs.h"
#include "msu_commonutils/parametermap.h"
#include "msu_commonutils/constants.h"
#include "msu_commonutils/sf.h"
#include "gsl/gsl_sf.h"
#include "msu_commonutils/log.h"

using namespace std;
namespace NMSUPratt{

	// These are defined outside the class because they are sometimes used by other programs
	namespace WaveFunctionRoutines{
		void getphaseshift_pipi(int I,int ell,double q,double *delta,double *ddeltadq);
		void getphaseshift_kpi(int twoI,int ell,double q,double *delta,double *ddeltadq);
		// These are used for the Reid potential
		void creid93(const double *r, const char *pname, const char *type, double *v11, double *v12, double *v22);
		double r93_w(const double *n, const double *m, const double *x);
		double r93_y(const double *n, const double *m, const double *x);
		double r93_p(const double *n, const double *m, const double *x);
		double r93_z(const double *n, const double *m, const double *x);
	};

	class CPlaneWave;
	class CPartWave;

	class CWaveFunction{
	public:
		CparameterMap parameters;
		int GetNQMAX();
		int GetNCHANNELS();
		double GetDELTA(int ichannel,int iq);
		double GetQ(int iq);
		double GetDELQ();
		bool GetIDENTICAL();
		void PrintCdelta(double Rx,double Ry,double Rz);
		double GetPsiSquared(double *pa,double *xa,double *pb,double *xb);
		double GetPsiSquared(double q,double r,double ctheta);
		virtual double CalcPsiSquared(int iq,double r,double ctheta);
		CWaveFunction();
		virtual ~CWaveFunction();
		void PrintPhaseShifts();
		double **Wepsilon,**delta,**ddeltadq,*eta,*channelweight;
		double GetIW(int ell,double epsilon,double q,int q1q2,double eta,double delta);
		double m1,m2;
		CPlaneWave **planewave;
		CPartWave ***partwave;
		void getqrctheta(double *pa,double *xa,double *pb,double *xb,double *q,double *r,double *ctheta);
		void getqrctheta(FourVector &pa,FourVector &xa,FourVector &pb,FourVector &xb,double &q,double &r,double &ctheta);
		char message[300];
	public:
		void ParsInit(string parsfilename);
		complex<double> ci;
		double MPI,MKAON,MPROTON,MLAMBDA,MNEUTRON;
		int *ell;
		//double **Wepsilon,**delta,**ddeltadq,*eta,*channelweight;
		int nchannels,q1q2;
		bool generic;
		double mu,muscale,symmweight,epsilon;
		int q1q2scale;
		bool STRONG,COULOMB,IDENTICAL;
		//CPlaneWave **planewave;
		//CPartWave ***partwave;
		// These last two are redefined (overwritten) in inherited classes
		void InitArrays();
		void InitWaves();
		void EffectiveRange(int ichannel,double scattlength,double Reff);
		//double GetIW(int ell,double epsilon,double q,int q1q2,
		//     double eta,double delta);
		void phaseshift_CoulombCorrect(int ell,double q,double eta,
		double &delta,double &ddeltadq);
		int ellmax,nqmax;
		double delq;
		double *qarray;
		// These are used for square-well solutions
		int *nwells;
		complex<double> ***A;
		complex<double> **cgsqwell;
		complex<double> ***DelPhiArray;
		int DelPhiArray_NRMAX;
		double DelPhiArray_DELR;
		void SquareWell_CalcDelPhi(int iq,double r,complex<double> *DelPhi);
		double **a;
		double **V0;
		void SquareWell_Init();
		void SquareWell_GetDelPhi(int iq,double r,complex<double> *DelPhi);
		void SquareWell_MakeArrays();
		void SquareWell_DeleteArrays();

		// This is done to account for the fact that there is a momentum dependence in the Coulomb potential
		double RelativisticCorrection(double r,int iq);
		static Crandy *randy;
	};

	class CPlaneWave{
	public:
		CPlaneWave(double etaset,int Q1Q2,double qset);
		complex<double> planewave(double r,double ctheta);
		complex<double> hyper(complex<double> a,complex<double> b,complex<double> cz); // This is Kummer's function M(a,b,z)
	public:
		complex<double>ci;
		double delx;
		int q1q2,nxmax;
		complex<double> couly;
		complex<double> cfact1;
		complex<double> cfact2;
		complex<double> chype[61];
		complex<double> chype1[11];
		complex<double> chype2[11];
		double q,eta; // q is the reduced mom., (p1-p2)/2
		complex<double> *hyperarray;
		char message[300];
	};

	class CPartWave{
	public:
		CPartWave(double etaset,int q1q2,double qset,int ell,double epsilonset);
		CPartWave(double etaset,int q1q2set,double qset,int ellset,
		double epsilonset,int nrmaxset,double delrset);
		~CPartWave();
		complex<double> GetPhiIncoming(double r);
		complex<double> GetPhiOutgoing(double r);
		double sigma;
	public:
		complex<double> ci;
		double epsilon;
		int ell,q1q2;
		double q,eta;
		complex<double> *phi; /* An incoming Coulomb partial wave
			which behaves as e^{-i(kr-eta*ln(2kr)+sigma)} */ 
			int nrmax;
		double delr;
		void phi_init();
		char message[300];
	};

	// Parent class is CWaveFunction

	class CWaveFunction_pipluspiplus_nostrong : public CWaveFunction {
	public:
		double CalcPsiSquared(int iq,double r,double ctheta);
		CWaveFunction_pipluspiplus_nostrong(string parsfilename);
	};

	class CWaveFunction_pipluspiminus_nostrong : public CWaveFunction {
	public:
		double CalcPsiSquared(int iq,double r,double ctheta);
		CWaveFunction_pipluspiminus_nostrong(string parsfilename);
	};

	class CWaveFunction_pkplus_phaseshift : public CWaveFunction {
	public:
		double CalcPsiSquared(int iq,double r,double ctheta);
		void read_phaseshifts();
		CWaveFunction_pkplus_phaseshift(string parsfilename);
	};
	// Parent class is CWaveFunction

	class CWaveFunction_ppiplus_phaseshift : public CWaveFunction {
	public:
		double CalcPsiSquared(int iq,double r,double ctheta);
		void read_phaseshifts();
		CWaveFunction_ppiplus_phaseshift(string parsfilename);
	};

	class CWaveFunction_generic : public CWaveFunction {
	public:
		double CalcPsiSquared(int iq,double r,double ctheta);
		void reset(int q1q2,double m1,double m2,double symmweight);
		CWaveFunction_generic(string parsfilename,int q1q2,double m1,double m2,double symmweight);
		bool KILLSYM;
	};

	class CWaveFunction_Xipi_phaseshift : public CWaveFunction {
	public:
		double CalcPsiSquared(int iq,double r,double ctheta);
		void get_phaseshifts_Xistar(double q,double &delta,double &ddeltadq);
		CWaveFunction_Xipi_phaseshift(string parsfilename);
	};

	class CWaveFunction_pn_phaseshift : public CWaveFunction {
	public:
		double CalcPsiSquared(int iq,double r,double ctheta);
		//void get_phaseshifts_pn(double q,double &delta,double &ddeltadq);
		void read_phaseshifts();
		CWaveFunction_pn_phaseshift(string parsfilename);
	};

	class CWaveFunction_plambda_phaseshift : public CWaveFunction {
	public:
		double CalcPsiSquared(int iq,double r,double ctheta);
		//void get_phaseshifts_plambda(double q,double &delta,double &ddeltadq);
		void get_phaseshifts();
		CWaveFunction_plambda_phaseshift(string parsfilename);
	};

	class CWaveFunction_kpluspiminus_phaseshift : public CWaveFunction {
	public:
		double CalcPsiSquared(int iq,double r,double ctheta);
		void get_phaseshifts_kpluspiminus();
		CWaveFunction_kpluspiminus_phaseshift(string parsfilename);
	};

	class CWaveFunction_pp_phaseshift : public CWaveFunction {
	public:
		double CalcPsiSquared(int iq,double r,double ctheta);
		//void get_phaseshifts_pp(double q,double &delta,double &ddeltadq);
		void read_phaseshifts();
		CWaveFunction_pp_phaseshift(string parsfilename);
	};

	class CWaveFunction_nn_phaseshift : public CWaveFunction {
	public:
		double CalcPsiSquared(int iq,double r,double ctheta);
		//void get_phaseshifts_nn(double q,double &delta,double &ddeltadq);
		void read_phaseshifts();
		CWaveFunction_nn_phaseshift(string parsfilename);
	};

	class CWaveFunction_pipluspiminus_phaseshift : public CWaveFunction {
	public:
		double CalcPsiSquared(int iq,double r,double ctheta);
		CWaveFunction_pipluspiminus_phaseshift(string parsfilename);
	};

	class CWaveFunction_pipluspiplus_phaseshift : public CWaveFunction {
	public:
		double CalcPsiSquared(int iq,double r,double ctheta);
		CWaveFunction_pipluspiplus_phaseshift(string parsfilename);
	};

	class CWaveFunction_lambdalambda_phaseshift : public CWaveFunction {
	public:
		double CalcPsiSquared(int iq,double r,double ctheta);
		CWaveFunction_lambdalambda_phaseshift(string parsfilename);
		void GetPhaseshifts();
	};

	class CWaveFunction_lambdalambdaantiparspin_phaseshift : public CWaveFunction {
	public:
		double CalcPsiSquared(int iq,double r,double ctheta);
		CWaveFunction_lambdalambdaantiparspin_phaseshift(string parsfilename);
		void GetPhaseshifts();
	};

	class CWaveFunction_lambdalambdaparspin_phaseshift : public CWaveFunction {
	public:
		double CalcPsiSquared(int iq,double r,double ctheta);
		CWaveFunction_lambdalambdaparspin_phaseshift(string parsfilename);
		void GetPhaseshifts();
	};

	class CWaveFunction_pipluspiminus_sqwell : public CWaveFunction{
	public:
		double CalcPsiSquared(int iq,double r,double ctheta);
		CWaveFunction_pipluspiminus_sqwell(string parsfilename);

		~CWaveFunction_pipluspiminus_sqwell();

	};

	class CWaveFunction_kpluspiminus_sqwell : public CWaveFunction{
	public:
		double CalcPsiSquared(int iq,double r,double ctheta);
		CWaveFunction_kpluspiminus_sqwell(string parsfilename);

		~CWaveFunction_kpluspiminus_sqwell();

	};

	class CWaveFunction_pipluspiplus_sqwell : public CWaveFunction{
	public:
		double CalcPsiSquared(int iq,double r,double ctheta);
		CWaveFunction_pipluspiplus_sqwell(string parsfilename);

		~CWaveFunction_pipluspiplus_sqwell();

	};

	class CWaveFunction_kpluspiplus_sqwell : public CWaveFunction{
	public:
		double CalcPsiSquared(int iq,double r,double ctheta);
		CWaveFunction_kpluspiplus_sqwell(string parsfilename);
		~CWaveFunction_kpluspiplus_sqwell();
	};

	class CWaveFunction_pkplus_sqwell : public CWaveFunction{
	public:
		double CalcPsiSquared(int iq,double r,double ctheta);
		CWaveFunction_pkplus_sqwell(string parsfilename);
		~CWaveFunction_pkplus_sqwell();
	};

	class CWaveFunction_ppiplus_sqwell : public CWaveFunction{
	public:
		double CalcPsiSquared(int iq,double r,double ctheta);
		CWaveFunction_ppiplus_sqwell(string parsfilename);
		~CWaveFunction_ppiplus_sqwell();

	};

	class CWaveFunction_ppiminus_sqwell : public CWaveFunction{
	public:
		double CalcPsiSquared(int iq,double r,double ctheta);
		CWaveFunction_ppiminus_sqwell(string parsfilename);
		~CWaveFunction_ppiminus_sqwell();
	};

	class	CWaveFunction_pp_schrod : public CWaveFunction{
	public:
		double CalcPsiSquared(int iq,double r,double ctheta);
		CWaveFunction_pp_schrod(string parsfilename);

		~CWaveFunction_pp_schrod();

		int nrmax_schrod;
		double rmax_schrod;
		complex<double> ***delpsi;

	public:

		void schrodinger(int iq,int ichannel);
		double Vreid(double r,int ichannel);
	};

	class	CWaveFunction_nn_schrod : public CWaveFunction{
	public:
		double CalcPsiSquared(int iq,double r,double ctheta);
		CWaveFunction_nn_schrod(string parsfilename);

		~CWaveFunction_nn_schrod();

		int nrmax_schrod;
		double rmax_schrod;
		complex<double> ***delpsi;

	public:

		void schrodinger(int iq,int ichannel);
		double Vreid(double r,int ichannel);
	};

	class	CWaveFunction_classical{
	public:
		double delq;
		double CalcPsiSquared(int iq,double r,double ctheta,double m1,double m2,int q1q2);
		double CalcPsiSquared(double q,double r,double ctheta,double m1,double m2,int q1q2);
		CWaveFunction_classical();
		~CWaveFunction_classical(); 
		char message[200];
	};

	class CWaveFunction_ppbar_nocoulomb : public CWaveFunction{
	public:
		double CalcPsiSquared(int iq,double r,double ctheta);
		CWaveFunction_ppbar_nocoulomb(string parsfilename);
		~CWaveFunction_ppbar_nocoulomb();
	public:
		vector<int> llmax;
		double VR,VI,a;
		vector<vector<complex<double>>> A,B;
		vector<complex<double>> qinside;
		void GetBessel(complex<double> x,int iq,vector<complex<double>> &jl,vector<complex<double>> &jlprime);
		void GetHankel(double x,int iq,vector<complex<double>> &hl,vector<complex<double>> &hlprime);
	};


	class CWaveFunction_optical : public CWaveFunction{
	public:
		double CalcPsiSquared(int iq,double r,double ctheta);
		CWaveFunction_optical(string parsfilename);
		~CWaveFunction_optical();
		void Reset(int q1q2,double m1,double m2,double VR,double VI,double Rset);
		void ClearInfo();
		void GetCLsigmaL();
		void GetExpansionCoefficients();
		void GetF_Complex(int iq,int L,complex<double> rho,complex<double> &F,complex<double> &Fprime);
		void GetAB();
	
		vector<int> llmax;
		int kmax;
		double VR,VI,R;
		vector<vector<vector<complex<double>>>> a;
		vector<vector<complex<double>>> A,B,CL,sigmaL;
		vector<complex<double>> qinside,etavec;
		vector<double> sigma_annihilation;
	};

	class CWaveFunction_pd_sqwell : public CWaveFunction{
	public:
		double CalcPsiSquared(int iq,double r,double ctheta);
		CWaveFunction_pd_sqwell(string parsfilename);
		~CWaveFunction_pd_sqwell();
	};
	
	class CWaveFunction_pd_tune : public CWaveFunction{
	public:
		double CalcPsiSquared(int iq,double r,double ctheta);
		CWaveFunction_pd_tune(string parsfilename);
		~CWaveFunction_pd_tune();
	};

	class CWaveFunction_pHe3_sqwell : public CWaveFunction{
	public:
		double CalcPsiSquared(int iq,double r,double ctheta);
		CWaveFunction_pHe3_sqwell(string parsfilename);
		~CWaveFunction_pHe3_sqwell();
		void CountNodes();
	};

	class CWaveFunction_Kminusd_sqwell : public CWaveFunction{
	public:
		double CalcPsiSquared(int iq,double r,double ctheta);
		CWaveFunction_Kminusd_sqwell(string parsfilename);
		~CWaveFunction_Kminusd_sqwell();
	};

	class CWaveFunction_pid_sqwell : public CWaveFunction{
	public:
		double CalcPsiSquared(int iq,double r,double ctheta);
		CWaveFunction_pid_sqwell(string parsfilename);
		~CWaveFunction_pid_sqwell();
	};

}

#endif

