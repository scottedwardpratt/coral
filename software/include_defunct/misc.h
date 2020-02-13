#ifndef __INCLUDE_MISC_H_
#define __INCLUDE_MISC_H_

#include <complex>
#include <cstdlib>
#include <cstring>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_poly.h>
#include "coralutils.h"

using namespace std;

namespace Misc{
	double DotProduct(FourVector &p,FourVector &q);
	void lorentz(FourVector &u,FourVector &p,FourVector &pprime);
	void BoostToCM(FourVector &u,FourVector &p,FourVector &ptilde);
	void Boost(FourVector &u,FourVector &ptilde,FourVector &p);
	
	//void lorentz(double *u,double *p1,double *p1prime);
	// These find tensor/vector in cm frame (same as above)
	void BoostToCM(double *u,double **Pi,double **PiTilde);
	//void BoostToCM(double *u,double *p,double *ptilde);

	// Given tensor/vector in cm frame, this finds value in lab frame
	void Boost(double *u,double **PiTilde,double **Pi);
	//void Boost(double *u,double *ptilde,double *p);
	void Boost(FourVector &u,FourVector &p);

	double cgc(double j1,double m1,double j2,double m2,double j,double m);
	double cgc_edmonds(double j1,double m1,double j2,double m2,double j,double m);
	bool comparestrings(char *s1,char *s2);
	double triangle(double m0,double m1,double m2);
	double triangle2(double m0squared,double m1squared,double m2squared);
	//double PCMS(double M, double particle_i, double particle_j);
	int Sign(int a);
	void outsidelong(double *pa,double *pb,double &qinv,double &qout,double &qside,double &qlong);
	void outsidelong(double *pa,double *pb, double &qinv, double &qout, double &qside, double &qlong,double &deleta,double &dely,double &delphi);
	double GetQinv(double *pa,double *pb);
	//double GetQinv(double *pa,double *pb);
	double GetRapidity(double *p);
	double GetDely(double *pa,double *pb);

	complex<double> cexp(complex<double> z);
	complex<double> ceiphi(double phi);
	complex<double> cpow(complex<double> z,complex<double> a);

	int iround(double x);

	int cgc_delta (int x, int y);
	double cgc_factorial (double n);
	double cgc_fractorial (double n,double m);
	double oldcgc(double j1,double m1,double j2,double m2,double j,double m);

	void Cubic(double a0,double a1,double a2,double a3,
	complex<double>& z1,complex<double>& z2,complex<double>& z3);

	// Quartic routines provided by Andrew Steiner
	double signswitch(double a, double b);
	void Quartic(double a0,double a1,double a2,double a3,double a4,complex<double> &z1, complex<double> &z2,complex<double> &z3,complex<double> &z4);
	void Quartic(double a0,double a1,double a2,double a3,double a4,complex<double> *z);
	// this corresponds to z0=x[0], z1=x[1]+ix[2], z2=x[1]-ix[2]
	void CubicResolvant(double r,double s,double t,double x[],double &d);

	int CubicReal(double a0,double a1,double a2,double a3,double *x);
	void CubicComplex(double a0,double a1,double a2,double a3,complex<double> &z1,complex<double> &z2,complex<double> &z3);
	void Pause();
	void Pause(int seconds);
	double CalcDelta_FromSqWells(int ell,double mu,int nwells,double q,double *V0,double *r);
	bool file_exists( const char* file_name );
	bool file_exists( const std::string& file_name );
	double mod(double x, double y);

};

#endif
