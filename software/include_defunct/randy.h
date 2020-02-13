#ifndef __RANDY_H__
#define __RANDY_H__

#include <cstdlib>
#include <cmath>
#include <array>
#include <random>
#include "coralutils.h"

using namespace std;

class CRandy{
public:
	CRandy(int iseed);
	int seed;
	double ran();
	double ran_gauss();
	void ran_gauss2(double &g1,double &g2);
	double ran_exp();
	void reset(int iseed);
	void generate_boltzmann(double mass,double T,FourVector &p);
	void generate_boltzmann_alt(double mass,double T,FourVector &p);
	int poisson();
	void set_mean(double mu);  // For Poisson Dist
	
private:
	std::mt19937 mt;
	std::uniform_real_distribution<double> ranu;
	std::normal_distribution<double> rang;
	std::poisson_distribution<int> ranp;
};

#endif