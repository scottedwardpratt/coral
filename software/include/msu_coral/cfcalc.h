#ifndef __INCLUDE_CFCALC_H__
#define __INCLUDE_CFCALC_H__

#include "msu_commonutils/commondefs.h"

using namespace std;

namespace CFCalc{
  double GetChiSquared(C3DArray *CFexp,C3DArray *Error,C3DArray *CFtheory);
  double GetChiSquared(int lx,int ly,int lz,CCHArray *CFexp,CCHArray *Error,CCHArray *CFtheory);
	double GetChiSquared(CCHArray *CFexp,CCHArray *CFtheory,CCHArray *error);
};

#endif
