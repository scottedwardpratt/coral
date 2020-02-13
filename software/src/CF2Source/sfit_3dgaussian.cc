#ifndef __INCLUDE_SFIT_3DGAUSSIAN_CC__
#define __INCLUDE_SFIT_3DGAUSSIAN_CC__

#include "sfit.h"
#include "sourcecalc.h"

using namespace std;

CCF2SFit_3DGaussian::CCF2SFit_3DGaussian(CSourceCalc *scset,C3DArray *cexpset,C3DArray *cerrorset,C3DArray *ctheory3Dset,CCHArray *ctheoryset,CCHArray *sourceset,CKernel *kernelset){
	// npars is different for different subclasses
	//
	sourcecalc=scset;
	cexp3D=cexpset;
	cerror3D=cerrorset;
	ctheoryCH=ctheoryset;
	ctheory3D=ctheory3Dset;
	sourceCH=sourceset;
	kernel=kernelset;
	Init();

	// initialization of pars is also unique to given subclass

	AddPar("lambda",(sourcecalc->spars).getD("lambda",0.5),
	0.02,0.0,1.5);
	AddPar("Rx",(sourcecalc->spars).getD("Rx",5),
	0.2,1.0,12.0);
	AddPar("Ry",(sourcecalc->spars).getD("Ry",5),
	0.2,1.0,12.0);
	AddPar("Rz",(sourcecalc->spars).getD("Rz",5),
	0.2,1.0,12.0);
	AddPar("Xoff",(sourcecalc->spars).getD("Xoff",0.0),
	0.2,-10.0,10.0);

	InitErrorMatrix();

}
#endif

