#include "msu_coral/sfit.h"

using namespace std;

CCF2SFit_GX1D::CCF2SFit_GX1D(CSourceCalc *scset,
			     CCHArray *cexpset,
			     CCHArray *cerrorset,
			     CCHArray *ctheoryset,
			     CCHArray *sourceset,
			     CKernel *kernelset){
  int i,j;

  //
  sourcecalc=scset;
  cexpCH=cexpset;
  cerrorCH=cerrorset;
  ctheoryCH=ctheoryset;
  sourceCH=sourceset;
  kernel=kernelset;
  Init();

  // initialization of pars is also unique to given subclass

  AddPar("lambdaG",(sourcecalc->spars).getD("lambdaG",0.3),
	     0.02,0.0,1.5);
  AddPar("R",(sourcecalc->spars).getD("R",5),
	     0.2,1.0,12.0);
  AddPar("lambdaX",(sourcecalc->spars).getD("lambdaX",0.3),
	     0.02,0.0,1.5);
  AddPar("X",(sourcecalc->spars).getD("X",10.0),
	     0.4,1.0,25.0);
  AddPar("a",(sourcecalc->spars).getD("a",5.0),
	     0.2,1.0,20.0);

  for(i=0;i<nfreepars;i++){
    for(j=0;j<nfreepars;j++){
      StepMatrix[i][j]=0.0;
      ErrorMatrix[i][j]=0.0;
      if(i==j){
	StepMatrix[i][j]=par[i]->error;
	ErrorMatrix[i][j]=par[i]->error*par[i]->error;
      }
    }
  }

}
