#ifndef __INCLUDE_SFIT_MINUIT_GX1D_CC__
#define __INCLUDE_SFIT_MINUIT_GX1D_CC__

CCF2S_Minuit_GX1D::CCF2S_Minuit_GX1D(CSourceCalc *scset,CCHArray *cexpset,
				 CCHArray *cerrorset,CCHArray *ctheoryset,
				 CCHArray *sourceset,CKernel *kernelset){
  ndim=1;

  // npars is different for different subclasses
  npars=5;
  //
  sourcecalc=scset;
  cexp=cexpset;
  cerror=cerrorset;
  ctheory=ctheoryset;
  source=sourceset;
  kernel=kernelset;
  if(pars!=NULL) delete [] pars;
  pars=new CMNPars[npars];
  if(xval!=NULL) delete [] xval;
  xval=new double[npars];


  // initialization of pars is also unique to given subclass
  pars[0].Set("lambda",(sourcecalc->spars).getD("lambdaG",0.6),1.0,0,0.0);
  pars[1].Set("Xfrac",(sourcecalc->spars).getD("Xfrac",0.5),1.0,0,0);
  pars[2].Set("R",(sourcecalc->spars).getD("R",5),1.0,0.0,20.0);
  pars[3].Set("X",(sourcecalc->spars).getD("X",10),1.0,0.0,20.0);
  pars[4].Set("a",(sourcecalc->spars).getD("a",5),1.0,0,0);

  InitMinuit();
  FixPar(4);
}
#endif
