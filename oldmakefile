#enter make, then make install, make clean, make uninstall

include ../../makefile_defs.mk

#################################################

INC= -I include -I ${MADAI_GSLPATH}/include  -I ${MADAI_X11PATH}/include -I ${MADAI_EIGEN_HOME} -I ../commonutils/include -I ../b3d2/include -I ../hydrobalance/include

HFILES_CORAL = include/parametermap.h include/xgraph.h include/sf.h include/misc.h include/wavefunction.h include/gslmatrix.h include/coralutils.h include/minimization.h include/constants.h include/chebyshev_fit.h include/arrays.h include/cfcalc.h include/kernel.h include/randy.h include/sourcecalc.h include/source2cf.h include/sfit_minuit.h include/sfit.h include/chebyshev.h ../b3d2/include/part.h

INSTALL_HFILES_CORAL = ${MADAI_INSTALLDIR}/include/minimization.h ${MADAI_INSTALLDIR}/include/sfit.h	${MADAI_INSTALLDIR}/include/cfcalc.h ${MADAI_INSTALLDIR}/include/kernel.h 	${MADAI_INSTALLDIR}/include/source2cf.h 	${MADAI_INSTALLDIR}/include/kernel.h ${MADAI_INSTALLDIR}/include/source2cf.h 	${MADAI_INSTALLDIR}/include/sourcecalc.h	${MADAI_INSTALLDIR}/include/wavefunction.h

install : ${MADAI_INSTALLDIR}/lib/libcoral.a 

${MADAI_INSTALLDIR}/lib/libcoral.a : commonutils lib/libcoral.a
	cp -f lib/libcoral.a ${MADAI_INSTALLDIR}/lib/

HFILES_CORAL = include/minimization.h include/sfit.h include/cfcalc.h include/kernel.h include/source2cf.h include/sourcecalc.h include/wavefunction.h

INSTALL_HFILES_CORAL = ${MADAI_INSTALLDIR}/include/minimization.h ${MADAI_INSTALLDIR}/include/sfit.h ${MADAI_INSTALLDIR}/include/cfcalc.h  ${MADAI_INSTALLDIR}/include/kernel.h ${MADAI_INSTALLDIR}/include/source2cf.h ${MADAI_INSTALLDIR}/include/sourcecalc.h ${MADAI_INSTALLDIR}/include/wavefunction.h

OBJFILES = build/minimization.o build/sfit.o build/sfit_3dgaussian.o build/sfit_blast.o build/sfit_GX1d.o build/cfcalc.o build/kernel.o build/source2cf.o build/sourcecalc.o build/sourcecalc_getgausspars.o build/sourcecalc_ellipticblast.o build/sourcecalc_gauss.o build/sourcecalc_GX1d.o build/sourcecalc_OSCAR.o build/sourcecalc_OSCAR_multibin.o build/sourcecalc_blast.o build/wf_generic.o build/wf_kpluspiminus_phaseshift.o build/wf_kpluspiminus_sqwell.o build/wf_kpluspiplus_sqwell.o build/wf_lambdalambda_phaseshift.o build/wf_lambdalambdaantiparspin_phaseshift.o build/wf_lambdalambdaparspin_phaseshift.o build/wf_nn_phaseshift.o build/wf_pipluspiminus_nostrong.o build/wf_pipluspiminus_phaseshift.o build/wf_pipluspiminus_sqwell.o build/wf_pipluspiplus_nostrong.o build/wf_pipluspiplus_phaseshift.o build/wf_pipluspiplus_sqwell.o build/wf_pkplus_phaseshift.o build/wf_pkplus_sqwell.o build/wf_plambda_phaseshift.o build/wf_pn_phaseshift.o build/wf_pp_phaseshift.o build/wf_pp_schrod.o build/wf_ppiminus_sqwell.o build/wf_ppiplus_phaseshift.o build/wf_ppiplus_sqwell.o build/wf_Xipi_phaseshift.o build/kpi_phaseshifts.o build/partwave.o build/pipi_phaseshifts.o build/planewave.o build/reid93.o build/wfcommon.o build/wfsquarewell.o

lib/libcoral.a : commonutils ${HFILES_CORAL} ${INSTALL_HFILES_CORAL} ${OBJFILES}
	mkdir -p lib;\
	rm -f lib/libcoral.a;\
	ar -ru lib/libcoral.a ${OBJFILES};

commonutils :
	make -C ../commonutils install

############ MAKE OBJECT FILES ###############

build/minimization.o : ${HFILES_CORAL} src/CF2Source/minimization.cc
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INC} src/CF2Source/minimization.cc -o build/minimization.o

build/sfit.o : 	${HFILES_CORAL} src/CF2Source/sfit.cc
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INC} src/CF2Source/sfit.cc -o build/sfit.o

build/sfit_3dgaussian.o : 	${HFILES_CORAL} src/CF2Source/sfit_3dgaussian.cc
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INC} src/CF2Source/sfit_3dgaussian.cc -o build/sfit_3dgaussian.o

build/sfit_blast.o : 	${HFILES_CORAL} src/CF2Source/sfit_blast.cc
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INC} src/CF2Source/sfit_blast.cc -o build/sfit_blast.o

build/sfit_GX1d.o : 	${HFILES_CORAL} src/CF2Source/sfit_GX1d.cc
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INC} src/CF2Source/sfit_blast.cc -o build/sfit_GX1d.o

build/cfcalc.o : ${HFILES_CORAL} src/CFCalc/cfcalc.cc
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INC} src/CFCalc/cfcalc.cc -o build/cfcalc.o

build/kernel.o : 	${HFILES_CORAL} src/Kernel/kernel.cc
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INC} src/Kernel/kernel.cc -o build/kernel.o

build/source2cf.o : ${HFILES_CORAL} src/Source2CF/source2cf.cc
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INC} src/Source2CF/source2cf.cc -o build/source2cf.o

build/sourcecalc.o : ${HFILES_CORAL} src/SourceCalc/sourcecalc.cc
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INC} src/SourceCalc/sourcecalc.cc -o build/sourcecalc.o

build/sourcecalc_getgausspars.o : ${HFILES_CORAL} src/SourceCalc/sourcecalc_getgausspars.cc
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INC} src/SourceCalc/sourcecalc_getgausspars.cc -o build/sourcecalc_getgausspars.o

build/sourcecalc_ellipticblast.o : ${HFILES_CORAL} src/SourceCalc/sourcecalc_ellipticblast.cc
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INC} src/SourceCalc/sourcecalc_ellipticblast.cc -o build/sourcecalc_ellipticblast.o

build/sourcecalc_gauss.o : ${HFILES_CORAL} src/SourceCalc/sourcecalc_gauss.cc
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INC} src/SourceCalc/sourcecalc_gauss.cc -o build/sourcecalc_gauss.o

build/sourcecalc_GX1d.o : ${HFILES_CORAL} src/SourceCalc/sourcecalc_GX1d.cc
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INC} src/SourceCalc/sourcecalc_GX1d.cc -o build/sourcecalc_GX1d.o

build/sourcecalc_OSCAR.o : ${HFILES_CORAL} src/SourceCalc/sourcecalc_OSCAR.cc
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INC} src/SourceCalc/sourcecalc_OSCAR.cc -o build/sourcecalc_OSCAR.o

build/sourcecalc_OSCAR_multibin.o : ${HFILES_CORAL} src/SourceCalc/sourcecalc_OSCAR_multibin.cc
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INC} src/SourceCalc/sourcecalc_OSCAR_multibin.cc -o build/sourcecalc_OSCAR_multibin.o

build/sourcecalc_blast.o : ${HFILES_CORAL} src/SourceCalc/sourcecalc_blast.cc
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INC} src/SourceCalc/sourcecalc_blast.cc -o build/sourcecalc_blast.o

build/wf_generic.o : ${HFILES_CORAL} src/WaveFunctions/wf/wf_generic.cc
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INC} src/WaveFunctions/wf/wf_generic.cc -o build/wf_generic.o

build/wf_kpluspiminus_phaseshift.o : ${HFILES_CORAL} src/WaveFunctions/wf/wf_kpluspiminus_phaseshift.cc
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INC} src/WaveFunctions/wf/wf_kpluspiminus_phaseshift.cc -o build/wf_kpluspiminus_phaseshift.o

build/wf_kpluspiminus_sqwell.o : ${HFILES_CORAL} src/WaveFunctions/wf/wf_kpluspiminus_sqwell.cc
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INC} src/WaveFunctions/wf/wf_kpluspiminus_sqwell.cc -o build/wf_kpluspiminus_sqwell.o

build/wf_kpluspiplus_sqwell.o : ${HFILES_CORAL} src/WaveFunctions/wf/wf_kpluspiplus_sqwell.cc
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INC} src/WaveFunctions/wf/wf_kpluspiplus_sqwell.cc -o build/wf_kpluspiplus_sqwell.o

build/wf_lambdalambda_phaseshift.o : ${HFILES_CORAL} src/WaveFunctions/wf/wf_lambdalambda_phaseshift.cc
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INC} src/WaveFunctions/wf/wf_lambdalambda_phaseshift.cc -o build/wf_lambdalambda_phaseshift.o

build/wf_lambdalambdaantiparspin_phaseshift.o : ${HFILES_CORAL} src/WaveFunctions/wf/wf_lambdalambdaantiparspin_phaseshift.cc
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INC} src/WaveFunctions/wf/wf_lambdalambdaantiparspin_phaseshift.cc -o build/wf_lambdalambdaantiparspin_phaseshift.o

build/wf_lambdalambdaparspin_phaseshift.o : ${HFILES_CORAL} src/WaveFunctions/wf/wf_lambdalambdaparspin_phaseshift.cc
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INC} src/WaveFunctions/wf/wf_lambdalambdaparspin_phaseshift.cc -o build/wf_lambdalambdaparspin_phaseshift.o

build/wf_nn_phaseshift.o : ${HFILES_CORAL} src/WaveFunctions/wf/wf_nn_phaseshift.cc
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INC} src/WaveFunctions/wf/wf_nn_phaseshift.cc -o build/wf_nn_phaseshift.o

build/wf_pipluspiminus_nostrong.o : ${HFILES_CORAL} src/WaveFunctions/wf/wf_pipluspiminus_nostrong.cc
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INC} src/WaveFunctions/wf/wf_pipluspiminus_nostrong.cc -o build/wf_pipluspiminus_nostrong.o

build/wf_pipluspiminus_phaseshift.o : ${HFILES_CORAL} src/WaveFunctions/wf/wf_pipluspiminus_phaseshift.cc
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INC} src/WaveFunctions/wf/wf_pipluspiminus_phaseshift.cc -o build/wf_pipluspiminus_phaseshift.o

build/wf_pipluspiminus_sqwell.o : ${HFILES_CORAL} src/WaveFunctions/wf/wf_pipluspiminus_sqwell.cc
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INC} src/WaveFunctions/wf/wf_pipluspiminus_sqwell.cc -o build/wf_pipluspiminus_sqwell.o

build/wf_pipluspiplus_nostrong.o : ${HFILES_CORAL} src/WaveFunctions/wf/wf_pipluspiplus_nostrong.cc
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INC} src/WaveFunctions/wf/wf_pipluspiplus_nostrong.cc -o build/wf_pipluspiplus_nostrong.o

build/wf_pipluspiplus_phaseshift.o : ${HFILES_CORAL} src/WaveFunctions/wf/wf_pipluspiplus_phaseshift.cc
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INC} src/WaveFunctions/wf/wf_pipluspiplus_phaseshift.cc -o build/wf_pipluspiplus_phaseshift.o

build/wf_pipluspiplus_sqwell.o : ${HFILES_CORAL} src/WaveFunctions/wf/wf_pipluspiplus_sqwell.cc
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INC} src/WaveFunctions/wf/wf_pipluspiplus_sqwell.cc -o build/wf_pipluspiplus_sqwell.o

build/wf_pkplus_phaseshift.o : ${HFILES_CORAL} src/WaveFunctions/wf/wf_pkplus_phaseshift.cc
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INC} src/WaveFunctions/wf/wf_pkplus_phaseshift.cc -o build/wf_pkplus_phaseshift.o

build/wf_pkplus_sqwell.o : ${HFILES_CORAL} src/WaveFunctions/wf/wf_pkplus_sqwell.cc
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INC} src/WaveFunctions/wf/wf_pkplus_sqwell.cc -o build/wf_pkplus_sqwell.o

build/wf_plambda_phaseshift.o : ${HFILES_CORAL} src/WaveFunctions/wf/wf_plambda_phaseshift.cc
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INC} src/WaveFunctions/wf/wf_plambda_phaseshift.cc -o build/wf_plambda_phaseshift.o

build/wf_pn_phaseshift.o : ${HFILES_CORAL} src/WaveFunctions/wf/wf_pn_phaseshift.cc
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INC} src/WaveFunctions/wf/wf_pn_phaseshift.cc -o build/wf_pn_phaseshift.o

build/wf_pp_phaseshift.o : ${HFILES_CORAL} src/WaveFunctions/wf/wf_pp_phaseshift.cc
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INC} src/WaveFunctions/wf/wf_pp_phaseshift.cc -o build/wf_pp_phaseshift.o

build/wf_pp_schrod.o : ${HFILES_CORAL} src/WaveFunctions/wf/wf_pp_schrod.cc
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INC} src/WaveFunctions/wf/wf_pp_schrod.cc -o build/wf_pp_schrod.o

build/wf_ppiminus_sqwell.o : ${HFILES_CORAL} src/WaveFunctions/wf/wf_ppiminus_sqwell.cc
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INC} src/WaveFunctions/wf/wf_ppiminus_sqwell.cc -o build/wf_ppiminus_sqwell.o

build/wf_ppiplus_phaseshift.o : ${HFILES_CORAL} src/WaveFunctions/wf/wf_ppiplus_phaseshift.cc
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INC} src/WaveFunctions/wf/wf_ppiplus_phaseshift.cc -o build/wf_ppiplus_phaseshift.o

build/wf_ppiplus_sqwell.o : ${HFILES_CORAL} src/WaveFunctions/wf/wf_ppiplus_sqwell.cc
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INC} src/WaveFunctions/wf/wf_ppiplus_sqwell.cc -o build/wf_ppiplus_sqwell.o

build/wf_Xipi_phaseshift.o : ${HFILES_CORAL} src/WaveFunctions/wf/wf_Xipi_phaseshift.cc
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INC} src/WaveFunctions/wf/wf_Xipi_phaseshift.cc -o build/wf_Xipi_phaseshift.o

build/kpi_phaseshifts.o : ${HFILES_CORAL} src/WaveFunctions/wfcommon/kpi_phaseshifts.cc
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INC} src/WaveFunctions/wfcommon/kpi_phaseshifts.cc -o build/kpi_phaseshifts.o

build/partwave.o : ${HFILES_CORAL} src/WaveFunctions/wfcommon/partwave.cc
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INC} src/WaveFunctions/wfcommon/partwave.cc -o build/partwave.o

build/planewave.o : ${HFILES_CORAL} src/WaveFunctions/wfcommon/planewave.cc
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INC} src/WaveFunctions/wfcommon/planewave.cc -o build/planewave.o

build/pipi_phaseshifts.o : ${HFILES_CORAL} src/WaveFunctions/wfcommon/pipi_phaseshifts.cc
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INC} src/WaveFunctions/wfcommon/pipi_phaseshifts.cc -o build/pipi_phaseshifts.o

build/reid93.o : ${HFILES_CORAL} src/WaveFunctions/wfcommon/reid93.cc
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INC} src/WaveFunctions/wfcommon/reid93.cc -o build/reid93.o

build/wfcommon.o : ${HFILES_CORAL} src/WaveFunctions/wfcommon/wfcommon.cc
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INC} src/WaveFunctions/wfcommon/wfcommon.cc -o build/wfcommon.o

build/wfsquarewell.o : ${HFILES_CORAL} src/WaveFunctions/wfcommon/wfsquarewell.cc
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INC} src/WaveFunctions/wfcommon/wfsquarewell.cc -o build/wfsquarewell.o

installdirs :
	mkdir -p ${MADAI_INSTALLDIR}/progdata;\

############## INSTALL HEADER FILES ###########

${MADAI_INSTALLDIR}/include/minimization.h : include/minimization.h
	cp -f include/minimization.h ${MADAI_INSTALLDIR}/include/

${MADAI_INSTALLDIR}/include/sfit.h : include/sfit.h
	cp -f include/sfit.h ${MADAI_INSTALLDIR}/include/

${MADAI_INSTALLDIR}/include/cfcalc.h : include/cfcalc.h
	cp -f include/cfcalc.h ${MADAI_INSTALLDIR}/include/

${MADAI_INSTALLDIR}/include/kernel.h : include/kernel.h
	cp -f include/kernel.h ${MADAI_INSTALLDIR}/include/

${MADAI_INSTALLDIR}/include/source2cf.h : include/source2cf.h
	cp -f include/source2cf.h ${MADAI_INSTALLDIR}/include/

${MADAI_INSTALLDIR}/include/sourcecalc.h : include/sourcecalc.h
	cp -f include/sourcecalc.h ${MADAI_INSTALLDIR}/include/

${MADAI_INSTALLDIR}/include/wavefunction.h : include/wavefunction.h
	cp -f include/wavefunction.h ${MADAI_INSTALLDIR}/include/

######################
uninstall:
	 rm -f ${MADAI_INSTALLDIR}/lib/libcoral.a ${MADAI_INSTALLDIR}/include/wavefunction.h ${MADAI_INSTALLDIR}/include/sourcecalc.h ${MADAI_INSTALLDIR}/include/source2cf.h ${MADAI_INSTALLDIR}/include/kernel.h ${MADAI_INSTALLDIR}/include/sfit.h ${MADAI_INSTALLDIR}/include/minimization.h ${MADAI_INSTALLDIR}/include/cfcalc.h

clean :
	rm -f build/*.o lib/*.a
