#include "msu_coral/wavefunction.h"
#include "msu_commonutils/constants.h"
#include "msu_commonutils/misc.h"
#include "msu_commonutils/sf.h"

CWaveFunction_nn_phaseshift::CWaveFunction_nn_phaseshift(string  parsfilename) : CWaveFunction() {
  int iq,ichannel;
  double q;
  ParsInit(parsfilename);
	
  m1=MNEUTRON;
  m2=m1;
	IDENTICAL=1;
  q1q2=0;
  nchannels=4;
	
  ellmax=1;
  InitArrays();
  CLog::Info("Arrays Initialized\n");
	
  ell[0]=0;
  ell[1]=ell[2]=ell[3]=1;
	
  InitWaves();
  CLog::Info("Partial Waves Initialized\n");
	
  // Channel weight is (2J+1)/[(2s1+1)*(2s2+1)]
  channelweight[0]=2*0.25;
  channelweight[1]=2*0.25;
  channelweight[2]=2*0.75;
  channelweight[3]=2*1.25;
  read_phaseshifts();
  //EffectiveRange(0,-16.75,2.7);
	
  for(ichannel=0;ichannel<nchannels;ichannel++){
    for(iq=0;iq<nqmax;iq++){
      q=qarray[iq];
      Wepsilon[ichannel][iq]=ddeltadq[ichannel][iq]
			-GetIW(ell[ichannel],epsilon,q,q1q2,eta[iq],delta[ichannel][iq])
			+GetIW(ell[ichannel],epsilon,q,q1q2,0.0,0.0);
      Wepsilon[ichannel][iq]=3.0*Wepsilon[ichannel][iq]
			/(4.0*PI*pow(epsilon,3));
    }
  }
}

double CWaveFunction_nn_phaseshift::CalcPsiSquared(int iq,double r,double ctheta){
  double psisquared,x,dpsi2,q,theta;
  double delta_1s0,delta_3p0,delta_3p1,delta_3p2;
  complex<double> psi,psisymm,psianti,psia,psib,hstar0,hstar1;
  complex<double> Xlm00,Xlm10,Xlm11;
  int ichannel;
	
  q=qarray[iq];
  if(iq>=nqmax){
    CLog::Fatal("iq too large! inside CWaveFunction_nn_phaseshift::CalcPsiSquared\n");
  }
  psia=planewave[iq]->planewave(r,ctheta);
  psib=planewave[iq]->planewave(r,-ctheta);
  psisymm=(1.0/sqrt(2.0))*(psia+psib);
  psianti=(1.0/sqrt(2.0))*(psia-psib);
  
  if(STRONG==1){
    if(r<epsilon){
      psisquared=0.25*real(psisymm*conj(psisymm))
			+0.75*real(psianti*conj(psianti));
      for(ichannel=0;ichannel<nchannels;ichannel++){
				dpsi2=channelweight[ichannel]*2.0*PI*Wepsilon[ichannel][iq]
				*pow(HBARC,3)/(q*q);
				//dpsi2=0.0;
				psisquared+=dpsi2;
      }

    }
    else{
      theta=acos(ctheta);
      x=q*r/HBARC;
      // Notation is (2S+1)-L-J
      delta_1s0=delta[0][iq];
      delta_3p0=delta[1][iq];
      delta_3p1=delta[2][iq];
      delta_3p2=delta[3][iq];
			//delta_3p0=delta_3p1=delta_3p2=0.0;
      hstar0=partwave[0][iq]->GetPhiIncoming(r)/x;
      hstar1=partwave[1][iq]->GetPhiIncoming(r)/x;
      // XlmLM 9s Y_{LM}*(1/2)*i^L*sqrt(4*PI*(2*L+1))*hstar_L
      Xlm00=0.5*sqrt(2.0)*sqrt(4.0*PI)*SpherHarmonics::Ylm(0,0,theta,0.0)*hstar0;
      Xlm10=ci*0.5*sqrt(2.0)*sqrt(12.0*PI)*SpherHarmonics::Ylm(1,0,theta,0.0)*hstar1;
      Xlm11=ci*0.5*sqrt(2.0)*sqrt(12.0*PI)*SpherHarmonics::Ylm(1,1,theta,0.0)*hstar1;
      
      // First do the case for S=0;
      psi=psisymm;
      // this refers to S=0, L=0, J=0 channel
      psi+=Xlm00*(Misc::ceiphi(-2.0*delta_1s0)-1.0);
      // S=0, L=1, J=1
      psisquared=0.25*real(psi*conj(psi));
      
      // Now let's do the case for S=1, M_S=+1
      psi=psianti;
      // S=1, L=1 and J=1,2,
      psi+=Xlm10*(0.5*Misc::ceiphi(-2.0*delta_3p1)
				+0.5*Misc::ceiphi(-2.0*delta_3p2)-1.0);
      psia=Xlm11*0.5*(Misc::ceiphi(-2.0*delta_3p1)-Misc::ceiphi(-2.0*delta_3p2));
      psisquared+=0.5*real(psi*conj(psi)+psia*conj(psia));
      // Term is doubled to account for M_S=-1
      
      // Now let's do the case with S=1, M_S=0;
      psi=psianti;
      // S=1, L=1, J=0,2
      psi+=Xlm10*((2.0/3.0)*Misc::ceiphi(-2.0*delta_3p2)
				+(1.0/3.0)*Misc::ceiphi(-2.0*delta_3p0)-1.0);
      psia=Xlm11*(Misc::ceiphi(-2.0*delta_3p2)-Misc::ceiphi(-2.0*delta_3p0))/3.0;
      psisquared+=0.25*real(psi*conj(psi)+2.0*psia*conj(psia));
    }
  }
  else psisquared=0.25*real(psisymm*conj(psisymm))
		+0.75*real(psianti*conj(psianti));
	
	psisquared*=RelativisticCorrection(r,iq);
  return psisquared;
	
}
void CWaveFunction_nn_phaseshift::read_phaseshifts(){
	//#include "nn_phaseshiftdat.cc"
	double delqdata=0.5;
	int nqdata=400;
	double data_delta[4][401]={{0,0.0414692,0.0827615,0.123704,0.164132,0.203893,0.242848,0.280875,0.317869,0.353745,0.388436,0.42189,0.454075,0.484973,0.514577,0.542893,0.569939,0.595737,0.620319,0.643718,0.665975,0.687131,0.70723,0.726315,0.744432,0.761624,0.777934,0.793407,0.808081,0.821999,0.835197,0.847712,0.85958,0.870833,0.881504,0.891622,0.901216,0.910312,0.918937,0.927113,0.934865,0.942213,0.949177,0.955777,0.96203,0.967954,0.973564,0.978876,0.983903,0.98866,0.993158,0.99741,1.00143,1.00522,1.0088,1.01218,1.01535,1.01835,1.02116,1.0238,1.02628,1.0286,1.03078,1.0328,1.03469,1.03645,1.03808,1.03959,1.04098,1.04225,1.04342,1.04449,1.04545,1.04632,1.04709,1.04778,1.04838,1.0489,1.04933,1.04969,1.04998,1.05019,1.05034,1.05042,1.05043,1.05038,1.05027,1.05011,1.04988,1.04961,1.04928,1.0489,1.04847,1.04799,1.04747,1.04691,1.0463,1.04565,1.04496,1.04423,1.04346,1.04266,1.04182,1.04095,1.04004,1.0391,1.03813,1.03713,1.0361,1.03504,1.03396,1.03285,1.03171,1.03055,1.02936,1.02815,1.02691,1.02566,1.02438,1.02308,1.02176,1.02042,1.01906,1.01768,1.01629,1.01488,1.01345,1.012,1.01054,1.00906,1.00756,1.00606,1.00453,1.003,1.00145,0.999884,0.998308,0.996719,0.995119,0.993507,0.991883,0.990249,0.988603,0.986948,0.985282,0.983606,0.981921,0.980227,0.978523,0.976811,0.97509,0.973361,0.971624,0.969879,0.968127,0.966367,0.9646,0.962826,0.961045,0.959258,0.957464,0.955664,0.953859,0.952047,0.95023,0.948407,0.946579,0.944745,0.942907,0.941064,0.939216,0.937363,0.935506,0.933644,0.931779,0.929909,0.928035,0.926158,0.924276,0.922391,0.920503,0.918611,0.916715,0.914817,0.912915,0.911011,0.909103,0.907193,0.905279,0.903363,0.901445,0.899524,0.8976,0.895674,0.893746,0.891816,0.889883,0.887948,0.886011,0.884073,0.882132,0.880189,0.878245,0.876299,0.874351,0.872401,0.87045,0.868498,0.866544,0.864588,0.862631,0.860673,0.858713,0.856752,0.85479,0.852826,0.850862,0.848896,0.846929,0.844961,0.842992,0.841022,0.839051,0.837079,0.835107,0.833133,0.831158,0.829183,0.827207,0.82523,0.823252,0.821274,0.819294,0.817315,0.815334,0.813353,0.811371,0.809389,0.807406,0.805422,0.803438,0.801453,0.799468,0.797482,0.795496,0.793509,0.791522,0.789535,0.787547,0.785558,0.78357,0.78158,0.779591,0.777601,0.775611,0.77362,0.771629,0.769638,0.767646,0.765654,0.763662,0.76167,0.759677,0.757684,0.755691,0.753698,0.751704,0.74971,0.747716,0.745722,0.743728,0.741733,0.739738,0.737744,0.735749,0.733753,0.731758,0.729763,0.727767,0.725772,0.723776,0.72178,0.719784,0.717788,0.715792,0.713796,0.7118,0.709804,0.707808,0.705812,0.703816,0.701819,0.699823,0.697827,0.695831,0.693835,0.691839,0.689843,0.687847,0.685851,0.683855,0.68186,0.679864,0.677869,0.675873,0.673878,0.671883,0.669888,0.667893,0.665899,0.663904,0.66191,0.659916,0.657922,0.655928,0.653935,0.651942,0.649949,0.647956,0.645964,0.643971,0.641979,0.639988,0.637996,0.636005,0.634015,0.632024,0.630034,0.628045,0.626055,0.624066,0.622078,0.62009,0.618102,0.616115,0.614128,0.612141,0.610155,0.60817,0.606184,0.6042,0.602216,0.600232,0.598249,0.596266,0.594284,0.592302,0.590321,0.588341,0.586361,0.584382,0.582403,0.580425,0.578448,0.576471,0.574494,0.572519,0.570544,0.56857,0.566596,0.564623,0.562651,0.56068,0.558709,0.556739,0.55477,0.552801,0.550833,0.548866,0.5469,0.544935,0.54297,0.541006,0.539043,0.537081,0.53512,0.533159,0.531199,0.529241,0.527283,0.525326,0.52337,0.521414,0.51946,0.517507,0.515554,0.513603,0.511652,0.509702,0.507754,0.505806,0.503859,0.501913,0.499969,0.498025,0.496082,0.49414,0.4922,0.49026,0.488321,0.486384},
		{0,3.74217e-08,2.99374e-07,1.01039e-06,2.39499e-06,4.67771e-06,8.08309e-06,1.28356e-05,1.91599e-05,2.72804e-05,3.74217e-05,4.98083e-05,6.46647e-05,8.22155e-05,0.000102685,0.000126298,0.000153279,0.000183853,0.000218243,0.000256676,0.000299374,0.000346562,0.000398466,0.00045531,0.000517318,0.000584714,0.000657724,0.000736572,0.000821481,0.000912678,0.00101039,0.00111483,0.00122623,0.00134482,0.00147082,0.00160446,0.00174595,0.00189552,0.0020534,0.00221982,0.00239499,0.00257914,0.0027725,0.00297529,0.00318773,0.00341005,0.00364248,0.00388523,0.00413854,0.00440263,0.00467771,0.00496403,0.00526179,0.00557123,0.0059081,0.00624709,0.00658609,0.00692509,0.00726408,0.00760308,0.00794208,0.00828107,0.00867434,0.00909151,0.00950868,0.00992585,0.010343,0.0107602,0.0111774,0.011627,0.0121152,0.0126033,0.0130914,0.0135795,0.0140676,0.0145557,0.0151102,0.0156707,0.0162312,0.0167918,0.0173523,0.0179128,0.0185109,0.0191132,0.0197154,0.0203177,0.0209199,0.0215421,0.0222105,0.022879,0.0235475,0.024216,0.024886,0.0255928,0.0262995,0.0270063,0.027713,0.0284231,0.0291799,0.0299367,0.0306934,0.0314502,0.0322188,0.0330108,0.0338028,0.0345948,0.0353868,0.0362085,0.0370368,0.0378651,0.0386935,0.0395374,0.0403983,0.0412592,0.0421201,0.0429886,0.0438743,0.04476,0.0456457,0.0465399,0.0474557,0.0483714,0.0492872,0.0502088,0.0511397,0.0520705,0.0530013,0.0539492,0.0549078,0.0558664,0.0568251,0.0577901,0.0587556,0.059721,0.0606963,0.0616875,0.0626786,0.0636697,0.0646532,0.0656359,0.0666187,0.0676073,0.0686012,0.069595,0.0705973,0.0716288,0.0726603,0.0736918,0.0747289,0.075766,0.0768031,0.0778171,0.0788266,0.079836,0.0808485,0.0818621,0.0828758,0.0839157,0.0849668,0.0860179,0.0870691,0.0881205,0.0891718,0.0901924,0.0912053,0.0922182,0.0932296,0.0942408,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437,0.0951437},
		{0,-2.2388e-08,-1.79104e-07,-6.04475e-07,-1.43283e-06,-2.79849e-06,-4.8358e-06,-7.67907e-06,-1.14626e-05,-1.63208e-05,-2.2388e-05,-2.97984e-05,-3.86864e-05,-4.91863e-05,-6.14325e-05,-7.55593e-05,-9.1701e-05,-0.000109992,-0.000130567,-0.000153559,-0.000179104,-0.000207335,-0.000238387,-0.000272394,-0.000309491,-0.000349812,-0.000393491,-0.000440662,-0.00049146,-0.00054602,-0.000604475,-0.000666959,-0.000733608,-0.000804556,-0.000879936,-0.000959883,-0.00104453,-0.00113402,-0.00122847,-0.00132803,-0.00143283,-0.001543,-0.00165868,-0.00178,-0.0019071,-0.0020401,-0.00217915,-0.00232438,-0.00247593,-0.00263392,-0.00279849,-0.00296978,-0.00314792,-0.00333305,-0.00353098,-0.00372981,-0.00392864,-0.00412748,-0.00432631,-0.00452514,-0.00472398,-0.00492281,-0.00515174,-0.00539392,-0.00563611,-0.00587829,-0.00612048,-0.00636266,-0.00660484,-0.00686471,-0.00714551,-0.0074263,-0.0077071,-0.00798789,-0.00826869,-0.00854948,-0.00886502,-0.00918374,-0.00950246,-0.00982118,-0.0101399,-0.0104586,-0.0107987,-0.0111411,-0.0114835,-0.011826,-0.0121684,-0.0125208,-0.0128963,-0.0132719,-0.0136475,-0.014023,-0.0143995,-0.0147966,-0.0151936,-0.0155907,-0.0159878,-0.0163865,-0.0168083,-0.0172301,-0.0176518,-0.0180736,-0.0185019,-0.0189433,-0.0193847,-0.0198261,-0.0202674,-0.0207236,-0.0211832,-0.0216428,-0.0221024,-0.0225706,-0.0230483,-0.0235259,-0.0240036,-0.0244853,-0.024976,-0.0254668,-0.0259575,-0.026453,-0.0269604,-0.0274678,-0.0279752,-0.0284862,-0.0290028,-0.0295194,-0.0300361,-0.0305622,-0.0310943,-0.0316263,-0.0321584,-0.0326963,-0.0332345,-0.0337727,-0.0343164,-0.034869,-0.0354215,-0.0359741,-0.0365274,-0.0370807,-0.0376341,-0.0381917,-0.0387532,-0.0393147,-0.0398803,-0.0404603,-0.0410403,-0.0416203,-0.042207,-0.0427937,-0.0433804,-0.0439624,-0.0445435,-0.0451245,-0.0457097,-0.0462964,-0.0468831,-0.0474822,-0.0480865,-0.0486908,-0.0492987,-0.0499079,-0.0505172,-0.0511198,-0.0517207,-0.0523216,-0.0529261,-0.0535308,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706,-0.0540706},
	{0,4.68005e-09,3.74404e-08,1.26361e-07,2.99523e-07,5.85006e-07,1.01089e-06,1.60526e-06,2.39618e-06,3.41175e-06,4.68005e-06,6.22914e-06,8.08712e-06,1.02821e-05,1.2842e-05,1.57952e-05,1.91695e-05,2.29931e-05,2.7294e-05,3.21004e-05,3.74404e-05,4.33419e-05,4.98331e-05,5.69421e-05,6.4697e-05,7.31257e-05,8.22565e-05,9.21173e-05,0.000102736,0.000114142,0.000126361,0.000139423,0.000153356,0.000168187,0.000183945,0.000200657,0.000218352,0.000237058,0.000256803,0.000277616,0.000299523,0.000322553,0.000346735,0.000372096,0.000398665,0.000426469,0.000455537,0.000485896,0.000517576,0.000550603,0.000585006,0.000620813,0.000658052,0.000696751,0.00074287,0.000789662,0.000836454,0.000883246,0.000930038,0.00097683,0.00102362,0.00107041,0.00112682,0.00118745,0.00124809,0.00130873,0.00136936,0.00143,0.00149064,0.00155764,0.00163218,0.00170672,0.00178126,0.0018558,0.00193034,0.00200488,0.00209526,0.00218708,0.0022789,0.00237072,0.00246254,0.00255436,0.00265233,0.00275098,0.00284963,0.00294828,0.00304693,0.00315183,0.00327127,0.0033907,0.00351014,0.00362958,0.0037493,0.00387554,0.00400178,0.00412803,0.00425427,0.00438188,0.00452852,0.00467515,0.00482179,0.00496843,0.00511737,0.00527088,0.00542439,0.00557789,0.0057314,0.00590107,0.0060744,0.00624773,0.00642107,0.00659768,0.00677786,0.00695804,0.00713822,0.00732431,0.00752375,0.0077232,0.00792264,0.00812398,0.00833019,0.0085364,0.0087426,0.00895606,0.00918094,0.00940581,0.00963069,0.0098597,0.0100913,0.0103229,0.0105545,0.0108032,0.0110528,0.0113023,0.0115545,0.0118107,0.012067,0.0123232,0.0125972,0.0128727,0.0131483,0.0134304,0.0137184,0.0140064,0.0142947,0.0145843,0.0148738,0.0151634,0.0154636,0.0157638,0.016064,0.0163814,0.0167022,0.017023,0.0173526,0.0176852,0.0180179,0.0183519,0.0186865,0.0190211,0.0193628,0.0197073,0.0200518,0.0204113,0.0207744,0.0211376,0.0215109,0.0218851,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191,0.0222191}};
	double data_ddeltadq[4][401]={{0.0833486,0.08282,0.0822913,0.081424,0.0802379,0.0787593,0.0770194,0.0750527,0.0728959,0.0705859,0.0681586,0.0656482,0.0630863,0.060501,0.0579171,0.0553558,0.0528348,0.0503686,0.0479685,0.0456432,0.043399,0.0412401,0.039169,0.0371869,0.0352937,0.0334883,0.031769,0.0301335,0.0285791,0.0271028,0.0257015,0.0243717,0.0231103,0.0219138,0.0207789,0.0197025,0.0186815,0.0177128,0.0167935,0.015921,0.0150925,0.0143056,0.0135579,0.0128473,0.0121715,0.0115287,0.010917,0.0103345,0.00977967,0.00925095,0.00874688,0.00826609,0.00780729,0.00736929,0.00695094,0.00655118,0.00616903,0.00580355,0.00545384,0.0051191,0.00479854,0.00449144,0.0041971,0.00391488,0.00364419,0.00338443,0.00313508,0.00289563,0.0026656,0.00244454,0.00223203,0.00202766,0.00183105,0.00164186,0.00145973,0.00128435,0.00111541,0.000952639,0.000795753,0.000644501,0.000498637,0.000357931,0.000222165,9.11302e-05,-3.53706e-05,-0.000157523,-0.000275506,-0.00038949,-0.000499634,-0.000606091,-0.000709008,-0.000808523,-0.000904768,-0.00099787,-0.00108795,-0.00117512,-0.00125949,-0.00134116,-0.00142024,-0.00149682,-0.00157098,-0.00164283,-0.00171244,-0.00177988,-0.00184525,-0.0019086,-0.00197001,-0.00202955,-0.00208728,-0.00214326,-0.00219756,-0.00225022,-0.0023013,-0.00235086,-0.00239894,-0.0024456,-0.00249088,-0.00253482,-0.00257746,-0.00261886,-0.00265904,-0.00269805,-0.00273593,-0.0027727,-0.0028084,-0.00284307,-0.00287674,-0.00290943,-0.00294118,-0.00297202,-0.00300197,-0.00303106,-0.00305932,-0.00308677,-0.00311343,-0.00313933,-0.00316449,-0.00318893,-0.00321267,-0.00323574,-0.00325814,-0.00327991,-0.00330106,-0.00332161,-0.00334157,-0.00336096,-0.0033798,-0.0033981,-0.00341588,-0.00343315,-0.00344993,-0.00346623,-0.00348207,-0.00349745,-0.0035124,-0.00352691,-0.00354102,-0.00355471,-0.00356802,-0.00358094,-0.00359349,-0.00360568,-0.00361752,-0.00362902,-0.00364018,-0.00365103,-0.00366156,-0.00367178,-0.00368171,-0.00369135,-0.00370071,-0.0037098,-0.00371862,-0.00372719,-0.0037355,-0.00374357,-0.00375141,-0.00375901,-0.00376639,-0.00377356,-0.00378051,-0.00378726,-0.00379381,-0.00380016,-0.00380633,-0.00381231,-0.00381811,-0.00382375,-0.00382921,-0.00383451,-0.00383965,-0.00384464,-0.00384948,-0.00385417,-0.00385873,-0.00386314,-0.00386742,-0.00387158,-0.0038756,-0.00387951,-0.0038833,-0.00388697,-0.00389053,-0.00389399,-0.00389734,-0.00390058,-0.00390373,-0.00390678,-0.00390974,-0.00391261,-0.00391539,-0.00391809,-0.00392071,-0.00392324,-0.0039257,-0.00392808,-0.00393039,-0.00393263,-0.0039348,-0.0039369,-0.00393894,-0.00394092,-0.00394284,-0.00394469,-0.00394649,-0.00394824,-0.00394993,-0.00395157,-0.00395316,-0.0039547,-0.0039562,-0.00395764,-0.00395905,-0.00396041,-0.00396173,-0.003963,-0.00396424,-0.00396544,-0.0039666,-0.00396773,-0.00396882,-0.00396988,-0.0039709,-0.00397189,-0.00397285,-0.00397378,-0.00397468,-0.00397555,-0.0039764,-0.00397721,-0.003978,-0.00397876,-0.0039795,-0.00398021,-0.0039809,-0.00398156,-0.0039822,-0.00398282,-0.00398341,-0.00398399,-0.00398454,-0.00398507,-0.00398558,-0.00398607,-0.00398654,-0.00398699,-0.00398742,-0.00398783,-0.00398822,-0.00398859,-0.00398895,-0.00398928,-0.0039896,-0.0039899,-0.00399018,-0.00399045,-0.00399069,-0.00399092,-0.00399113,-0.00399133,-0.0039915,-0.00399167,-0.00399181,-0.00399193,-0.00399204,-0.00399213,-0.00399221,-0.00399226,-0.0039923,-0.00399233,-0.00399233,-0.00399232,-0.00399229,-0.00399224,-0.00399218,-0.0039921,-0.003992,-0.00399188,-0.00399175,-0.0039916,-0.00399143,-0.00399124,-0.00399103,-0.00399081,-0.00399057,-0.00399031,-0.00399003,-0.00398973,-0.00398941,-0.00398908,-0.00398872,-0.00398835,-0.00398796,-0.00398755,-0.00398711,-0.00398666,-0.00398619,-0.0039857,-0.0039852,-0.00398467,-0.00398412,-0.00398355,-0.00398296,-0.00398235,-0.00398172,-0.00398107,-0.0039804,-0.00397971,-0.00397899,-0.00397826,-0.00397751,-0.00397673,-0.00397594,-0.00397512,-0.00397428,-0.00397343,-0.00397255,-0.00397165,-0.00397072,-0.00396978,-0.00396882,-0.00396783,-0.00396682,-0.0039658,-0.00396475,-0.00396367,-0.00396258,-0.00396147,-0.00396033,-0.00395918,-0.003958,-0.0039568,-0.00395558,-0.00395434,-0.00395308,-0.0039518,-0.00395049,-0.00394917,-0.00394782,-0.00394645,-0.00394507,-0.00394366,-0.00394223,-0.00394078,-0.00393931,-0.00393782,-0.00393632,-0.00393479,-0.00393324,-0.00393167,-0.00393008,-0.00392847,-0.00392685,-0.0039252,-0.00392354,-0.00392185,-0.00392015,-0.00391843,-0.00391669,-0.00391493,-0.00391316,-0.00391137,-0.00390956,-0.00390773,-0.00390588,-0.00390402,-0.00390214,-0.00390025,-0.00389834,-0.00389641,-0.00389447,-0.00389251,-0.00389053,-0.00388855,-0.00388654,-0.00388452,-0.00388249,-0.00388044,-0.00387838,-0.00387631,-0.00387422},
		{0,2.2453e-07,8.98121e-07,2.02077e-06,3.59248e-06,5.61326e-06,8.08309e-06,1.1002e-05,1.43699e-05,1.8187e-05,2.2453e-05,2.71682e-05,3.23324e-05,3.79456e-05,4.40079e-05,5.05193e-05,5.74797e-05,6.48892e-05,7.27478e-05,8.10554e-05,8.98121e-05,9.90178e-05,0.000108673,0.000118777,0.000129329,0.000140331,0.000151782,0.000163683,0.000176032,0.00018883,0.000202077,0.000215774,0.000229919,0.000244513,0.000259557,0.00027505,0.000290991,0.000307382,0.000324222,0.000341511,0.000359248,0.000377435,0.000396071,0.000415156,0.000434691,0.000454674,0.000475106,0.000495987,0.000517318,0.000539097,0.000561326,0.000584003,0.00060713,0.000630705,0.000677994,0.000677994,0.000677994,0.000677994,0.000677994,0.000677994,0.000677994,0.000677994,0.000834343,0.000834343,0.000834343,0.000834343,0.000834343,0.000834343,0.000834343,0.000976233,0.000976233,0.000976233,0.000976233,0.000976233,0.000976233,0.000976233,0.00112106,0.00112106,0.00112106,0.00112106,0.00112106,0.00112106,0.0012045,0.0012045,0.0012045,0.0012045,0.0012045,0.00133694,0.00133694,0.00133694,0.00133694,0.00133694,0.00141352,0.00141352,0.00141352,0.00141352,0.00141352,0.00151352,0.00151352,0.00151352,0.00151352,0.00151352,0.00158403,0.00158403,0.00158403,0.00158403,0.00158403,0.00165667,0.00165667,0.00165667,0.00165667,0.0017218,0.0017218,0.0017218,0.0017218,0.00177138,0.00177138,0.00177138,0.00177138,0.00183156,0.00183156,0.00183156,0.00183156,0.00186163,0.00186163,0.00186163,0.00186163,0.00191724,0.00191724,0.00191724,0.00191724,0.00193082,0.00193082,0.00193082,0.00198224,0.00198224,0.00198224,0.00198224,0.00196553,0.00196553,0.00196553,0.00198767,0.00198767,0.00198767,0.00206298,0.00206298,0.00206298,0.00206298,0.00207422,0.00207422,0.00207422,0.00201882,0.00201882,0.00201882,0.0020273,0.0020273,0.0020273,0.0021022,0.0021022,0.0021022,0.00210269,0.00210269,0.00210269,0.00202581,0.00202581,0.00202581,0.00202246,0.00202246,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
		{0,-1.34328e-07,-5.37311e-07,-1.20895e-06,-2.14924e-06,-3.35819e-06,-4.8358e-06,-6.58206e-06,-8.59697e-06,-1.08805e-05,-1.34328e-05,-1.62537e-05,-1.93432e-05,-2.27014e-05,-2.63282e-05,-3.02237e-05,-3.43879e-05,-3.88207e-05,-4.35222e-05,-4.84923e-05,-5.37311e-05,-5.92385e-05,-6.50146e-05,-7.10594e-05,-7.73728e-05,-8.39548e-05,-9.08055e-05,-9.79249e-05,-0.000105313,-0.00011297,-0.000120895,-0.000129089,-0.000137552,-0.000146283,-0.000155283,-0.000164551,-0.000174089,-0.000183895,-0.000193969,-0.000204312,-0.000214924,-0.000225805,-0.000236954,-0.000248372,-0.000260058,-0.000272014,-0.000284237,-0.00029673,-0.000309491,-0.000322521,-0.000335819,-0.000349386,-0.000363222,-0.000377327,-0.000397667,-0.000397667,-0.000397667,-0.000397667,-0.000397667,-0.000397667,-0.000397667,-0.000397667,-0.000484368,-0.000484368,-0.000484368,-0.000484368,-0.000484368,-0.000484368,-0.000484368,-0.00056159,-0.00056159,-0.00056159,-0.00056159,-0.00056159,-0.00056159,-0.00056159,-0.000637441,-0.000637441,-0.000637441,-0.000637441,-0.000637441,-0.000637441,-0.000684848,-0.000684848,-0.000684848,-0.000684848,-0.000684848,-0.000751136,-0.000751136,-0.000751136,-0.000751136,-0.000751136,-0.000794161,-0.000794161,-0.000794161,-0.000794161,-0.000794161,-0.000843524,-0.000843524,-0.000843524,-0.000843524,-0.000843524,-0.000882743,-0.000882743,-0.000882743,-0.000882743,-0.000882743,-0.000919149,-0.000919149,-0.000919149,-0.000919149,-0.000955319,-0.000955319,-0.000955319,-0.000955319,-0.000981507,-0.000981507,-0.000981507,-0.000981507,-0.00101481,-0.00101481,-0.00101481,-0.00101481,-0.00103326,-0.00103326,-0.00103326,-0.00103326,-0.00106413,-0.00106413,-0.00106413,-0.00106413,-0.00107639,-0.00107639,-0.00107639,-0.00110511,-0.00110511,-0.00110511,-0.00110511,-0.0011067,-0.0011067,-0.0011067,-0.00112293,-0.00112293,-0.00112293,-0.00115999,-0.00115999,-0.00115999,-0.00115999,-0.00117346,-0.00117346,-0.00117346,-0.0011621,-0.0011621,-0.0011621,-0.00117339,-0.00117339,-0.00117339,-0.00120861,-0.00120861,-0.00120861,-0.00121855,-0.00121855,-0.00121855,-0.00120184,-0.00120184,-0.00120184,-0.00120939,-0.00120939,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
	{0,2.80803e-08,1.12321e-07,2.52722e-07,4.49284e-07,7.02007e-07,1.01089e-06,1.37593e-06,1.79714e-06,2.2745e-06,2.80803e-06,3.39771e-06,4.04356e-06,4.74557e-06,5.50373e-06,6.31806e-06,7.18855e-06,8.1152e-06,9.09801e-06,1.0137e-05,1.12321e-05,1.23834e-05,1.35909e-05,1.48545e-05,1.61742e-05,1.75502e-05,1.89823e-05,2.04705e-05,2.20149e-05,2.36155e-05,2.52722e-05,2.69851e-05,2.87542e-05,3.05794e-05,3.24608e-05,3.43983e-05,3.6392e-05,3.84419e-05,4.05479e-05,4.27101e-05,4.49284e-05,4.72029e-05,4.95336e-05,5.19204e-05,5.43634e-05,5.68626e-05,5.94179e-05,6.20293e-05,6.4697e-05,6.74207e-05,7.02007e-05,7.30368e-05,7.59291e-05,7.88775e-05,9.35838e-05,9.35838e-05,9.35838e-05,9.35838e-05,9.35838e-05,9.35838e-05,9.35838e-05,9.35838e-05,0.000121273,0.000121273,0.000121273,0.000121273,0.000121273,0.000121273,0.000121273,0.000149081,0.000149081,0.000149081,0.000149081,0.000149081,0.000149081,0.000149081,0.000183643,0.000183643,0.000183643,0.000183643,0.000183643,0.000183643,0.000197301,0.000197301,0.000197301,0.000197301,0.000197301,0.000238874,0.000238874,0.000238874,0.000238874,0.000238874,0.000252487,0.000252487,0.000252487,0.000252487,0.000252487,0.000293275,0.000293275,0.000293275,0.000293275,0.000293275,0.000307014,0.000307014,0.000307014,0.000307014,0.000307014,0.000346663,0.000346663,0.000346663,0.000346663,0.000360357,0.000360357,0.000360357,0.000360357,0.000398881,0.000398881,0.000398881,0.000398881,0.000412414,0.000412414,0.000412414,0.000412414,0.000449752,0.000449752,0.000449752,0.000449752,0.000463186,0.000463186,0.000463186,0.000463186,0.000499162,0.000499162,0.000499162,0.00051248,0.00051248,0.00051248,0.00051248,0.000551135,0.000551135,0.000551135,0.000575932,0.000575932,0.000575932,0.00057913,0.00057913,0.00057913,0.00057913,0.000600474,0.000600474,0.000600474,0.00064164,0.00064164,0.00064164,0.000665301,0.000665301,0.000665301,0.000669196,0.000669196,0.000669196,0.000688993,0.000688993,0.000688993,0.000726335,0.000726335,0.000726335,0.000748247,0.000748247,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}};
	//
  int iqdata,iq,ichannel;
  double w1,w2,q;
  for(ichannel=0;ichannel<4;ichannel++){
	  for(iq=0;iq<nqmax;iq++){
		  q=qarray[iq];
		  iqdata=int(floor(q)/delqdata);
		  if(iqdata<nqdata){
			  w1=(delqdata*double(iqdata+1)-q)/delqdata;
			  w2=1.0-w1;
			  delta[ichannel][iq]=w1*data_delta[ichannel][iqdata]
				  +w2*data_delta[ichannel][iqdata+1];
			  ddeltadq[ichannel][iq]=w1*data_ddeltadq[ichannel][iqdata]
				  +w2*data_ddeltadq[ichannel][iqdata+1];
		  }
		  else{
			  delta[ichannel][iq]=0.0;
			  ddeltadq[ichannel][iq]=0.0;
		  }			
	  }
  }
}
