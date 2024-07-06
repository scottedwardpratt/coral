#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <complex>
#include <string>
#include <cstring>
#include <boost/math/special_functions.hpp>

const double PI=4.0*atan(1.0);
const double HBARC=197.3269602;

using namespace std;
using namespace boost::math;

int main(int argc,char *argv[]){
	const int NQMAX=7;
	double qarray[NQMAX]={0.0};
	double elab[NQMAX]={0.0};
	double deltaS12[NQMAX]={0.0};
	double deltaS32[NQMAX]={0.0};
	const double MP=938.272,MN=939.565;
	double MD=MP+MN-2.224575;
	double q,delta;
	int iq;
  FILE *fptr;
	FILE *output=fopen("qarray.txt","w");
	fptr=fopen("delta_2D_32.txt","r");
	for(iq=0;iq<NQMAX;iq++){
		fscanf(fptr,"%lf %lf",&elab[iq],&delta);
		qarray[iq]=(MP*MD/(MP+MD))*sqrt(2.0*elab[iq]/MP);
		fprintf(output,"%lf\n",qarray[iq]);
		deltaS12[iq]+=(2.0/5.0)*delta;
	}
	fclose(fptr);
	fclose(output);
	
	fptr=fopen("delta_2D_52.txt","r");
	for(iq=0;iq<NQMAX;iq++){
		fscanf(fptr,"%lf %lf",&elab[iq],&delta);
		deltaS12[iq]+=(3.0/5.0)*delta;

	}
	fclose(fptr);
	
	fptr=fopen("delta_4D_12.txt","r");
	for(iq=0;iq<NQMAX;iq++){
		fscanf(fptr,"%lf %lf",&elab[iq],&delta);
		deltaS32[iq]+=(1.0/10.0)*delta;
	}
	fclose(fptr);
	
	fptr=fopen("delta_4D_32.txt","r");
	for(iq=0;iq<NQMAX;iq++){
		fscanf(fptr,"%lf %lf",&elab[iq],&delta);
		deltaS32[iq]+=(2.0/10.0)*delta;
	}
	fclose(fptr);
	
	fptr=fopen("delta_4D_52.txt","r");
	for(iq=0;iq<NQMAX;iq++){
		fscanf(fptr,"%lf %lf",&elab[iq],&delta);
		deltaS32[iq]+=(3.0/10.0)*delta;
	}
	fclose(fptr);
	
	fptr=fopen("delta_4D_72.txt","r");
	for(iq=0;iq<NQMAX;iq++){
		fscanf(fptr,"%lf %lf",&elab[iq],&delta);
		deltaS32[iq]+=(4.0/10.0)*delta;
	}
	fclose(fptr);
	
	output=fopen("deltabar.txt","w");
	fprintf(output,"# q(MeV/c)    delta_S=1/2    delta_S=3/2\n");
	for(iq=0;iq<NQMAX;iq++){
		fprintf(output,"%6.3f %8.4f %8.4f\n",elab[iq],deltaS12[iq],deltaS32[iq]);
	}
	fclose(output);
	
  return 0;
}


