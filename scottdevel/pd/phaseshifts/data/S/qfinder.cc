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
	const int NQMAX=18;
	const double MP=938.272,MN=939.565;
	double MD=MP+MN-2.224575;
	double elab,q,delta,mu,vrel;
	int iq;
  FILE *fptr=fopen("delta_2S_12.txt","r");
	FILE *output=fopen("qarray.txt","w");
	
	fprintf(output,"%d\n",NQMAX);
	for(iq=0;iq<NQMAX;iq++){
		fscanf(fptr,"%lf %lf",&elab,&delta);
		mu=MP*MD/(MP+MD);
		q=mu*sqrt(2.0*elab/MP);
		vrel=sqrt(2.0*elab/MP);
		fprintf(output,"%lf\n",q);
		printf("q=%g, elab=%g=?%g\n",q,0.5*MP*vrel*vrel,elab);
	}
	fclose(fptr);
	fclose(output);
	
  return 0;
}


