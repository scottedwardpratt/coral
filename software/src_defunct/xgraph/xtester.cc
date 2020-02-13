//#include <stream>
#include <cstdlib>
#include <cmath>
//#include <dlfcn.h>
#include <X11/Xlib.h>
#include "xgraph.h"

int main(){
  double xmin,ymin,xmax,ymax;
  double x,y;
  char dummy[100];
  int window_xoff=20,window_yoff=100,window_width=500,window_height=500;
  double x1,y1,x2,y2;
  double symbolsize=0.02;
  char title[20];
  int i,npts;

  xmin=0.0;
  ymin=-1.1;
  xmax=12.0;
  ymax=1.1;
  //printf("Enter width, height, xoff, yoff :");
  //scanf("%d %d %d %d",&window_width,&window_height,&window_xoff,&window_yoff);
  CXGraph xgraph(window_width,window_height,window_xoff,window_yoff);
  window_xoff+=window_width+20;
  CXGraph ygraph(window_width,window_height,window_xoff,window_yoff);

  xgraph.setaxes(xmin,ymin,xmax,ymax);
  xgraph.drawaxes();

  x=0.5*(xmax+xmin);
  y=0.5*(ymin+ymax);
  xgraph.drawtext(title,x,y);
  xgraph.setcolor("cyan");

  
  ygraph.setaxes(xmin,ymin,xmax,ymax);
  ygraph.drawaxes();
  x1=xmin+0.2*(xmax-xmin);
  y1=ymin+0.2*(ymax-ymin);
  x2=xmin+0.4*(xmax-xmin);
  y2=ymin+0.6*(ymax-ymin);
  ygraph.drawarrow(x1,y1,x2,y2,.03);

  npts=50;
  for(i=0;i<=npts;i++){
    x=double(i)*xmax/double(npts);
    y=sin(x);
    xgraph.drawpoint(x,y);
    xgraph.drawcircle(x,y,symbolsize);
    xgraph.drawsquare(x,y,symbolsize);
    xgraph.drawdiamond(x,y,symbolsize);
    xgraph.drawuptriangle(x,y,symbolsize);
    xgraph.drawdowntriangle(x,y,symbolsize);
  }


  //pause();
  printf("Enter anything :");
  scanf("%s",&dummy);
  xgraph.closedisplay();
  ygraph.closedisplay();
  return 0;
}

#include "xgraph.cc"
