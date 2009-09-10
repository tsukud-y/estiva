#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include "estiva.h"
#include "confary.h"
#define is(type) if(siz(argv[0])==sizeof(type))
#define eqn(str) if(!strcmp(fmt,str))
#define illegaltype(fmt) \
else{fprintf(stderr,"matf: illegal type %s\n",fmt); exit(1);}

static void args(void **argv,...)
{ int i; va_list ap;
  va_start(ap,argv);
  forall(0,i,dim1(argv)) *va_arg(ap,void **) = argv[i];
  va_end(ap);
}

int matutil_matf(void *argv0,char *fmt,...)
{ static void **argv;

  { void *buf[9];char c;int i,argc;static double tmp;va_list ap;
    va_start(ap,fmt);
    
    for(i=0,argc=0;'\0'!=(c=fmt[i]);i++)if(c=='%')switch(fmt[i+1]){
    case 'p':
      buf[++argc]=va_arg(ap,void *);
      break;
    case 'f':
      tmp=va_arg(ap,double);
      buf[++argc]= &tmp;
      break;
    default:
      fprintf(stderr,"matf: illegal format %s\n",fmt); exit(1);
    }
    va_end(ap);
    if(ary1(argv,argc+1)==NULL){
      fprintf(stderr,"matf: not enough memory\n"); exit(1);
    }
    argv[0]=argv0;forall(1,i,argc) argv[i]=buf[i];
  }
  eqn("%p+%f*%p"){
    is(float){
      float **A,*Ai,**K,*Ki,*M,Mi; 
#include "matutil/eqn1.c"
    }else is(double){
      double **A,*Ai,**K,*Ki,*M,Mi; 
#include "matutil/eqn1.c"
    }illegaltype(fmt);
  }
  else eqn("%p*(%p+%f*%p)"){
    is(float){
      float *R, *M, *F, *U; 
#include "matutil/eqn2.c"
    }else is(double){
      double *R, *M, *F, *U; 
#include "matutil/eqn2.c"
    }illegaltype(fmt);
  }
  else eqn("%pt*%p"){
    is(float){
      float **T,*Ti,Tij,**A,**B;
#include "matutil/eqn3.c"
    }else is(double){
      double **T,*Ti,Tij,**A,**B;
#include "matutil/eqn3.c"
    }illegaltype(fmt);
  }
  else eqn("-(%p*%p+%p*%p)"){
    is(float){
      float **B,*Bi,Bij,**HXtAX,*HXtAXi,**HX,**HYtAY,*HYtAYi,**HY;
#include "matutil/eqn4.c"
    }else is(double){
      double **B,*Bi,Bij,**HXtAX,*HXtAXi,**HX,**HYtAY,*HYtAYi,**HY;
#include "matutil/eqn4.c"
    }illegaltype(fmt);
  }
  else eqn("%p*%p+%p*%p"){
    is(float){
      float *P,Pi,**HXtAX,*HXtAXi,*Rx,**HYtAY,*HYtAYi,*Ry;
#include "matutil/eqn5.c"
    }else is(double){
      double *P,Pi,**HXtAX,*HXtAXi,*Rx,**HYtAY,*HYtAYi,*Ry;
#include "matutil/eqn5.c"
    }illegaltype(fmt);
  }
  else eqn("%p*(%p+%p*%p)"){
    is(float){
      float *U,Ui,**A,*Ai,*R,**H,*Hj,*P; static float *w,wj;
#include "matutil/eqn6.c"
    }else is(double){
      double *U,Ui,**A,*Ai,*R,**H,*Hj,*P; static float *w,wj;
#include "matutil/eqn6.c"
    }illegaltype(fmt);
  }
  else{ fprintf(stderr,"matf: not supported %s\n",fmt); exit(1);}
  return 1;
}
