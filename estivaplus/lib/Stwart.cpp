#include "estivaplus.h"
#include "stwart.h"
#include <f2c.h>

void Stwart(vector<long>&IA,vector<long>&MA,
            vector<long>&JROW,vector<long>&JCOL)
{
  static long M = JROW.size();
  vector<long> IR(M),IC(M),R(M),C(M),IP(M),JP(M),IW(M);
  static long KERNS, MEND, LG, IER, L = MA[M];
  // printf("M = %ld\n",M);
  // forVector(MA,i) printf("MA[%ld] = %ld\n",i,MA[i]);
  // forVector(IA,i) printf("IA[%ld] = %ld\n",i,IA[i]);
  // printf("L=%ld\n",L);
  // forVector(R,i) printf("R[%ld] = %ld\n",i,R[i]);
  // forVector(C,i) printf("C[%ld] = %ld\n",i,C[i]);
  // forVector(IR,i) printf("IR[%ld] = %ld\n",i,IR[i]);
  // forVector(IC,i) printf("IC[%ld] = %ld\n",i,IC[i]);
  // forVector(JROW,i) printf("JROW[%ld] = %ld\n",i,JROW[i]);
  // forVector(JCOL,i) printf("JCOL[%ld] = %ld\n",i,JCOL[i]);
  // forVector(IP,i) printf("IP[%ld] = %ld\n",i,IP[i]);
  // forVector(JP,i) printf("JP[%ld] = %ld\n",i,JP[i]);
  // forVector(IW,i) printf("IW[%ld] = %ld\n",i,IW[i]);
  printf("KERNS = %ld, MEND = %ld, LG = %ld, IER = %ld\n",
	 KERNS, MEND, LG, IER);
  

  stwart_(&IA[0],&L,&MA[0],&M,&R[0],&C[0],&IR[0],&IC[0],
          &JROW[0],&JCOL[0],&IP[0],&JP[0],&KERNS,&MEND,&IW[0],&LG,&IER);
}
