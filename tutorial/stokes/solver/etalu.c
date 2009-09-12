#include <stdio.h>
#include "ary.h"
#include "spm.h"

#define A(i,j) (*spm_double(pA,i,j))

typedef long*    I_  ;
typedef long**   I__ ;
typedef double*  D_  ;
typedef double** D__ ;

#define D_  static D_
#define D__ static D__
#define I_  static I_
#define I__ static I__


int etalu(void* pA,double* b)
{
  /* (i)  引数の型と種類                                                */
  D_     A, E      ;/*    D  一次元配列  A(L), E(L2)                    */
  I_     IA, IE    ;/*    I  一次元配列  IA(L), IE(L2)                  */
  D_     B, X, WE  ;/*    D  一次元配列  B(M), X(M), WE(M)              */
  I_     NW, IP    ;/*    I  一次元配列  NWORK(0:M), IP(M)              */
  I_     NETA      ;/*    I  一次元配列  NETA(-1:2*M)                   */
  static long   L, L2, M, NUM, IS, IER                                  ;
  static double EPS                                                     ;

  static long   i, j, k                                                 ;

  M   = dim1(pA) ;

  ary1( NW, M+1 );

  for(j=1;j<=M;j++)
    for(i=1;i<=M;i++)
      if(A(i,j) != 0.0) NW[j]++;
  NW[0] = 0;


  for(L=0,j=1;j<=M;j++) L += NW[j];

  L2 = L*2;
  

  k=0;
  for(i=1;i<=M;i++) for(j=1;j<=M;j++) if( A(i,j) != 0.0) k++;
  fprintf(stderr,"M L L2 k %d %d %d %d\n",M,L,L2,k); fflush(stderr);  

  ary1(  A, L ); ary1(  E, L2 );
  ary1( IA, L ); ary1( IE, L2 );
  ary1(  B, M ); ary1(  X, M  ); ary1( WE, M );
                 ary1( IP, M  );
  ary1(NETA, 2*M+2 );

  /* (ii) 主プログラム → サブルーチン                                  */
  /*      サブルーチンをCALLするときには, つぎの値を与える.             */
  /*  A   : 配列Aの各列の非ゼロ要素を列ごとに並べておく.                */
  /*  IA  : 配列Aの各要素に対応した行番号を列ごとに並べておく.          */
  /*  L   : 主プログラムにおける配列Aの大きさ                           */
  /*  L2  : 主プログラムにおける配列Eの大きさ  L2>L                     */
  /*  B   : 右辺ベクトルbを入れておく.                                  */
  /*  NW  : 行列Aの各列の非ゼロ要素数の累和を入れる. NW(0)=0とする.     */
  /*  M   : 行列Aの列数を入れておく.                                    */
  /*  EPS : 特異性の判定値を与える. ふつうは 1.0×10^(-14)              */

  k=0;
  for(j=1;j<=M;j++)
    for(i=1;i<=M;i++)
      if(A(i,j) != 0.0){
	A[k] = A(i,j); 
	IA[k] = i;
	k++;
      }

  B   = &b[1]    ;

  EPS = 1.0e-14  ;

  /* (a) リンクの方法                                                   */
  /* CALL ETALU(A,IA,L,L2,B,NW,M,EPS,NUM,X,E,IE,WE,NETA,IP,IS,IER)      */

  etalu_(A,IA,&L,&L2,B,NW,&M,&EPS,&NUM,X,E,IE,WE,NETA,IP,&IS,IER);

  for(i=0; i<M; i++) B[i] = X[i+1];
  return IER;
}

