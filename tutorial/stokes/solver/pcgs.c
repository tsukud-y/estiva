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

static int count_NL(void *pA)
{
  long i, j, k, N, NL ;

  N = dim1(pA);

  for(NL=0, i=1; i<=N; i++){
    for(k=0, j=1; j<=N; j++) if(A(i,j) != 0.0) k++;
    if(NL<k) NL = k;
  }
  return NL;
}


int pcgs(void* pA, double* b)
{

  /* (i)  引数の型と種類                                                */
  D_     D, B         ;/*    D  一次元配列 D(N), B(N)                   */
  D__    A            ;/*    D  二次元配列 A(N1, 2 * NL)                */
  I__    IA           ;/*    I  二次元配列 IA(N1, 2 * NL)               */
  D_     R            ;/*    D  一次元配列 R(N)                         */
  long   NL, N1, N, ITR, IER                                            ;
  double EPS, S                                                         ;
  D_     X, DD, P, Q  ;/*    D  一次元配列で, 要素は 0〜N               */
  I_     M            ;/*    I  一次元配列  M(2 * N)  作業用            */

  /*--    DからRまでと, EPSからSまでの引数は, サブルーチンPCGと同じ   --*/  

  D_     R0, E, H     ;/*    D  一次元配列  要素数はN                   */
  D_     W            ;/*    D  一次元配列  W(0:N)                      */
  
  long   i, j, k      ;
  
  N   = dim1(pA)      ;
  NL  = count_NL(pA)  ;
  
  ary1(  D, N   )         ;
  ary2(  A, 2*NL, N+2*NL );
  ary2( IA, 2*NL, N+2*NL );
  ary1(  R, N   )         ;
  ary1(  X, N+1 )         ;  ary1( DD, N+1 );  ary1( P, N+1 );  ary1( Q, N+1 ); 
  ary1(  M, 2*N )         ;

  ary1( R0, N   )         ;  ary1(  E, N   );  ary1( H, N   );
  ary1(  W, N+1 )         ;
  

  /* (ii) 主プログラム → サブルーチン                                  */
  /*      サブルーチンをCALLするときには, つぎの値を与える.             */
  /*  D   : 配列Dの第1〜第n位置に行列Aの対角要素を入れておく            */
  /*  N   : 行列Aの行数を入れておく.                                    */
  /*  N1  : 配列Aの行数を入れておく. N1≧N+2*NL でないといけない.       */
  /*  NL  : 行列Aの各行における非ゼロ要素数の最大値を入れておく.        */
  /*  B   : 連立一次方程式の右辺を入れておく.                           */
  /* EPS  : 収束判定置を入れておく. ふつうは 1.×10^(-7)                */
  /* ITR  : 打切りまでの最大繰返し回数を入れておく.                     */

  /*--    D, N, N1, NL, B, EPS, ITR については, サブルーチンPCGと同じ --*/

  /*  A   : 配列Aの各行1〜NL要素は, 行列Aの下三角部分の各行の非ゼロ     */
  /*        要素を入れる. また, 各行のNL+1〜2*NL要素は, 上三角部分      */
  /*        の各行の非ゼロ要素を入れる. ただし, 各行の対角要素は配列Dに */
  /*        入れる.                                                     */
  /* IA   : 配列Aに入れた要素の列番号を, 対応する位置に入れておく.      */
  /*  S   : LUCGS法のとき 0., MLUCG法のときσ(>0)を入れておく.          */  

  for(i=1; i<=N; i++) D[i-1] = A(i,i); 

  for(i=1; i<=N; i++)    for(k=0, j=1; j<i; j++)    if(A(i,j) != 0.0){
    A[k][i-1]  = A(i,j)  ;
    IA[k][i-1] = j       ;
    k++                  ;
  }

  for(i=1; i<=N; i++)    for(k=NL, j=i+1; j<=N; j++) if(A(i,j) != 0.0){
    A[k][i-1]  = A(i,j)  ;
    IA[k][i-1] = j       ;
    k++                  ;
  }
      

  N1  =  N+2*NL;

  B   = &b[1]  ;
  EPS = 1.0e-7 ;
  ITR = N      ;
  S   = 0.     ;
  
  /* (a) リンクの方法                                                   */
  /* CALL PCGS(D,A,IA,N,N1,NL,B,EPS,ITR,S,X,DD,P,Q,R,R0,E,H,W,M,IER)    */

  pcgs_(D,A[0],IA[0],&N,&N1,&NL,B,&EPS,&ITR,&S,X,DD,P,Q,R,R0,E,H,W,M,&IER);


  for(i=0; i<N; i++) B[i] = X[i+1];
  return IER;
}
