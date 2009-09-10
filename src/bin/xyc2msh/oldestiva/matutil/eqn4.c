{ static int i,j,k,m,n,**Q,*Qj,Qi0;
  args(argv,&B,&HXtAX,&HX,&HYtAY,&HY);

  m=dim2(HX);n=dim1(HX);ary2(Q,n+1,6);forall(1,j,n){
    k=1;
    forall(1,i,m)
      if(HX[i][j] != 0.0||HY[i][j]!=0.0) Q[j][k++] = i;
    Q[j][0] = k-1;
  }
  m=dim2(B);n=dim1(B);forall(1,i,m){
    HXtAXi=HXtAX[i]; HYtAYi=HYtAY[i]; Qi0=Q[i][0]; Bi=B[i]; 
    forall(1,j,n){
      Bij=0.0; Qj=Q[j];
      forall(1,k,Qi0) 
	Bij -= HXtAXi[Qj[k]]*HX[Qj[k]][j] + HYtAYi[Qj[k]]*HY[Qj[k]][j];
      Bi[j]=Bij;
    }
  }
}


