{ int i,j,m,n;
  args(argv,&U,&A,&R,&H,&P);
      
  m=dim2(H);n=dim1(H);

  ary1(w,m+1);forall(1,j,m){
    wj = R[j]; Hj=H[j];
    forall(1,i,n)if(Hj[i]!=0.0) wj += Hj[i]*P[i];
    w[j] = wj;
  }
  forall(1,i,m){
    Ui = 0.0; Ai=A[i];
    forall(1,j,m) Ui += Ai[j]*w[j];
    U[i] = Ui;
  }
}
