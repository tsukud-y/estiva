{ int i,j,m; double t,*tp;
  args(argv,&A,&K,&tp,&M);

  t= *tp;
  m=dim2(A);
  forall(1,i,m){
    Ai=A[i]; Ki=K[i]; Mi=M[i];
    forall(1,j,m) Ai[j] = i==j? Ki[j]+t*Mi: Ki[j];
  }
}
