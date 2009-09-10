{ int i,m; double t,*tp;    
  args(argv,&R,&M,&F,&tp,&U);
  t = *tp;
  m=dim1(R);
  forall(1,i,m) R[i] = M[i]*(F[i]+t*U[i]);
}
