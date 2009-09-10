{ static int i,j,k,m,n,**Q,*Qi,Qi0;
  args(argv,&T,&A,&B);

  m=dim2(A);n=dim1(A);ary2(Q,n+1,3);forall(1,j,n){
    k=1;
    forall(1,i,m)if(A[i][j] != 0.0) Q[j][k++] = i; 
    Q[j][0] = k-1;
  }
  m=dim2(T);n=dim1(T);forall(1,i,m){
    Ti=T[i]; Qi=Q[i]; Qi0=Qi[0]; 
    forall(1,j,n){
      Tij=0.0;
      forall(1,k,Qi0) Tij += A[Qi[k]][i]*B[Qi[k]][j];
      Ti[j]=Tij;
    }
  }
}
