{ int *p,i,j,k,m,n,w,pi,pilimit; 
  cp(A,a);cp(X,x);n=dim2(a);w=halfbw(a);m=dim2(x);

  if(NULL==(p=LU(a))) return 0;

  forall(1,k,m){
    b = x[k];
    forall(1,i,n){
      pilimit=less(i+w,n); bpi=b[p[i]];
      forall(i+1,j,pilimit) b[j] += a[j][i]*bpi; 
    }
    for(i=n;i>=1;i--){
      s = b[(pi=p[i])]; api=a[pi]; pilimit=less(pi+w,n);
      forall(i+1,j,pilimit) s -= api[j]*b[j];
      b[pi] = s/api[i];
    }
  }
  return 1;
}
