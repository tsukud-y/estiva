{ static int i,j,k,w,n,*p,pk,pklimit,klimit; 
  cp(A,a);n=dim2(a);w=halfbw(a);
  
  if(ary1(p,n+1) != NULL)forall(1,k,n){
    pk=k; nrm=absv(a[k][k]);
    
    forall(k+1,i,less(k+w,n))if(nrm<absv(a[i][k])){
      pk=i; nrm=absv(a[i][k]);
    }
    if(nrm == 0.0) return NULL;
    p[k]=pk;
    klimit=less(k+w,n); apkk= -a[pk][k];
    forall(k+1,i,klimit){
      a[i][k] /= apkk; 
      pklimit=less(pk+w,n); aik=a[i][k]; apk=a[pk]; ai=a[i];
      forall(k+1,j,pklimit) ai[j] += aik*apk[j];
    }
  }
  return p;
}
