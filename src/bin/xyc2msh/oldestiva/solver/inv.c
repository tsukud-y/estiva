{ int i,j,m; cp(A,a);m=dim2(a);

  ary2(e,m+1,m+1);
  if(e  == NULL) return 0;
  reset(e); forall(1,i,m) e[i][i] = 1.0;
  if(gauss(a,e) == 0) return 0;
  forall(1,i,m)forall(1,j,m) a[i][j] = e[j][i];
  return 1;
}
