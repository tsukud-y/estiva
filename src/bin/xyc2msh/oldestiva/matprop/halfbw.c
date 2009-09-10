{ int n,i,j;
  n = dim2(a);
  forall(1,i,n)forall(1,j,i)
    if(a[j][n-(i-j)] != 0.0||a[n-(i-j)][j] != 0.0) return n-i;
}
