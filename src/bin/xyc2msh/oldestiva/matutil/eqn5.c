{ int i,j,e,m; 
  args(argv,&P,&HXtAX,&Rx,&HYtAY,&Ry);

  e=dim1(P);m=dim1(Rx);forall(1,i,e){
    Pi = 0.0; HXtAXi=HXtAX[i]; HYtAYi=HYtAY[i];
    forall(1,j,m) Pi += HXtAXi[j]*Rx[j] + HYtAYi[j]*Ry[j];
    P[i] = Pi;
  }
}
