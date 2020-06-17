/* kineticenergy.c */
/* a function to calculate the kinetic energy and momentum  */
/* edited 6/19/2013 by MAK */

double velocity(){
     
  int kk,jj,ii;
  double ktot=0,kcont=0,kcontp=0;

  for (kk=0 ; kk<DIM ; kk++){
    ptotal[kk] = 0;
  }
  /* calculating the total momentom of the particles */
  for (kk = 0; kk < DIM ; kk++){
    ptotal[kk]=0;
    for ( jj = 0; jj < N ; jj++) {
      p[jj][kk]= mass[jj]*vel[jj][kk];
      ptotal[kk]=ptotal[kk]+p[jj][kk];
    }
    for ( ii = 0; ii < M ; ii++) {
      pp[ii][kk]= massp[ii]*velp[ii][kk];
      ptotal[kk]=ptotal[kk]+pp[ii][kk];
    }
    
  }   
   

  ktot=0;

  for (kk = 0; kk < DIM ; kk++){
    for ( jj = 0; jj < N ; jj++) {
      kcont =0.5*mass[jj]*(vel[jj][kk]*vel[jj][kk]);
      ktot = ktot + kcont;
    }
    for ( ii = 0; ii < M ; ii++) {
      kcontp =0.5*massp[ii]*(velp[ii][kk]*velp[ii][kk]);
      ktot = ktot + kcontp;
    }
          
  }
  
  return ktot;

}

/* -------------------------- end of function velocity  --------  */
