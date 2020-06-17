/* chainandmotor.c */
/* a function to calculate the total energy and force between */
/* the motor proteinand chain */
/* edited 6/16/2013 by MAK */


double morse(){

  int ii,jj,kk;
  double Rsq,R1,ucont,fcont,dist;
  double emo=0;
  double diff[DIM]; 
  
  
  bm = psb + log(2)/kappa;

  for (jj = 0; jj < N; jj++){
    for (kk=0 ; kk<DIM ; kk++){
      Fmoc[jj][kk] = 0;
    }
  }
  
  for (jj = 0; jj < M; jj++){
    for (kk=0 ; kk<DIM ; kk++){
       Fmop[jj][kk] = 0;
    }
  }

  for( ii = 0 ; ii < M; ii++){ /* loop over protein sites */
    for (jj =0 ; jj < N; jj++){ /* loop over chain beads */
      Rsq=0;
      for (kk=0 ; kk<DIM ; kk++){
	diff[kk] = posp[ii][kk]-pos[jj][kk]; 
	Rsq=Rsq + diff[kk]*diff[kk];
      }
      
      dist=sqrt(Rsq);
      R1 = 1-exp(-kappa*(dist-bm));
      ucont = (De*R1*R1)-De;
      emo =emo+ucont;


      for (kk=0 ; kk<DIM ; kk++){
	fcont = -2*kappa*De*R1*(exp(-kappa*(dist-bm)))
	  *(posp[ii][kk]-pos[jj][kk])/dist;
	Fmop[ii][kk]=Fmop[ii][kk]+ fcont;
	Fmoc[jj][kk]=Fmoc[jj][kk]- fcont;}
    }
   
  }
  
  return emo;
}
