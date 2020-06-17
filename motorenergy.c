/* motorenergy.c */
/* a function to calculate the total energy and force in the motor protein */
/* edited 6/14/2013 by MAK */


double pspring(){

  int ii,jj,kk;
  double diff[DIM];
  double dist = 0,Rsq = 0;
  double ucont =0, fcont=0;
  double pesp = 0;
 
   
  for (jj = 0; jj < M; jj++){
    for (kk=0 ; kk<DIM ; kk++){
      Fspp[jj][kk] = 0;
    }
  }
     
  for (ii = 0 ; ii < M -1; ii++){
    Rsq = 0;
    for (kk=0 ; kk<DIM ; kk++){
      diff[kk] = posp[ii][kk]-posp[ii+1][kk]; 
      Rsq=Rsq + diff[kk]*diff[kk];
    }
   
    dist=sqrt(Rsq);
    ucont = 0.5 *psc* (dist-psb)*(dist-psb);
    pesp= pesp + ucont;

    for (kk=0 ; kk<DIM ; kk++){
      fcont = -psc*(dist-psb)*(posp[ii][kk]-posp[ii+1][kk])/dist;
      Fspp[ii][kk]=Fspp[ii][kk]+ fcont;
      Fspp[ii+1][kk]=Fspp[ii+1][kk]- fcont;}
  }
   

  return pesp;

}
