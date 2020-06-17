/* chainenergy.c */
/* a function to calculate the total energy and force in the chaine */
/* edited 6/5/2013 by MAK */



double lj(){

  int ii,jj,kk;
  double Rsq,flj,ulj,R1,R2,S1,S2,Ssq;
  double Ulenna,Flenna;
  double elj=0;
  double epsilon;
  double sigma;
  double diff[DIM]; 
  double ushift=0;


  for (jj = 0; jj < N; jj++){
    for (kk=0 ; kk<DIM ; kk++){
      Flj[jj][kk] = 0;
    }
  }

  for( ii = 0 ; ii < N-1; ii++){
    for (jj =ii+1 ; jj < N; jj++){
      Rsq=0;
      for (kk=0 ; kk<DIM ; kk++){
	diff[kk] = pos[ii][kk]-pos[jj][kk]; 
	Rsq=Rsq + diff[kk]*diff[kk];
      }
           
      sigma = 0.5*(sig[ii]+sig[jj]);
      epsilon = sqrt(eps[ii]*eps[jj]);
      R1 = (sigma)*(sigma)/Rsq; 
      R2 = R1*R1*R1;
      S1 = (sigma)*(sigma)/(rcut*rcut);
      S2 = S1*S1*S1;
      ushift = 4*epsilon*((S2*S2)-S2);
      ulj = 4*epsilon*((R2*R2)-R2);
      Ssq = rcut * rcut ;
      if ( Rsq < Ssq ){
        Ulenna = ulj - ushift;}
      else if (Rsq > Ssq){
	Ulenna = 0 ;}
     
      elj= elj+Ulenna;
     
      for (kk=0 ; kk<DIM ; kk++){
	flj = (48*epsilon/(sigma*sigma))
	  *((R2*R2*R1)-(0.5*R2*R1))*(pos[ii][kk]-pos[jj][kk]);
        if ( Rsq < Ssq ){
          Flenna = flj;}
        else if (Rsq > Ssq){
          Flenna = 0 ;}
   
	Flj[ii][kk]=Flj[ii][kk]+ Flenna;
	Flj[jj][kk]=Flj[jj][kk]- Flenna;

      }
    }
  }
  
  return elj;
}

/* -------------------------- end of function lj  --------  */


double spring(){

  int ii,jj,kk;
  double diff[DIM];
  double dist = 0,Rsq = 0;
  double ucont =0, fcont=0;
  double esp = 0;
 
   
  for (jj = 0; jj < N; jj++){
    for (kk=0 ; kk<DIM ; kk++){
      Fsp[jj][kk] = 0;
    }
  }
     
  for (ii = 0 ; ii < N -1; ii++){
    Rsq = 0;
    for (kk=0 ; kk<DIM ; kk++){
      diff[kk] = pos[ii][kk]-pos[ii+1][kk]; 
      Rsq=Rsq + diff[kk]*diff[kk];
    }
   
    dist=sqrt(Rsq);
    ucont = 0.5*sc*(dist-sb)*(dist-sb);
    esp= esp + ucont;

    for (kk=0 ; kk<DIM ; kk++){
      fcont = -sc*(dist-sb)*(pos[ii][kk]-pos[ii+1][kk])/dist;
      Fsp[ii][kk]=Fsp[ii][kk]+ fcont;
      Fsp[ii+1][kk]=Fsp[ii+1][kk]- fcont;}
  }
   

  return esp;
}

/* -------------------------- end of function spring --------------- */ 

double external(){

  int ii,jj,kk;
  double ext = 0;
 
  for (jj = 0; jj < N; jj++){
    for (kk=0 ; kk<DIM ; kk++){
      Fext[jj][kk] = 0;
    }
  }
    

  Fext[0][0]=-F0;
  Fext[0][1]=0;
  Fext[0][2]=0;
    
  Fext[N-1][0]=F0;
  Fext[N-1][1]=0;
  Fext[N-1][2]=0;

  for(ii =1 ; ii<N-1 ; ii++)
    {
      Fext[ii][0]=0;
      Fext[ii][1]=0;
      Fext[ii][2]=0;
    }
    
  return ext;
}


/* -------------------------- end of function external -------------------- */ 

double force(){

  int ii,jj,kk;
  double epot =0;
  double elj =0 ;
  double esp =0 ;
  double ext =0;
  double pesp =0 ;
  double emo =0;
  
  for (jj = 0; jj < N; jj++){
    for (kk=0 ; kk<DIM ; kk++){
      Fint[jj][kk] = 0;
    }
  }

 for (jj = 0; jj < M; jj++){
    for (kk=0 ; kk<DIM ; kk++){
      Fintp[jj][kk] = 0;
    }
  }
  
  elj = lj();
  esp = spring();
  ext = external();
  pesp = pspring();
  emo = morse();
    
  for (ii = 0 ; ii< N ; ii++){
    for (jj =0 ; jj<DIM ; jj++) {
     
      Fint[ii][jj] =Flj[ii][jj]+ Fsp[ii][jj]+Fext[ii][jj]+Fmoc[ii][jj];

    }}

  for (ii = 0 ; ii< M ; ii++){
    for (jj =0 ; jj<DIM ; jj++) {
     
      Fintp[ii][jj] =Fspp[ii][jj]+Fmop[ii][jj];

    }}

  epot = esp+elj+pesp+emo;

  return epot;
}

/* -------------------------- end of function force -------- */ 
