/* file: mandc.c */
/* usage: Motor Protein and chain */
/* edited: June 13,2013 by Maral */

# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <string.h>
# include <time.h>

# include "parameters.h"

/* single Particle 3 dimensional velocity and position */

double pos[N][DIM];
double posp[M][DIM];/* motor*/
double posini[N][DIM]; /* initial position */
double posinip[M][DIM]; /* initial position */
double vel[N][DIM];
double velp[M][DIM];
double Fint[N][DIM]; /* total force on each chain */
double Fintp[M][DIM]; /* total force on each protein */
double Flj[N][DIM];  /* Lennard Jones force on each prticle */
double Fsp[N][DIM];   /* Spring force on each prticle */
double Fspp[M][DIM];   /* Spring force on each motor */
double Fext[N][DIM];  /* External force on each prticle */
double Fmoc[N][DIM];  /* morse force on each chain */
double Fmop[M][DIM];  /* Morse force on each protein */
double acc[N][DIM];
double accp[M][DIM];
double mass[N];
double massp[M];
double p[N][DIM];
double pp[M][DIM];
double ptotal[DIM];
double ptraj[NTRAJ][DIM];
double ktraj[NTRAJ],pottraj[NTRAJ];
double ttraj[NTRAJ];
double sig[N];
double eps[N];
double length[DIM];

/* the function that calculate the tolal momentum and kinetic energy */
# include "kineticenergy.c"  
/* the function that calculate the tolal energy of protein by pesp and Fspp */
# include "motorenergy.c"  
/* the function that calculate the  energy of protein and chain by emo and Fmo */ 
# include "chainandmotor.c"
/* the function that calculate the tolal energy of chain by epot and Fint */  
# include "chainenergy.c" 



FILE *g_out;

int main()
{
  int ii,jj,kk,nn,ntt;
  double epot=0,ktot=0,emo=0;
  double ediff,ediffrel;
  double mtot=0,pmtot=0;
  double a ;
  double deltax = 0.01 ;
  int nt;
  double Etotini = 0; /* initial total energy */
  double dN; /* particle number as double precision variable */

  /* initial value */

  printf(" #index  chain  mass\n");
  for (ii=0 ; ii < N ; ii++){
    mass[ii]= chm;
    sig[ii] =chsig ;
    eps[ii] = cheps;
    printf ( " %3d    %6.2f \n",ii,mass[ii]);
  }
  printf(" #index  protein  mass\n");
  for (ii=0 ; ii < M ; ii++){
    massp[ii]= pm;
    printf ( " %3d    %6.2f \n",ii,massp[ii]);
  }   
 
  

 
  /* calculating the position of the chain particles */
   
  dN = (double)N;
  a = sb + deltax/(dN-1);
  for ( nn= 0 ; nn < N ; nn++){
    for (ii = 0; ii < Nx ; ii++){
      for (jj = 0; jj < Ny ; jj++){
	for ( kk = 0; kk < Nz ; kk++) {
	  pos[nn][0]=ii*a;
	  pos[nn][1]=jj*a;
	  pos[nn][2]=kk*a;
	  nn = nn+1;
	}
      }
    }
  }


  /* calculating the position of the motor particles */
  for ( nn= 0 ; nn < M ; nn++){
    for (ii = 0; ii < Mx ; ii++){
      for (jj = 0; jj < My ; jj++){
	for ( kk = 0; kk < Mz ; kk++) {
	  posp[nn][0]=a*(N-1)/2 - psb/2+ ii*psb;
	  posp[nn][1]=(jj+1)*h;
	  posp[nn][2]=0;
	  nn = nn+1;
	}
      }
    }
  }

  printf("#index    xc        yc       zc \n");
  for( ii = 0 ; ii < N; ii++){
    printf (" %3d ",ii); 
    for (kk=0 ; kk<DIM ; kk++){
      printf (" %6.2f  ",pos[ii][kk]); 
    }
    printf ("\n");
  }

  printf("#index   px        py      pz \n");
  for( ii = 0 ; ii < M; ii++){
    printf (" %3d ",ii); 
    for (kk=0 ; kk<DIM ; kk++){
      printf (" %6.2f  ",posp[ii][kk]); 
    }
    printf ("\n");
  }

  g_out = fopen("inipos.txt","w");
 
  fprintf(g_out,"#index   chx        chy      chz \n");
  for( ii = 0 ; ii < N; ii++){
    fprintf (g_out," %3d ",ii); 
    for (kk=0 ; kk<DIM ; kk++){
      fprintf (g_out," %6.2f  ",pos[ii][kk]); 
    }
    fprintf (g_out,"\n");
  }
  fclose(g_out);
 
  


  /* calculating the velocity of the particles */
 
  for (jj = 0; jj < N ; jj++){
    mtot=mtot + mass[jj];
    for ( kk = 0; kk < DIM ; kk++) {
      vel[jj][kk]= 0 ;

    }     
  }
  printf("#index   chvx        chvy      chvz" );
  printf("\n");
  for( ii = 0 ; ii < N; ii++){
    printf (" %3d ",ii); 
    for (kk=0 ; kk<DIM ; kk++){
      printf (" %6.2f  ",vel[ii][kk]);
    }
    printf ("\n");
  }

  for (jj = 0; jj < M ; jj++){
    pmtot=pmtot + massp[jj];
    for ( kk = 0; kk < DIM ; kk++) {
      velp[jj][kk]= 0 ;

    }     
  }

  printf("#index    pvx      pvy      pvz" );
  printf("\n");
  for( ii = 0 ; ii < M; ii++){
    printf (" %3d ",ii); 
    for (kk=0 ; kk<DIM ; kk++){
      printf (" %6.2f  ",velp[ii][kk]);
    }
    printf ("\n");
  }

  
  /* calculating the force of particles */
  

  epot = force();
  emo = morse();

  printf("#index    chFx        chFy        chFz \n");
  for( ii = 0 ; ii < N; ii++){
    printf (" %3d ",ii); 
    for (kk=0 ; kk<DIM ; kk++){
      printf (" %8.4f  ",Fint[ii][kk]); 
    }
      
    printf ("\n");
  }
  printf("#index    pFx        pFy        pFz \n");
  for( ii = 0 ; ii < M; ii++){
    printf (" %3d ",ii); 
    for (kk=0 ; kk<DIM ; kk++){
      printf (" %8.4f  ",Fintp[ii][kk]); 
    }
      
    printf ("\n");
  }

  /* calculating the accelaration of the particles */
  printf("#index   CHax       CHay        CHaz \n");
  for( ii = 0 ; ii < N; ii++){
    printf (" %3d ",ii); 
    for (kk=0 ; kk<DIM ; kk++){
      acc[ii][kk]=Fint[ii][kk]/mass[ii];
      printf (" %8.4f  ",acc[ii][kk]); 
    }
      
    printf ("\n");
  }
  printf("#index    Pax        Pay        Paz \n");
  for( ii = 0 ; ii < M; ii++){
    printf (" %3d ",ii); 
    for (kk=0 ; kk<DIM ; kk++){
      accp[ii][kk]=Fintp[ii][kk]/massp[ii];
      printf (" %8.4f  ",accp[ii][kk]); 
    }
      
    printf ("\n");
  }

  /* finished setting initial conditions */
  /* write to trajectory */
  epot = force();
  ktot = velocity();

  for(kk=0 ; kk<DIM; kk++){
    ptraj[0][kk] = ptotal[kk];
  }
  ktraj[0] = ktot;
  pottraj[0] = epot;

 
  Etotini = pottraj[0]+ktraj[0]; /* initial total energy */

  ntt = 0;
  for ( nt=1; nt<NSTEP ;nt++){  /* loop over time steps */
   
    epot = force();
     
    for (jj = 0; jj < N ; jj++){
      for ( kk = 0; kk < DIM ; kk++){
	pos[jj][kk] = pos[jj][kk] + vel[jj][kk]*delt+ 
	  0.5*delt*delt*Fint[jj][kk]/mass[jj];             
	vel[jj][kk] = vel[jj][kk] + 0.5*delt*Fint[jj][kk]/mass[jj];
      }
    }  
       
     for (jj = 0; jj < M ; jj++){
    for ( kk = 0; kk < DIM ; kk++){
    posp[jj][kk] = posp[jj][kk] + velp[jj][kk]*delt+ 
      0.5*delt*delt*Fintp[jj][kk]/massp[jj];             
    velp[jj][kk] = velp[jj][kk] + 0.5*delt*Fintp[jj][kk]/massp[jj];
    }
    }
    epot = force();
    for (jj = 0; jj < N ; jj++){
      for ( kk = 0; kk < DIM ; kk++) {
	vel[jj][kk] = vel[jj][kk] + 0.5*delt*Fint[jj][kk]/mass[jj];
      }
    }
     for (jj = 0; jj < M ; jj++){
    for ( kk = 0; kk < DIM ; kk++) {
    	velp[jj][kk] = velp[jj][kk] + 0.5*delt*Fintp[jj][kk]/massp[jj];
     }
    }
     
      epot = force();
      ktot = velocity();
    if (nt%NGAP == 0){
     
      ntt = ntt + 1;
      ttraj[ntt] = ttraj[ntt-1] + NGAP*delt ;
      pottraj[ntt] = epot;
      ktraj[ntt]= ktot;
      for (kk=0 ; kk<DIM ; kk++){
	ptraj[ntt][kk]= ptotal[kk];
      }
    }

  }/* the end of time steps */

  g_out = fopen("pandc.dat","w");
  fprintf (g_out,"#time  potential-energy kinetik-energy total-energy  ptotx ptoty ptotz\n");
  for( nt=0 ; nt<NTRAJ ; nt++){
    fprintf(g_out," %8.4e   %8.4e  %8.4e %8.4e  %8.4e  %8.4e   %8.4e  \n",
	    ttraj[nt],pottraj[nt],ktraj[nt],pottraj[nt]+ktraj[nt],
	     ptraj[nt][0],ptraj[nt][1],ptraj[nt][2]);
  }
  fclose(g_out);
  ediff = pottraj[NTRAJ-1]+ktraj[NTRAJ-1]-Etotini;
  ediffrel = 2.0*ediff/(pottraj[NTRAJ-1]+ktraj[NTRAJ-1]+Etotini);
 
  printf(" Efinal - Einitial = %12.4e relative deviation %12.4e \n",
	 ediff,ediffrel);     

  return 0;
}

/* -------------------------- end of main program  --------  */


