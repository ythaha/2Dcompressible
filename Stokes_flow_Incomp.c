/*   Functions which solve for the velocity and pressure fields using Uzawa-type iteration loop.  */

#include <math.h>
#include <malloc.h>
#include <sys/types.h>
#include "element_definitions.h"
#include "global_defs.h"

extern int Emergency_stop;

/* Master loop for pressure and (hence) velocity field */

void solve_constrained_flow_iterative(E)
     struct All_variables *E;

{ 
    double *D1;
    double *u;
    double *R,*Bp;
    double residual_ddash;
    double vmag;
    double vdot(),pdot();

    static int been_here = 0;
   
    float solve_Ahat_p_fhat();
    float solve_Ahat_p_fhat_compressible();
    void v_from_vector();
    void dp_to_nodes();   
    void thermal_buoyancy();
    void assemble_forces();
    void lagr_multi_compressible();
   
    int steps,cycles;
    int i,j,k,doff,vel_cycles_previous,vel_calls_previous;
  
    double time,CPU_time0();
 
    const int npno = E->mesh.npno;
    const int gnpno = E->mesh.npno;
    const int nno = E->mesh.nno;
    const int dims = E->mesh.nsd;
    const int neq = E->mesh.neq;
    const int gneq = E->mesh.neq;

    double *r_compressible;
    double relative_err_v, relative_err_p;
    double *old_v, *old_p;
    double *v2_temp,*grad_v2;
    double temp1,temp2,temp3,temp4;
    int num_of_loop;

    old_v = (double *)malloc((neq+1)*sizeof(double));
    old_p = (double *)malloc((npno+1)*sizeof(double));

    //v2_temp and grad_v2 are used to deal with Lagrangian multiplier
    v2_temp = (double *)malloc((npno+1)*sizeof(double));
    grad_v2 = (double *)malloc((neq+1)*sizeof(double));
    for(i=0;i<neq;i++)  {
	    old_v[i] = 0.0;
	    grad_v2[i] = 0.0;
    }
    for(i=1;i<=npno;i++)  {
	    old_p[i] = 0.0;
	    v2_temp[i] = 0.0;
    }
    

    time=CPU_time0();
    
    cycles=E->control.p_iterations;
   
    /* Solve for velocity and pressure, correct for bc's */

    if(E->control.compressibility == 0) {
    	residual_ddash=solve_Ahat_p_fhat(E,E->U,E->P,E->F,E->control.accuracy,&cycles);
    }
    else {
      relative_err_v = 1.0;
      relative_err_p = 1.0;
      num_of_loop = 0;

      //we don't use guessed U and P
//      for(i=0;i<neq;i++) {
//	  E->U[i]=0.0;
//      }
//      for(i=1;i<=npno;i++) {
//	  E->P[i]=0.0;
//      }


      //add incompressibility constaints can improve the accuracy
      //like the imbalance between two heating terms may be several percent
      //at the beginning if you don't use tole_comp as a criteria 
      //while((relative_err_v >= E->control.relative_err_accuracy || relative_err_p >= E->control.relative_err_accuracy || E->monitor.incompressibility >= E->control.tole_comp) 
      while((relative_err_v >= E->control.relative_err_accuracy || relative_err_p >= E->control.relative_err_accuracy) 
		      && num_of_loop <= E->control.compress_iter_maxstep) {

	if(E->control.dyn_p_in_buoyancy) {
             thermal_buoyancy(E);
             assemble_forces(E,0);
        }

	//deal with the Lagrangian multiplier with compressibility
	if (E->control.augmented_Lagr) {
	    lagr_multi_compressible(E,v2_temp,grad_v2,E->F);
	}

	for(i=0;i<neq;i++) {
		old_v[i] = E->U[i];
	}
	for(i=1;i<=npno;i++) {
		old_p[i] = E->P[i];
	}
	
    	residual_ddash=solve_Ahat_p_fhat(E,E->U,E->P,E->F,E->control.accuracy,&cycles);

	temp1 = temp2 = 0.0;
        for(i=0;i<neq;i++)  {
             temp1 += (E->U[i]-old_v[i])*(E->U[i]-old_v[i]);
             temp2 += E->U[i]*E->U[i];
        }
        temp1 = sqrt(temp1/neq);
        temp2 = sqrt(temp2/neq);
//        fprintf(stderr,"|| Delta V || = %lf, || V || = %lf \n", temp1, temp2);
        relative_err_v = temp1/temp2;
	temp3 = temp4 = 0.0;
	for(i=1;i<=npno;i++)  {
	     temp3 += (E->P[i]-old_p[i])*(E->P[i]-old_p[i]);
	     temp4 += E->P[i]*E->P[i];
	}
        temp3 = sqrt(temp3/npno);
        temp4 = sqrt(temp4/npno);
//        fprintf(stderr,"|| Delta P || = %lf, || P || = %lf \n", temp3, temp4);
        relative_err_p = temp3/temp4;
	
	num_of_loop++;
	fprintf(stderr, "Relative error err_v / v = %f and err_p / p = %f after %d loops\n", relative_err_v, relative_err_p, num_of_loop);
	fprintf(E->fp, "Relative error err_v / v = %f and err_p / p = %f after %d loops\n\n", relative_err_v, relative_err_p, num_of_loop);
	fflush(E->fp); //using fflush here is very very important, or the program may not run as you want

      } //end for while
  
    } //end for else

    been_here++;

    v_from_vector(E,E->V,E->U); 
    dp_to_nodes(E,E->P,E->NP,E->mesh.levmax);   

    free((void *) old_v);
    free((void *) old_p);

    free((void *) v2_temp);
    free((void *) grad_v2);
    
  return; 
}
/*-------------------------------------------------*/
void solve_constrained_flow_iterative_anotherway(E)
     struct All_variables *E;

{ 
    double *D1;
    double *u;
    double *R,*Bp;
    double residual_ddash;
    double vmag;
    double vdot(),pdot();

    static int been_here = 0;
   
    float solve_Ahat_p_fhat();
    float solve_Ahat_p_fhat_compressible();
    void v_from_vector();
    void dp_to_nodes();   
    void thermal_buoyancy();
    void assemble_forces();
   
    int steps,cycles;
    int i,j,k,doff,vel_cycles_previous,vel_calls_previous;
  
    double time,CPU_time0();
 
    const int npno = E->mesh.npno;
    const int gnpno = E->mesh.npno;
    const int nno = E->mesh.nno;
    const int dims = E->mesh.nsd;
    const int neq = E->mesh.neq;
    const int gneq = E->mesh.neq;

    double *r_compressible;
    double relative_err_v, relative_err_p;
    double *old_v, *old_p;
    double temp1,temp2,temp3,temp4;
    int num_of_loop_out,num_of_loop_in;

  old_v = (double *)malloc((neq+1)*sizeof(double));
  old_p = (double *)malloc((npno+1)*sizeof(double));
  for(i=0;i<neq;i++)  {
	    old_v[i] = 0.0;
  }
  for(i=1;i<=npno;i++)  {
	    old_p[i] = 0.0;
  }
    

  time=CPU_time0();
    
  cycles=E->control.p_iterations;
   
    /* Solve for velocity and pressure, correct for bc's */

  if(E->control.compressibility == 0) {
    	residual_ddash=solve_Ahat_p_fhat(E,E->U,E->P,E->F,E->control.accuracy,&cycles);
  }
  else {
    relative_err_p = 1.0;
    num_of_loop_out = 0;
    //outside loop for P
    while(relative_err_p >= E->control.relative_err_accuracy
		                            && num_of_loop_out <= E->control.compress_iter_maxstep) {

      for(i=1;i<=npno;i++) {
                old_p[i] = E->P[i];
      }
      if(E->control.dyn_p_in_buoyancy) {
	      thermal_buoyancy(E);
	      assemble_forces(E,0);
      }

      relative_err_v = 1.0;
      num_of_loop_in = 0;
      //inner loop for V
      while(relative_err_v >= E->control.relative_err_accuracy 
		      && num_of_loop_in <= E->control.compress_iter_maxstep) {

	for(i=0;i<neq;i++) {
		old_v[i] = E->U[i];
	}
	
    	residual_ddash=solve_Ahat_p_fhat(E,E->U,E->P,E->F,E->control.accuracy,&cycles);

	temp1 = temp2 = 0.0;
        for(i=0;i<neq;i++)  {
             temp1 += (E->U[i]-old_v[i])*(E->U[i]-old_v[i]);
             temp2 += E->U[i]*E->U[i];
        }
        temp1 = sqrt(temp1/neq);
        temp2 = sqrt(temp2/neq);
        relative_err_v = temp1/temp2;
	num_of_loop_in++;
	fprintf(stderr, "Relative error err_v / v = %f after %d inner loops\n", relative_err_v, num_of_loop_in);
      } //end for inner loop while
      
      temp3 = temp4 = 0.0;
      for(i=1;i<=npno;i++)  {
	     temp3 += (E->P[i]-old_p[i])*(E->P[i]-old_p[i]);
	     temp4 += E->P[i]*E->P[i];
      }
      temp3 = sqrt(temp3/npno);
      temp4 = sqrt(temp4/npno);
      relative_err_p = temp3/temp4;
      fprintf(stderr,"|| Delta P || = %lf, || P || = %lf \n", temp3, temp4);
      num_of_loop_out++;

      fprintf(stderr, "Relative error err_p / p = %f after %d loops\n\n", relative_err_p, num_of_loop_out);
      //using fflush for write to file is very very important, or the program may not run as you want

      } //end for outside loop while
  
    } //end for else

    been_here++;

    v_from_vector(E,E->V,E->U); 
    dp_to_nodes(E,E->P,E->NP,E->mesh.levmax);   

    free((void *) old_v);
    free((void *) old_p);
    
  return; 
}

/*  ==========================================================================  */

float solve_Ahat_p_fhat(E,V,P,F,imp,steps_max)

     struct All_variables *E;
     double *V,*P,*F;
     double imp;
     int *steps_max;
    
{ 
  int i,j,k,ii,count,convergent,valid,problems,lev,lev_low,npno,neq,steps;
  int gnpno,gneq;
  
  static int been_here = 0;
  double *p1,*r1,*u;
  double *r0,*r2,*z0,*z1,*s1,*s2,*Ah,*u1;
  double *shuffle, *R;
  double alpha,delta,s2dotAhat,r0dotr0,r1dotz1;
  double residual, initial_residual, last_residual,v_res;
  double dpressure,dvelocity;
  
  double vdot(),pdot();

  float CPU_time();
  double time0,time,CPU_time0();

  void assemble_div_u();
  void assemble_del2_u();
  void assemble_grad_p();
  void strip_bcs_from_residual();
  int  solve_del2_u();

  double *r_compressible;
   
  const int dims=E->mesh.nsd;
  const int n=loc_mat_size[E->mesh.nsd];

  npno=E->mesh.npno;
  neq=E->mesh.neq;

  gnpno=E->mesh.npno;
  gneq=E->mesh.neq;

  r0 = (double *)malloc((npno+1)*sizeof(double));
  r1 = (double *)malloc((npno+1)*sizeof(double));
  r2 = (double *)malloc((npno+1)*sizeof(double));
  z0 = (double *)malloc((npno+1)*sizeof(double));
  z1 = (double *)malloc((npno+1)*sizeof(double));
  s1 = (double *)malloc((npno+1)*sizeof(double));
  s2 = (double *)malloc((npno+1)*sizeof(double));
  p1 = (double *)malloc((npno+1)*sizeof(double));
  Ah = (double *)malloc((neq+1)*sizeof(double));
  u1 = (double *)malloc((neq+1)*sizeof(double));

  r_compressible = (double *)malloc((npno+1)*sizeof(double));

  problems=0;
  time0=time=CPU_time0();

  been_here ++;
 

  /* calculate the velocity residual, note there are tricks involved here */

  lev=E->mesh.levmax;

  /* The assemble_C_u must be put here. If move this after the modification of V[i] (V[i] += u1[i]), 
   * it will affect the convergence of Conjugate gradient. The reason? still don't know */
  if(E->control.compressibility != 0.0) {
  	assemble_C_u(E,V,r_compressible,lev);
  }

  for(i=0;i<neq;i++) 
	V[i] = 0.0;
  for(i=1;i<=npno;i++) 
	P[i] =0.0;

  assemble_grad_p(E,P,Ah,lev);
  assemble_del2_u(E,V,u1,lev,1);

  for(i=0;i<neq;i++) 
      Ah[i] = F[i] - Ah[i] - u1[i]; 

  v_res=sqrt(vdot(E,F,F,lev)/gneq);

  fprintf(stderr,"initial residue of momentum equation F %.8e %d\n",v_res,gneq);

  strip_bcs_from_residual(E,Ah,lev);

  valid=solve_del2_u(E,u1,Ah,imp*v_res,E->mesh.levmax);
  strip_bcs_from_residual(E,u1,lev);

  if(!valid) problems++; 
  
  for(i=0;i<neq;i++)  {
      V[i] += u1[i];
     }

  assemble_div_u(E,V,r1,lev);
//for(i=1;i<=npno;i++) {
//        printf("%d %g\n",i,r1[i]);
//  }

  if(E->control.compressibility != 0.0) {
    for(i=1;i<=npno;i++) {
          r1[i] += r_compressible[i];
    }
  }

  residual = initial_residual = sqrt(pdot(E,r1,r1,lev)/gnpno);

  E->monitor.vdotv = sqrt(vdot(E,V,V,lev)/gneq);

  E->monitor.incompressibility = residual/E->monitor.vdotv;
         
   for(i=1;i<=npno;i++)
        p1[i] = 0.0;
     
   count = 0;
   convergent=0;

   if (E->control.print_convergence)  {
         fprintf(E->fp,"AhatP (%03d) after %g seconds with div/v=%.3e for step %d\n",count,CPU_time0()-time0,E->monitor.incompressibility,E->monitor.solution_cycles); /**/
         fflush(E->fp);
         }         

   dpressure = 1.0;
   dvelocity = 1.0;

   //while( count==0 || ((count < *steps_max) && (dvelocity >= imp || dpressure >=imp || E->monitor.incompressibility >= E->control.tole_comp) ))  { 
   while( count==0 || ((count < *steps_max) && (dvelocity >= imp || dpressure >=imp) ) )  {

     for(j=1;j<=npno;j++)
       z1[j] = E->BPI[lev][j]*r1[j];
     
     r1dotz1 = pdot(E,r1,z1,lev);

     if ((count == 0))
       for(j=1;j<=npno;j++)
            s2[j] = z1[j];
     else {
       r0dotr0=pdot(E,r0,z0,lev);
       assert(r0dotr0 != 0.0  /* Division by zero in head of incompressibility iteration */);
       delta = r1dotz1/r0dotr0;
       for(j=1;j<=npno;j++)
            s2[j] = z1[j] + delta * s1[j];
       }
      
     assemble_grad_p(E,s2,Ah,lev); 

     valid=solve_del2_u(E,u1,Ah,imp*v_res,lev);  
     strip_bcs_from_residual(E,u1,lev);

     if(!valid) problems++;
      
     assemble_div_u(E,u1,Ah,lev);

     s2dotAhat=pdot(E,s2,Ah,lev);
      
	                 /* alpha defined this way is the same as R&W */
     alpha = r1dotz1/s2dotAhat; 
     
     for(j=1;j<=npno;j++)   {
       r2[j] = r1[j] - alpha * Ah[j];
//       p1[j] += alpha * s2[j];
       P[j] += alpha * s2[j];       
       }
     
     for(j=0;j<neq;j++)
       V[j] -= alpha * u1[j];
      

     //compute the div/v for compressible or incompressible
     assemble_div_u(E,V,Ah,lev);
     //////
     if(E->control.compressibility != 0.0) {
  	assemble_C_u(E,V,r_compressible,lev);
        for(j=1;j<=npno;j++) {
          Ah[j] += r_compressible[j];
        }
     }
     //////
     E->monitor.vdotv = vdot(E,V,V,E->mesh.levmax);
     E->monitor.incompressibility = sqrt((gneq/gnpno)*(1.0e-32+pdot(E,Ah,Ah,lev)/(1.0e-32+E->monitor.vdotv)));

     dpressure=alpha*sqrt(pdot(E,s2,s2,lev)/(1.0e-32+pdot(E,P,P,lev)));
     dvelocity=alpha*sqrt(vdot(E,u1,u1,lev)/(1.0e-32+E->monitor.vdotv));
         
     count++;
     if (E->control.print_convergence )  {
       fprintf(E->fp,"AhatP (%03d) after %g seconds with div/v=%.7e for step %d dv=%g dp=%g\n",count,CPU_time0()-time0,E->monitor.incompressibility,E->monitor.solution_cycles,dvelocity,dpressure); /**/
       fflush(E->fp);
       }         

     shuffle=s1;s1=s2;s2=shuffle;
     shuffle=r0;r0=r1;r1=r2;r2=shuffle;
     shuffle=z0;z0=z1;z1=shuffle;

     }       /* end loop for conjugate gradient   */

    if(problems) {
      fprintf(E->fp,"Convergence of velocity solver may affect continuity\n");
      fprintf(E->fp,"Consider running with the `see_convergence=on' option\n");
      fprintf(E->fp,"To evaluate the performance of the current relaxation parameters\n");
      fflush(E->fp);
      }

//  for(j=1;j<=npno;j++) 
//      P[j] += p1[j];

    
  free((void *) r0);
  free((void *) r1);       
  free((void *) r2);
  free((void *) z0);
  free((void *) z1);
  free((void *) s1);
  free((void *) s2);
  free((void *) u1);
  free((void *) Ah);
  free((void *) p1);

  free((void *) r_compressible);

  
    *steps_max=count;

    return(residual);
 }


/*  ==========================================================================  */



void v_from_vector(E,V,F)
     struct All_variables *E;
     float **V;
     double *F;
{
  int node,d;
  unsigned int type;

  const int addi_dof = additional_dof[E->mesh.nsd];
  const int nno = E->mesh.nno;
  const int dofs = E->mesh.dof;

  for(node=1;node<=nno;node++)     {
      if(E->node[node] & OFFSIDE) continue;
     
      V[1][node] = F[E->id[node].doff[1]]; 
      V[2][node] = F[E->id[node].doff[2]]; 
      if(dofs==3) V[3][node] = F[E->id[node].doff[3]];
      if (E->node[node] & VBX)
             V[1][node] = E->VB[1][node]; 
      if (E->node[node] & VBZ)
             V[2][node] = E->VB[2][node]; 
      if (dofs==3 && E->node[node] & VBY)
             V[3][node] = E->VB[3][node]; 

    }
  return;
}
