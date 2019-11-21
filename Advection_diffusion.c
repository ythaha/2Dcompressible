/*****************************************
 *   CC  III  TTTTT   CC   OO   MM MM    *
 *  C     I     T    C    O  O  M M M    *
 *  C     I     T    C    O  O  M   M    *
 *   CC  III    T     CC   OO   M   M    *
 *                                       *  
 * Developed at CIT for COnvection in    *
 * the Mantle by Louis Moresi 1992-today *
 *                                       *
 * You are free to use this code but it  * 
 * is distrubuted as BeWare i.e. it does *
 * not carry any guarantees or warranties *
 * of reliability.                       *
 *                                       *
 * Please respect all the time and work  *
 * that went into the development of the *
 * code.                                 *  
 *                                       *
 * LM                                    *
 *****************************************/
/*   Functions which solve the heat transport equations using Petrov-Galerkin
     streamline-upwind methods. The process is basically as described in Alex
     Brooks PhD thesis (Caltech) which refers back to Hughes, Liu and Brooks.  */

#include <malloc.h>
#include <sys/types.h>
#include <math.h>
#include "element_definitions.h"
#include "global_defs.h"

extern int Emergency_stop;

struct el { double gpt[9]; };

/* ============================================
   Generic adv-diffusion for temperature field.
   ============================================ */


void advection_diffusion_parameters(E)
     struct All_variables *E;

{

    void std_timestep();
    /* Set intial values, defaults & read parameters*/
    
    E->advection.temp_iterations = 2; /* petrov-galerkin iterations: minimum value. */
    E->advection.total_timesteps = 1; 
    E->advection.sub_iterations = 1;
    E->advection.last_sub_iterations = 1;
    E->advection.gamma = 0.5;
    E->advection.dt_reduced = 1.0;         

    E->monitor.T_maxvaried = 1e05;
 
    input_boolean("ADV",&(E->advection.ADVECTION),"on");
    E->advection.ADVECTION=1;
   
    input_int("visc_heating",&(E->control.visc_heating),"1");
    input_int("adi_heating",&(E->control.adi_heating),"1");
    input_int("latent_heating",&(E->control.latent_heating),"1");

    input_int("minstep",&(E->advection.min_timesteps),"1");
    input_int("maxstep",&(E->advection.max_timesteps),"1000");
    input_int("maxtotstep",&(E->advection.max_total_timesteps),"1000000");
    input_float("finetunedt",&(E->advection.fine_tune_dt),"0.9");
    input_float("fixed_timestep",&(E->advection.fixed_timestep),"0.0");
    input_int("adv_sub_iterations",&(E->advection.temp_iterations),"2,2,nomax");
    input_float("maxadvtime",&(E->advection.max_dimensionless_time),"10.0");
   
    input_float("sub_tolerance",&(E->advection.vel_substep_aggression),"0.005");  
    input_int("maxsub",&(E->advection.max_substeps),"25");
  
    input_float("liddefvel",&(E->advection.lid_defining_velocity),"0.01");
    input_float("sublayerfrac",&(E->advection.sub_layer_sample_level),"0.5");            
    input_int("markers_per_ele",&(E->advection.markers_per_ele),"0");
	      
  /* allocate memory */

  return;
}

void advection_diffusion_allocate_memory(E)
     struct All_variables *E;

{ int i;

  E->Tdot= (float *)malloc((E->mesh.nno+1)*sizeof(float));
  for(i=1;i<=E->mesh.nno;i++) 
    E->Tdot[i]=0.0;

  if (E->control.composition)   {
    E->Cdot= (float *)malloc((E->mesh.nno+1)*sizeof(float));
    for(i=1;i<=E->mesh.nno;i++)
      E->Cdot[i]=0.0;
    if (!(strcmp(E->control.comp_adv_method,"field")==0)) {
      E->advection.markers = E->advection.markers_per_ele*E->mesh.nel;
      for(i=1;i<=E->mesh.nsd;i++)   {
        E->VO[i] = (float *) malloc ((E->advection.markers+1)*sizeof(float));
        E->XMC[i] = (float *) malloc ((E->advection.markers+1)*sizeof(float));
        E->XMCpred[i] = (float *) malloc ((E->advection.markers+1)*sizeof(float));
        }
      E->C12 = (int *) malloc ((E->advection.markers+1)*sizeof(int));
      E->CElement = (int *) malloc ((E->advection.markers+1)*sizeof(int));
      }

    }

return;
}

void PG_timestep_particle(E)
     struct All_variables *E;
{
    void timestep();
    void predictor();
    void corrector();
    void pg_solver();
    void remove_horiz_ave();
    void std_timestep();
    void temperatures_conform_bcs();
    void thermal_buoyancy();
    void Runge_Kutta();
    void Euler();
    void get_fixed_temp();
    float Tmax(),T_interior1;
    int i,j,psc_pass,count,steps,iredo;
    int keep_going;

    float *T1, *Tdot1;
    static float *DTdot;
    FILE *fp;

    static int loops_since_new_eta = 0;
    static int been_here = 0;
    static int on_off = 0;

   if (been_here++==0)  {
      DTdot= (float *)malloc((E->mesh.nno+1)*sizeof(float));
      }

    if (on_off==0)  {
      E->advection.timesteps++;
      std_timestep(E);
      E->advection.total_timesteps++;
      }


    if (on_off==1)  {

         Runge_Kutta(E,E->XMC,E->XMCpred,E->C,E->V,E->VO);

        }

    else if (on_off==0)   {

       Euler(E,E->XMC,E->XMCpred,E->C,E->V,E->VO);

                                 /* update temperature    */

       predictor(E,E->T,E->Tdot,1);
       for(psc_pass=0;psc_pass<E->advection.temp_iterations;psc_pass++)   {
             pg_solver(E,E->T,E->Tdot,DTdot,E->V,1.0,E->TB,E->node);
             corrector(E,E->T,E->Tdot,DTdot,1);
             }

       temperatures_conform_bcs(E,E->T);

       }


    thermal_buoyancy(E);

    if( E->advection.timesteps < E->advection.max_timesteps)
        E->control.keep_going = 1;
    else
        E->control.keep_going = 0;

    on_off=(on_off==0)?1:0;

return;
}




void PG_timestep(E)
     struct All_variables *E;
{    
    void timestep();
    void predictor();
    void corrector();
    void pg_solver();
    void remove_horiz_ave();
    void std_timestep();
    void temperatures_conform_bcs();
    void thermal_buoyancy();
    float Tmax(),T_interior1;

    int i,j,psc_pass,count,steps,iredo;
    int keep_going;

    const int nno=E->mesh.nno;

    float *DTdot, *T1, *Tdot1;

    static int loops_since_new_eta = 0;
    static int been_here = 0;
   
    DTdot= (float *)malloc((E->mesh.nno+1)*sizeof(float));
    T1= (float *)malloc((E->mesh.nno+1)*sizeof(float));
    Tdot1= (float *)malloc((E->mesh.nno+1)*sizeof(float));

    if (been_here++ ==0)    {
	  E->advection.timesteps=0;
      }

    E->advection.timesteps++;
 
    std_timestep(E);

    for(i=1;i<=nno;i++) {
        T1[i] = E->T[i];
    }


   if (E->control.composition)         {
	     predictor(E,E->C,E->Cdot,0);
	     for(psc_pass=0;psc_pass<E->advection.temp_iterations;psc_pass++) {
	 	   pg_solver(E,E->C,E->Cdot,DTdot,E->V,E->control.comp_diff,E->CB,E->node);
		   corrector(E,E->C,E->Cdot,DTdot,0);
	       }	     
         }


    /* update temperature    */

       predictor(E,E->T,E->Tdot,1);
       for(psc_pass=0;psc_pass<E->advection.temp_iterations;psc_pass++)   {
         pg_solver(E,E->T,E->Tdot,DTdot,E->V,1.0,E->TB,E->node);
         corrector(E,E->T,E->Tdot,DTdot,1);
         }	     

    E->advection.total_timesteps++;

    temperatures_conform_bcs(E,E->T); 


    thermal_buoyancy(E);

    //compute CBF heat flux for boundary nodes, see Zhong_etal_PEPI_1992
//DUE to the WRONG operation, we delete Size_does_matter.c, so the MassForHF is lost, need to be re-program
//    heat_flux_CBF(E,T1,E->T,E->advection.timestep,1.0);
 

    if( E->advection.timesteps < E->advection.max_timesteps)
        E->control.keep_going = 1;
    else
 	    E->control.keep_going = 0;

    free((void *) DTdot );  
    free((void *) T1 );   
    free((void *) Tdot1 );
    
return;  
}


/* ==============================
   predictor and corrector steps.
   ============================== */

void predictor(E,field,fielddot,ic)
     struct All_variables *E;
     float *field,*fielddot;
     int ic;

{ 
    int node;
    float multiplier;

   multiplier = (1.0-E->advection.gamma) * E->advection.timestep;

  if (ic==1)
    for(node=1;node<=E->mesh.nno;node++)  {
	if(!(E->node[node] & (OFFSIDE | TBX | TBZ | TBY))) 
	    field[node] += multiplier * fielddot[node] ;
	fielddot[node] = 0.0;
       }
  else
    for(node=1;node<=E->mesh.nno;node++)  {
	if(!(E->node[node] & OFFSIDE )) 
	    field[node] += multiplier * fielddot[node] ;
	fielddot[node] = 0.0;
       }

   return; }

void corrector(E,field,fielddot,Dfielddot,ic)
     struct All_variables *E;
     float *field,*fielddot,*Dfielddot;
     int ic;
    
{  int node;
   float multiplier;

   multiplier = E->advection.gamma * E->advection.timestep;

 if (ic==1)
   for(node=1;node<=E->mesh.nno;node++) {
       if(!(E->node[node] & (OFFSIDE | TBX | TBZ | TBY)))
	   field[node] += multiplier * Dfielddot[node];
       fielddot[node] +=  Dfielddot[node]; 
       }
 else
   for(node=1;node<=E->mesh.nno;node++) {
       if(!(E->node[node] & OFFSIDE ))
	   field[node] += multiplier * Dfielddot[node];
       fielddot[node] +=  Dfielddot[node]; 
       }
   
   return;  
 }

/* ===================================================
   The solution step -- determine residual vector from
   advective-diffusive terms and solve for delta Tdot
   Two versions are available -- one for Cray-style 
   vector optimizations etc and one optimized for 
   workstations.
   =================================================== */


void pg_solver(E,T,Tdot,DTdot,V,diff,TBC,FLAGS)
     struct All_variables *E;
     float *T,*Tdot,*DTdot;
     float **V;
     float diff;
     float **TBC;
     unsigned int *FLAGS;
{
    void get_global_shape_fn();
    void pg_shape_fn();
    void element_residual();
    int el,e,a,i,a1;
    double xk[3][5],Eres[9];  /* correction to the (scalar) Tdot field */
    float tempdiff;

    struct Shape_function PG;
    struct Shape_function GN;
    struct Shape_function_dA dOmega;
    struct Shape_function_dx GNx;
 
    const int dims=E->mesh.nsd;
    const int dofs=E->mesh.dof;
    const int ends=enodes[dims];

    
    for(i=1;i<=E->mesh.nno;i++)
 	  DTdot[i] = 0.0;

    tempdiff = diff;
    for(el=1;el<=E->mesh.nel;el++)    {

	  get_global_shape_fn(E,el,&GN,&GNx,&dOmega,0,E->mesh.levmax);
	  pg_shape_fn(E,el,&PG,&GNx,V,tempdiff);
	  element_residual(E,el,PG,GNx,dOmega,V,T,Tdot,Eres,tempdiff);

      for(a=1;a<=ends;a++) {
	    a1 = E->ien[el].node[a];
	    	DTdot[a1] += Eres[a]; 
        }

      } /* next element */

    for(i=1;i<=E->mesh.nno;i++) {
 	  if(E->node[i] & OFFSIDE) continue;
	  DTdot[i] *= E->Mass[i];         /* lumped mass matrix */
      }
     
    return;    
}



/* ===================================================
   Petrov-Galerkin shape functions for a given element
   =================================================== */

void pg_shape_fn(E,el,PG,GNx,V,diffusion)
     struct All_variables *E;
     int el;
     struct Shape_function *PG;
     struct Shape_function_dx *GNx;
     float **V;
     float diffusion;

{ 
    int i,j,node;
    int *ienmatrix;

    double uc1,uc2,uc3;
    double u1,u2,u3,VV[4][9];
    double uxse,ueta,ufai,xse,eta,fai,dx1,dx2,dx3,adiff;

    double prod1,unorm,twodiff;
    
    const int dims=E->mesh.nsd;
    const int dofs=E->mesh.dof;
    const int lev=E->mesh.levmax;
    const int nno=E->mesh.nno;
    const int ends=enodes[E->mesh.nsd];
    const int vpts=vpoints[E->mesh.nsd];
  
    ienmatrix=E->ien[el].node;

    twodiff = 2.0*diffusion;
 
    uc1 =  uc2 = uc3 = 0.0;

    for(i=1;i<=ends;i++)   {
        node = ienmatrix[i];
        VV[1][i] = V[1][node];
        VV[2][i] = V[2][node];
      }
         
    for(i=1;i<=ENODES2D;i++) {
      uc1 +=  E->N.ppt[GNPINDEX(i,1)]*VV[1][i];
      uc2 +=  E->N.ppt[GNPINDEX(i,1)]*VV[2][i];
      }
    dx1 = 0.5*(E->X[1][ienmatrix[3]]+E->X[1][ienmatrix[4]]
              -E->X[1][ienmatrix[1]]-E->X[1][ienmatrix[2]]);
    dx2 = 0.5*(E->X[2][ienmatrix[3]]+E->X[2][ienmatrix[4]]
              -E->X[2][ienmatrix[1]]-E->X[2][ienmatrix[2]]);
    uxse = fabs(uc1*dx1+uc2*dx2); 


    dx1 = 0.5*(E->X[1][ienmatrix[2]]+E->X[1][ienmatrix[3]]
              -E->X[1][ienmatrix[1]]-E->X[1][ienmatrix[4]]);
    dx2 = 0.5*(E->X[2][ienmatrix[2]]+E->X[2][ienmatrix[3]]
              -E->X[2][ienmatrix[1]]-E->X[2][ienmatrix[4]]);
    ueta = fabs(uc1*dx1+uc2*dx2); 

    xse = (uxse>twodiff)? (1.0-twodiff/uxse):0.0;
    eta = (ueta>twodiff)? (1.0-twodiff/ueta):0.0;

    unorm = uc1*uc1 + uc2*uc2;

    adiff = (unorm>0.000001)?( (uxse*xse+ueta*eta)/(2.0*unorm) ):0.0;

    for(i=1;i<=VPOINTS2D;i++) {
          u1 = u2 = 0.0;
          for(j=1;j<=ENODES2D;j++)  /* this line heavily used */ {
            u1 += VV[1][j] * E->N.vpt[GNVINDEX(j,i)]; 
            u2 += VV[2][j] * E->N.vpt[GNVINDEX(j,i)];
	    }
	    
	  for(j=1;j<=ENODES2D;j++) {
             prod1 = (u1 * GNx->vpt[GNVXINDEX(0,j,i)] +
                      u2 * GNx->vpt[GNVXINDEX(1,j,i)]);
             PG->vpt[GNVINDEX(j,i)] = E->N.vpt[GNVINDEX(j,i)] + adiff * prod1;
	     }
    }

   return;
 }



/* ==========================================
   Residual force vector from heat-transport.
   Used to correct the Tdot term.
   =========================================  */

void element_residual(E,el,PG,GNx,dOmega,V,field,fielddot,Eres,diff)
     struct All_variables *E;
     int el;
     struct Shape_function PG;
     struct Shape_function_dA dOmega;
     struct Shape_function_dx GNx;
     float **V;
     float *field,*fielddot;
     double Eres[9];
     float diff;

{
    int i,j,a,k,node,nodes[4],d,aid,back_front,onedfns;
    double Q;
    double dT[9],VV[4][9];
    double tx1[9],tx2[9],tx3[9];
    double v1[9],v2[9],v3[9];
    double adv_dT,t2[4];
    double T,DT;
    static int been_here=0;
 
    register double prod,sfn;
    struct Shape_function1 GM;
    struct Shape_function1_dA dGamma;
    double temp;
    double rhoAtGauss[9];
    int topbot;

    void get_global_1d_shape_fn();
    void get_global_1d_shape_fn_simple();
 
    const int dims=E->mesh.nsd;
    const int dofs=E->mesh.dof;
    const int nno=E->mesh.nno;
    const int lev=E->mesh.levmax;
    const int ends=enodes[dims];
    const int vpts=vpoints[dims];
    const int diffusion = (diff != 0.0); 
    const int elz = E->mesh.elz;

    for(i=1;i<=vpts;i++)	{ 
      dT[i]=0.0;
      v1[i] = tx1[i]=  0.0;
      v2[i] = tx2[i]=  0.0;
      rhoAtGauss[i] = 0.0;
      }

    for(i=1;i<=ends;i++)   {
      node = E->ien[el].node[i];
        VV[1][i] = V[1][node];
        VV[2][i] = V[2][node];
      }
         
  
    for(j=1;j<=ends;j++)       {
      node = E->ien[el].node[j];
	  T = field[node];  
	  if(E->node[node] & (TBX | TBY | TBZ))
	    DT=0.0;
	  else
	    DT = fielddot[node];
	
	    for(i=1;i<=vpts;i++)  {
	 	  dT[i] += DT * E->N.vpt[GNVINDEX(j,i)];
		  tx1[i] +=  GNx.vpt[GNVXINDEX(0,j,i)] * T; 
		  tx2[i] +=  GNx.vpt[GNVXINDEX(1,j,i)] * T;   
		  sfn = E->N.vpt[GNVINDEX(j,i)];
		  v1[i] += VV[1][j] * sfn;
		  v2[i] += VV[2][j] * sfn;
		  rhoAtGauss[i] += E->rho_r[node]*E->N.vpt[GNVINDEX(j,i)];
	    }
      }



   Q=E->control.Q0;
   E->rad_heat.total=E->control.Q0;

/*
   if (diff>0.9)   {
     for(j=1;j<=ends;j++) 
       Q += E->C[E->ien[i].node[j]]; 
     Q = Q/ends;
     Q = Q*E->control.Q0;
     }
*/

//   Q = (E->rad_heat.total + E->heating_visc[el] - E->heating_adi[el])/E->heating_latent[el];

    /* construct residual from this information */


    if(diffusion){
       for(j=1;j<=ends;j++) {
          Eres[j]=0.0;
          for(i=1;i<=vpts;i++) 
		/*
	    if(E->control.compressibility == 0) 
              Eres[j] -= PG.vpt[GNVINDEX(j,i)] * dOmega.vpt[i] 
	             * (dT[i] + v1[i] * tx1[i] + v2[i] * tx2[i] - (E->rad_heat.total
				  +E->heating_visc[el] - E->heating_adi[el])/E->heating_latent[el]) +
                     diff/E->heating_latent[el]*dOmega.vpt[i] 
                     * (GNx.vpt[GNVXINDEX(0,j,i)] * tx1[i] +
		        GNx.vpt[GNVXINDEX(1,j,i)] * tx2[i] ); 
	    else
              Eres[j] -= PG.vpt[GNVINDEX(j,i)] * dOmega.vpt[i] 
	             * (dT[i] + v1[i] * tx1[i] + v2[i] * tx2[i] - (E->rad_heat.total
				  + E->heating_visc[el]/rhoAtGauss[i] - E->heating_adi[el])/E->heating_latent[el])
                  +   diff/E->heating_latent[el]*dOmega.vpt[i] 
                     * (GNx.vpt[GNVXINDEX(0,j,i)] * tx1[i] + GNx.vpt[GNVXINDEX(1,j,i)] * tx2[i] 
		   		 + E->control.compressibility * E->N.vpt[GNVINDEX(j,i)] * tx2[i])/rhoAtGauss[i]; 
		*/
              Eres[j] -= PG.vpt[GNVINDEX(j,i)] * dOmega.vpt[i] 
	             * (dT[i] + v1[i] * tx1[i] + v2[i] * tx2[i] - (E->rad_heat.total
				  + E->heating_visc[el].node[i]/rhoAtGauss[i] - E->heating_adi[el].node[i])/E->heating_latent[el])
                  +   diff/E->heating_latent[el]*dOmega.vpt[i] 
                     * (GNx.vpt[GNVXINDEX(0,j,i)] * tx1[i] + GNx.vpt[GNVXINDEX(1,j,i)] * tx2[i] 
		   		 + E->control.compressibility * E->N.vpt[GNVINDEX(j,i)] * tx2[i])/rhoAtGauss[i]; 
	    }
      }
/*
    else { // no diffusion term 
	    for(j=1;j<=ends;j++) {
    	  Eres[j]=0.0;
		  for(i=1;i<=vpts;i++) 
		    Eres[j] -= PG.vpt[GNVINDEX(j,i)] * dOmega.vpt[i] * (dT[i] - Q + v1[i] * tx1[i] + v2[i] * tx2[i]);
		  }
      }	
*/
	/* See brooks etc: the diffusive term is excused upwinding for 
	   rectangular elements  */
  
    /* include BC's for fluxes at (nominally horizontal) edges (X-Y plane) */
	onedfns=0;
	topbot = 0;
	for(a=1;a<=ends;a++)
	    if (2==dims && (E->NODE[lev][E->ien[el].node[a]] & FBZ)) {
		Eres[a] += 1/exp(E->control.compressibility)*E->TB[2][E->ien[el].node[a]]*0.5/elz;
		/*
		if (!onedfns++) get_global_1d_shape_fn(E,el,&GM,&dGamma,topbot);
		if(topbot==0 && a==1) 
		    for(j=1;j<=onedvpoints[dims];j++) {
			Eres[a] += 1/exp(E->control.compressibility)*dGamma.vpt[GMVGAMMA(topbot+1,j)] *
			    E->M.vpt[GMVINDEX(1,j)] * g_1d[j].weight[dims-1] *
			    E->TB[2][E->ien[el].node[a]] * E->M.vpt[GMVINDEX(1,j)];
		    }
		else if(topbot==0 && a==4)
		    for(j=1;j<=onedvpoints[dims];j++) {
			Eres[a] += 1/exp(E->control.compressibility)*dGamma.vpt[GMVGAMMA(topbot+1,j)] *
			    E->M.vpt[GMVINDEX(2,j)] * g_1d[j].weight[dims-1] *
			    E->TB[2][E->ien[el].node[a]] * E->M.vpt[GMVINDEX(2,j)];
		    }
		  */

	    }


//old one for fluxes
/*    if(FLAGS!=NULL) {
	onedfns=0;
	for(a=1;a<=ends;a++)
	    if (FLAGS[E->ien[el].node[a]] & FBZ) {
		if (!onedfns++) get_global_1d_shape_fn(E,el,&GM,&dGamma);
 
		nodes[1] = loc[loc[a].node_nebrs[0][0]].node_nebrs[2][0];
		nodes[2] = loc[loc[a].node_nebrs[0][1]].node_nebrs[2][0];
		nodes[4] = loc[loc[a].node_nebrs[0][0]].node_nebrs[2][2];
		nodes[3] = loc[loc[a].node_nebrs[0][1]].node_nebrs[2][2];
	  
		for(aid=0,j=1;j<=onedvpoints[E->mesh.nsd];j++)
		    if (a==nodes[j])
			aid = j;
		if(aid==0)  
		    printf("%d: mixed up in pg-flux int: looking for %d\n",el,a);

		if (loc[a].plus[1] != 0)
		    back_front = 0;
		else back_front = dims;

		for(j=1;j<=onedvpoints[dims];j++)
		    for(k=1;k<=onedvpoints[dims];k++)
			Eres[a] += dGamma.vpt[GMVGAMMA(1+back_front,j)] *
			    E->M.vpt[GMVINDEX(aid,j)] * g_1d[j].weight[dims-1] *
			    BC[2][E->ien[el].node[a]] * E->M.vpt[GMVINDEX(k,j)];
	    }
    } 
    
 */
    return; 
}




/* =====================================================
   Obtain largest possible timestep (no melt considered)
   =====================================================  */

//NOTE  timestep
void std_timestep(E)
     struct All_variables *E;

{ 
    static int been_here = 0;
    static float diff_timestep,root3,root2;
    int i,d,n,nel,el,node;

    float adv_timestep;
    float ts,uc1,uc2,uc3,uc,size,step,VV[4][9];
    
    const int dims=E->mesh.nsd;
    const int dofs=E->mesh.dof;
    const int nno=E->mesh.nno;
    const int lev=E->mesh.levmax;
    const int ends=enodes[dims];
    
	nel=E->mesh.nel;

    if(E->advection.fixed_timestep != 0.0) {
      E->advection.timestep = E->advection.fixed_timestep;
      return;
    }

    if (been_here == 0)  {
	  diff_timestep = 1.0e8; 
	  for(el=1;el<=nel;el++)  { 
	      ts = E->eco[el].size[1]*E->eco[el].size[1];
	      diff_timestep = min(diff_timestep,ts);
	      ts = E->eco[el].size[2]*E->eco[el].size[2];
	      diff_timestep = min(diff_timestep,ts);
	    }
      diff_timestep = 0.5 * diff_timestep;
      }


    
  adv_timestep = 1.0e8;
  for(el=1;el<=nel;el++) {

	for(i=1;i<=ends;i++)   {
      node = E->ien[el].node[i];
        VV[1][i] = E->V[1][node];
        VV[2][i] = E->V[2][node];
        if(dims==3) VV[3][i] = E->V[3][node];
      }
         
    uc=uc1=uc2=uc3=0.0;	  
    if(3==dims) {
      for(i=1;i<=ENODES3D;i++) {
        uc1 += E->N.ppt[GNPINDEX(i,1)]*VV[1][i];
        uc2 += E->N.ppt[GNPINDEX(i,1)]*VV[2][i];
        uc3 += E->N.ppt[GNPINDEX(i,1)]*VV[3][i];
        }
      uc = fabs(uc1)/E->eco[el].size[1] + fabs(uc2)/E->eco[el].size[2] + fabs(uc3)/E->eco[el].size[3];

	  step = (0.5/uc);
	  adv_timestep = min(adv_timestep,step);
      }
    else {
	  for(i=1;i<=ENODES2D;i++) {
        uc1 += E->N.ppt[GNPINDEX(i,1)]*VV[1][i];
        uc2 += E->N.ppt[GNPINDEX(i,1)]*VV[2][i];
        }
      uc = fabs(uc1)/E->eco[el].size[1] + fabs(uc2)/E->eco[el].size[2];
	  
	  step = (0.5/uc);
	  adv_timestep = min(adv_timestep,step);
      }
    }

    adv_timestep = E->advection.dt_reduced * adv_timestep;         

    adv_timestep =  1.0e-32+min(E->advection.fine_tune_dt*adv_timestep,diff_timestep);
    
    E->advection.timestep = adv_timestep;

//    fprintf(stderr,"time advance = %f\n",adv_timestep);

    return; 
  }


  void process_heating(E)
  struct All_variables *E;
  { 

  int e,i,j,ee;
  static int been=0;
  static float para1;
  double slope,temp1,temp2,temp3;
  FILE *fp;
  char filename[250];
  double rho_r,nodeHeating[9];
  void return_horiz_ave();
  

    const int dims = E->mesh.nsd;
    const int ends = enodes[dims];
    const int lev = E->mesh.levmax;
    const int nno = E->mesh.nno;
    const int vpts = vpoints[dims];
    const double two=2.0;

//  void strain_rate_2_inv();
  void strain_rate_2_inv_special();
  void heat_gint_to_nodes();



  return_horiz_ave(E,E->T,E->Have.T);

  E->rad_heat.total = E->control.Q0; 

  slope = (E->data.therm_exp_factor-1.0)/(E->sphere.ro-E->sphere.ri);

  temp1 = E->data.disptn_number/E->control.Ra_temp;

  for (e=1;e<=E->mesh.nel;e++)  {
    E->heating_latent[e] = 1.0;
    }

 if (E->control.visc_heating)  {
   strain_rate_2_inv_special(E,E->heating_visc,0);
  
   heat_gint_to_nodes(E,E->heating_visc,E->strain_rate_2,lev);
 
   for (e=1;e<=E->mesh.nel;e++) { 
     for (i=1;i<=vpts;i++) {
	    E->heating_visc[e].node[i] = temp1*E->EVi[(e-1)*vpts+i]*E->heating_visc[e].node[i];
     }
   }
   
 }

 if (E->control.adi_heating)  {
   for (e=1;e<=E->mesh.nel;e++)  {

    for (j=1;j<=ends;j++)    {
      i = E->ien[e].node[j];
      nodeHeating[j] = E->V[2][i]*(E->T[i]+E->control.Ts )*E->data.disptn_number;  
      }
    for (i=1;i<=vpts;i++) {
	    E->heating_adi[e].node[i] = 0.0;
	    for (j=1;j<=ends;j++) {
		    E->heating_adi[e].node[i] += nodeHeating[j] * E->N.vpt[GNVINDEX(j,i)];
	    }
	    
    }
    
    for(i=1;i<=vpts;i++) {
    	E->heating_adi[e].node[i] = E->heating_adi[e].node[i]*(slope*(E->eco[e].centre[2]-E->sphere.ri) + 1.0);
    }

    } //end for e
  } //end for if

//if(E->monitor.solution_cycles%(10*E->control.record_every)==0)  {
/*
if(0)  {
  sprintf(filename,"%s/heating.%d",E->control.data_file,E->monitor.solution_cycles);  
  fp=fopen(filename,"w");
  fprintf(fp,"QQ %g %g %g\n",E->control.Ra_temp,E->data.disptn_number,E->rad_heat.total);
  for (e=1;e<=E->mesh.nel;e++) {
	  temp1 = temp2 = temp3 = 0.0;
	  for(j=1;j<=vpts;j++) {
		  temp1 += E->heating_visc[e].node[j];
		  temp2 += E->heating_adi[e].node[j];
		  temp3 += E->EVi[(e-1)*vpts+j];
	  }
	  temp1 /= vpts;
	  temp2 /= vpts;
	  temp3 /= vpts;
    fprintf(fp,"%d %g %g %g\n",e,temp3,temp1,temp2);
    //fprintf(fp,"%g %g\n",temp1,temp2);
  }
  fclose(fp);
}
*/

/*
 if (E->control.visc_heating)  {
   strain_rate_2_inv(E,E->heating_visc,0);
 
   for (e=1;e<=E->mesh.nel;e++)  {
     temp2 = 0.0;
     for (i=1;i<=vpts;i++)
       temp2 += E->EVi[(e-1)*vpts+i];
     temp2 = temp2/vpts;

     E->heating_visc[e] = temp1*temp2*E->heating_visc[e];
     }
   }

 if (E->control.adi_heating)  {
   for (e=1;e<=E->mesh.nel;e++)  {
    ee = (e-1)%E->mesh.elz+1;

    temp2 = 0.0;
    for (i=1;i<=ends;i++)    {
      j = E->ien[e].node[i];
        temp2 = temp2 + E->V[2][j]*(E->T[j]+E->control.Ts )*E->data.disptn_number;  
      }
    temp2 = temp2/ends;
    E->heating_adi[e] = temp2*(slope*(E->eco[e].centre[2]-E->sphere.ri) + 1.0);

    }
  }

  if (E->control.Ra_670!=0.0)   {
    temp1 = 2.0*E->control.clapeyron670*E->control.Ra_670/(E->control.Ra_temp*E->control.width670);
    for (e=1;e<=E->mesh.nel;e++)  {
      temp2 = 0;
      temp3 = 0;
      for (i=1;i<=ends;i++)    {
        j = E->ien[e].node[i];
        temp2 = temp2 + temp1*(1.0-E->Fas670[j])*E->Fas670[j]
                   *E->V[2][j]*(E->T[j]+E->control.Ts)*E->data.disptn_number;
        temp3 = temp3 + temp1*E->control.clapeyron670
                        *(1.0-E->Fas670[j])*E->Fas670[j]
                        *(E->T[j]+E->control.Ts)*E->data.disptn_number;
        }
      temp2 = temp2/ends;
      temp3 = temp3/ends;
      E->heating_adi[e] += temp2;
      E->heating_latent[e] += temp3;
      }
    }

  if (E->control.Ra_410!=0.0)   {
    temp1 = 2.0*E->control.clapeyron410*E->control.Ra_410/(E->control.Ra_temp*E->control.width410);
    for (e=1;e<=E->mesh.nel;e++)  {
      temp2 = 0;
      temp3 = 0;
      for (i=1;i<=ends;i++)    {
        j = E->ien[e].node[i];
        temp2 = temp2 + temp1*(1.0-E->Fas410[j])*E->Fas410[j]
                   *E->V[2][j]*(E->T[j]+E->control.Ts)*E->data.disptn_number;
        temp3 = temp3 + temp1*E->control.clapeyron410
                       *(1.0-E->Fas410[j])*E->Fas410[j]
                       *(E->T[j]+E->control.Ts)*E->data.disptn_number;
        }
      temp2 = temp2/ends;
      temp3 = temp3/ends;
      E->heating_adi[e] += temp2;
      E->heating_latent[e] += temp3;
      }
    }

  fprintf(E->fp,"QQ %g \n",E->rad_heat.total);

if(E->monitor.solution_cycles%1000==0)  {
  sprintf(filename,"%s/heating.%d",E->control.data_file,E->monitor.solution_cycles);  
  fp=fopen(filename,"w");
  fprintf(fp,"QQ %g %g %g\n",E->control.Ra_temp,E->data.disptn_number,E->rad_heat.total);
  for (e=1;e<=E->mesh.nel;e++)
    fprintf(fp,"%d %g %g %g\n",e,E->EVi[(e-1)*vpts+1],E->heating_visc[e],E->heating_adi[e]);
 fclose(fp);
   }
  return; 
*/
  }


/*
  void process_heating_old(E)
  struct All_variables *E;
  { 

  int e,i,j,ee;
  static int been=0;
  static float para1;
  double slope,temp1,temp2,temp3;
  FILE *fp;
  char filename[250];
  void return_horiz_ave();

    const int dims = E->mesh.nsd;
    const int ends = enodes[dims];
    const int lev = E->mesh.levmax;
    const int nno = E->mesh.nno;
    const int vpts = vpoints[dims];
    const double two=2.0;

  void strain_rate_2_inv();

  if (been==0)  {
    para1 = E->sphere.ro_dim*E->sphere.ro_dim/(E->data.density*E->data.Cp*E->data.ref_temperature*E->data.therm_diff);
    E->data.disptn_number = E->data.therm_exp*E->data.grav_acc*E->sphere.ro_dim/E->data.Cp;
    been++;
    }

  return_horiz_ave(E,E->T,E->Have.T);

  if (E->rad_heat.num==0)  {
     E->rad_heat.total = para1*E->control.Q0; 
  }
  else {
     temp1 = 4.6e9-E->monitor.time_scale*E->monitor.elapsed_time;
     temp2 = E->rad_heat.percent[0]*E->rad_heat.concen[0]*E->rad_heat.heat_g[0]
	         *exp(temp1*log(two)/E->rad_heat.decay_t[0])
           + E->rad_heat.percent[1]*E->rad_heat.concen[1]*E->rad_heat.heat_g[1]
	         *exp(temp1*log(two)/E->rad_heat.decay_t[1])
           + E->rad_heat.percent[2]*E->rad_heat.concen[2]*E->rad_heat.heat_g[2]
	         *exp(temp1*log(two)/E->rad_heat.decay_t[2])
           + E->rad_heat.percent[3]*E->rad_heat.concen[3]*E->rad_heat.heat_g[3]
	         *exp(temp1*log(two)/E->rad_heat.decay_t[3]);
     E->rad_heat.total = para1*E->data.density*temp2; 
  }

  slope = (E->data.therm_exp_factor-1.0)/(E->sphere.ro-E->sphere.ri);

  temp1 = E->data.disptn_number/E->control.Ra_temp;

  for (e=1;e<=E->mesh.nel;e++)  {
    E->heating_latent[e] = 1.0;
    }

 if (E->control.visc_heating)  {
   strain_rate_2_inv(E,E->heating_visc,0);
 
   for (e=1;e<=E->mesh.nel;e++)  {
     temp2 = 0.0;
     for (i=1;i<=vpts;i++)
       temp2 += E->EVi[(e-1)*vpts+i];
     temp2 = temp2/vpts;

     E->heating_visc[e] = temp1*temp2*E->heating_visc[e];
     }
   }

 if (E->control.adi_heating)  {
   for (e=1;e<=E->mesh.nel;e++)  {
    ee = (e-1)%E->mesh.elz+1;

    temp2 = 0.0;
    for (i=1;i<=ends;i++)    {
      j = E->ien[e].node[i];
//      temp2 = temp2 + E->V[2][j]*(E->T[j]+E->control.Ts )*E->data.disptn_number;  
      temp2 = temp2 + E->V[2][j]*((E->Have.T[ee]+E->Have.T[ee+1])*0.5+E->control.Ts )*E->data.disptn_number;
      }
    temp2 = temp2/ends;
    E->heating_adi[e] = temp2*(slope*(E->eco[e].centre[2]-E->sphere.ri) + 1.0);

    }
  }

  if (E->control.Ra_670!=0.0)   {
    temp1 = 2.0*E->control.clapeyron670*E->control.Ra_670/(E->control.Ra_temp*E->control.width670);
    for (e=1;e<=E->mesh.nel;e++)  {
      temp2 = 0;
      temp3 = 0;
      for (i=1;i<=ends;i++)    {
        j = E->ien[e].node[i];
        temp2 = temp2 + temp1*(1.0-E->Fas670[j])*E->Fas670[j]
                   *E->V[2][j]*(E->T[j]+E->control.Ts)*E->data.disptn_number;
        temp3 = temp3 + temp1*E->control.clapeyron670
                        *(1.0-E->Fas670[j])*E->Fas670[j]
                        *(E->T[j]+E->control.Ts)*E->data.disptn_number;
        }
      temp2 = temp2/ends;
      temp3 = temp3/ends;
      E->heating_adi[e] += temp2;
      E->heating_latent[e] += temp3;
      }
    }

  if (E->control.Ra_410!=0.0)   {
    temp1 = 2.0*E->control.clapeyron410*E->control.Ra_410/(E->control.Ra_temp*E->control.width410);
    for (e=1;e<=E->mesh.nel;e++)  {
      temp2 = 0;
      temp3 = 0;
      for (i=1;i<=ends;i++)    {
        j = E->ien[e].node[i];
        temp2 = temp2 + temp1*(1.0-E->Fas410[j])*E->Fas410[j]
                   *E->V[2][j]*(E->T[j]+E->control.Ts)*E->data.disptn_number;
        temp3 = temp3 + temp1*E->control.clapeyron410
                       *(1.0-E->Fas410[j])*E->Fas410[j]
                       *(E->T[j]+E->control.Ts)*E->data.disptn_number;
        }
      temp2 = temp2/ends;
      temp3 = temp3/ends;
      E->heating_adi[e] += temp2;
      E->heating_latent[e] += temp3;
      }
    }

  fprintf(E->fp,"QQ %g \n",E->rad_heat.total);

if(E->monitor.solution_cycles%1000==0)  {
  sprintf(filename,"%s/heating.%d",E->control.data_file,E->monitor.solution_cycles);  
  fp=fopen(filename,"w");
  fprintf(fp,"QQ %g %g %g\n",E->control.Ra_temp,E->data.disptn_number,E->rad_heat.total);
  for (e=1;e<=E->mesh.nel;e++)
    fprintf(fp,"%d %g %g %g\n",e,E->EVi[(e-1)*vpts+1],E->heating_visc[e],E->heating_adi[e]);
 fclose(fp);
   }
  return; 
  }
*/

