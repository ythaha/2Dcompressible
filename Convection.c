/* Assumes parameter list is opened and reads the things it needs. 
   Variables are initialized etc, default values are set */


#include <math.h>
#include <malloc.h>
#include <sys/types.h>
#include "element_definitions.h"
#include "global_defs.h"
#include <stdlib.h> /* for "system" command */
#include <strings.h>

void set_convection_defaults(E)
     struct All_variables *E;
{
    void PG_timestep_with_melting();
    void PG_timestep();
    void PG_timestep_particle();
    void read_convection_settings();
    void convection_derived_values();
    void convection_allocate_memory();
    void convection_boundary_conditions();
    void node_locations();
    void convection_initial_fields();
    void twiddle_thumbs();

   input_int("composition",&(E->control.composition),"0");
   input_float("comp_diffusivity",&(E->control.comp_diff),"0");
   input_string("comp_adv_method",E->control.comp_adv_method,NULL);

   input_float("compressibility",&(E->control.compressibility),"0");
   input_float("relative_err_accuracy",&(E->control.relative_err_accuracy),"0.01");
   input_int("dyn_p_in_buoyancy",&(E->control.dyn_p_in_buoyancy),"0");



    if (E->control.composition && strcmp(E->control.comp_adv_method,"field")==0)
      E->next_buoyancy_field = PG_timestep;
    else if (E->control.composition && strcmp(E->control.comp_adv_method,"particle")==0)  {
      E->next_buoyancy_field = PG_timestep_particle;
      }
    else {
      E->next_buoyancy_field = PG_timestep;
      }

 
    E->special_process_new_buoyancy = twiddle_thumbs; 
    E->problem_settings = read_convection_settings;
    E->problem_derived_values = convection_derived_values;
    E->problem_allocate_vars = convection_allocate_memory;
    E->problem_boundary_conds = convection_boundary_conditions;
    E->problem_initial_fields = convection_initial_fields;
    E->problem_node_positions = node_locations;
    E->problem_update_node_positions = twiddle_thumbs;
    E->problem_update_bcs = twiddle_thumbs;

    sprintf(E->control.which_data_files,"Temp,Strf,Pres");
    sprintf(E->control.which_horiz_averages,"Temp,Visc,Vrms");
    sprintf(E->control.which_running_data,"Step,Time,");
    sprintf(E->control.which_observable_data,"Shfl");
 
return;
}

void read_convection_settings(E)
     struct All_variables *E;
    
{ 
    void advection_diffusion_parameters();
  float density_diff;
    
/* parameters */

    input_float("rayleigh",&(E->control.Ra_temp),"essential");
// NOTE E->data.ref_viscosity 是计算出来的，但也能输入，？？？
    E->data.ref_viscosity = E->data.grav_acc*E->data.density*E->data.therm_exp
                  *E->data.ref_temperature*E->sphere.ro_dim*E->sphere.ro_dim*E->sphere.ro_dim
                  /(E->control.Ra_temp*E->data.therm_diff);

    input_float("rayleigh_comp",&(E->control.Ra_comp),"essential");

    density_diff=E->control.Ra_comp*E->data.ref_viscosity*E->data.therm_diff/(E->data.grav_acc*E->sphere.ro_dim*E->sphere.ro_dim*E->sphere.ro_dim);

    fprintf(E->fp,"Ra_temp=%.5e Ra_comp=%.5e %.5e %.5e\n",E->control.Ra_temp,E->control.Ra_comp,E->data.ref_viscosity,density_diff);

    input_boolean("halfspace",&(E->convection.half_space_cooling),"off");
    input_float("halfspage",&(E->convection.half_space_age),"nodefault");
    
    input_int("temperature_blobs",&(E->convection.temp_blobs),"0");
    input_float_vector("temperature_blobx",E->convection.temp_blobs,E->convection.temp_blob_x);
    input_float_vector("temperature_bloby",E->convection.temp_blobs,E->convection.temp_blob_y);
    input_float_vector("temperature_blobz",E->convection.temp_blobs,E->convection.temp_blob_z);
    input_float_vector("temperature_blobsize",E->convection.temp_blobs,E->convection.temp_blob_radius);
    input_float_vector("temperature_blobDT",E->convection.temp_blobs,E->convection.temp_blob_T);
    input_float_vector("temperature_blobbg",E->convection.temp_blobs,E->convection.temp_blob_bg);
    input_int_vector("temperature_blobsticky",E->convection.temp_blobs,E->convection.temp_blob_sticky);
    
    input_int("temperature_zones",&(E->convection.temp_zones),"0");
    input_float_vector("temperature_zonex1",E->convection.temp_zones,E->convection.temp_zonex1);
    input_float_vector("temperature_zonex2",E->convection.temp_zones,E->convection.temp_zonex2);
    input_float_vector("temperature_zonez1",E->convection.temp_zones,E->convection.temp_zonez1);
    input_float_vector("temperature_zonez2",E->convection.temp_zones,E->convection.temp_zonez2);
    input_float_vector("temperature_zoney1",E->convection.temp_zones,E->convection.temp_zoney1);
    input_float_vector("temperature_zoney2",E->convection.temp_zones,E->convection.temp_zoney2);
    input_float_vector("temperature_zoney2",E->convection.temp_zones,E->convection.temp_zoney2);
    input_float_vector("temperature_zoney2",E->convection.temp_zones,E->convection.temp_zoney2);
    input_float_vector("temperature_zonehw",E->convection.temp_zones,E->convection.temp_zonehw);
    input_float_vector("temperature_zonemag",E->convection.temp_zones,E->convection.temp_zonemag);
    input_int_vector("temperature_zonesticky",E->convection.temp_zones,E->convection.temp_zone_sticky);
    
    input_int("num_perturbations",&(E->convection.number_of_perturbations),"0,0,32");
    input_float_vector("perturbmag",E->convection.number_of_perturbations,E->convection.perturb_mag);
    input_float_vector("perturbk",E->convection.number_of_perturbations,E->convection.perturb_k);
    
    advection_diffusion_parameters(E);
    
    if (E->control.restart)    {
       input_int("restart_timesteps",&(E->monitor.solution_cycles),"0");
       input_string("oldfile",E->convection.old_T_file,"initialize");
       }

    return;
}

/* =================================================================
   Any setup which relates only to the convection stuff goes in here
   ================================================================= */

void convection_derived_values(E)  
     struct All_variables *E;
 
{ 

return;
}

void convection_allocate_memory(E)
     struct All_variables *E;

{ void advection_diffusion_allocate_memory();

  advection_diffusion_allocate_memory(E);

return;
}

/* ============================================ */
    
void convection_initial_fields(E)
     struct All_variables *E;

{ 
    void convection_initial_temperature();
    void convection_initial_markers();
    void convection_initial_density();

    convection_initial_density(E);
    report(E,"convection, initial temperature");
    convection_initial_temperature(E);
    if (E->control.composition && !(strcmp(E->control.comp_adv_method,"field")==0))    {
       convection_initial_markers(E);
       }

  return; }

/* =========================================== */

void convection_boundary_conditions(E)
     struct All_variables *E;

{
    void velocity_boundary_conditions();
    void temperature_boundary_conditions();
    void temperatures_conform_bcs();
    void composition_boundary_conditions();
 
    velocity_boundary_conditions(E);      /* universal */
    temperature_boundary_conditions(E);

    temperatures_conform_bcs(E);

    composition_boundary_conditions(E);

    return;
}

/* ===============================
   Initialization of fields .....
   =============================== */

void convection_initial_density(E)
	     struct All_variables *E;
{
    int node;
    int nox,noy,noz;
    int i,j,k;

    noy=E->mesh.noy;
    noz=E->mesh.noz;
    nox=E->mesh.nox;


    for(i=1;i<=noy;i++)
        for(j=1;j<=nox;j++)
              for(k=1;k<=noz;k++)  {
                  node=k+(j-1)*noz+(i-1)*nox*noz;
		  E->rho_r[node] = exp( (1-E->X[2][node])*E->control.compressibility );

//		  printf("%f \n",E->rho_r[node]);
	      }
}

void convection_initial_temperature(E)
     struct All_variables *E;
{
    int i,j,k,p,node,ii,jj;
    double temp,base,radius,radius2;
    double modified_plgndr_a(),drand48();
    FILE *fp;
    void remove_horiz_ave();
    void temperatures_conform_bcs();
    void thermal_buoyancy();
    void process_restart_tc();
    
    int in1,in2,in3,instance,nox,noy,noz,nfz,ok,noz2,ll,mm;
    char output_file[255];
    double amm,tbase,tbase1,t1,r1,weight,para1,plate_velocity,delta_temp,age;
    double x00,x01,x02,slope,con;
    double bmax;

    int read_previous_field();
    int tempcount=0;

    const int dims=E->mesh.nsd;
    const float e_5=1.0e-5;

    noy=E->mesh.noy;  
    noz=E->mesh.noz;  
    nox=E->mesh.nox;  

    para1 = E->control.Ts*E->data.ref_temperature + 0.4*E->data.ref_temperature; 

    tbase = (para1 - E->control.Ts*E->data.ref_temperature)/E->data.ref_temperature;
    tbase1 = (para1 + 200 - E->control.Ts*E->data.ref_temperature)/E->data.ref_temperature;


        amm = E->convection.perturb_k[0];

/*
        noz2 = (noz-1)/2+1;
        con = (noz-1)/(E->sphere.ro-E->sphere.ri);
*/
        con = E->convection.perturb_mag[0];

          for(i=1;i<=noy;i++)
            for(j=1;j<=nox;j++)
              for(k=1;k<=noz;k++)  {
                node=k+(j-1)*noz+(i-1)*nox*noz;
                t1=E->X[1][node];
                r1=E->X[2][node];

                E->T[node] = 0.0;
                E->C[node] = 0.0;

//                E->T[node] = 1-r1 + 
//                   con*cos(M_PI*amm*t1)*
//                   sin(M_PI*r1);

		
		if(E->control.stokes) {
			con=noz-1;
			if(k==(1+noz)/2) {
//				tempcount++;
//				printf("%d k == %d, amm = %f, con = %f\n",tempcount,k,amm,con);
				//E->T[node] = con*cos(M_PI*amm*t1);
				E->T[node] = con*cos(M_PI*amm*t1);
			}
//			E->T[node] =
//                                  con*cos(M_PI*amm*t1)*
//                                  sin(M_PI*r1);
		}
		else {
                    //E->T[node] = (1-r1 + con*drand48());
		    E->T[node] = (1-r1 + con*cos(M_PI*amm*t1/E->mesh.layer[1]));
		    //E->T[node] = 100+(con*cos(M_PI*amm*t1/1.5));

		    //special initial condition for Scott King's timedependent ALA cases
//		    bmax = 1.0 - E->control.Ts*(exp(E->data.disptn_number) - 1.0);
//		    E->T[node] = bmax/2.0 + 0.001*cos(t1*M_PI)*sin(r1*M_PI);
//		    E->T[node] += E->control.Ts*exp(E->data.disptn_number*(1.0-r1));
		}

/* 
*/

              if (r1<=E->viscosity.zcrust1)
                 E->C[node] = 1.0;
              else
                 E->C[node] = 0.0;

              E->C[node] = 0.0;


              E->node[node] = E->node[node] | (INTX | INTZ | INTY);
			
              }    /* close the loop for node */


   if (E->control.restart==1 )
         process_restart_tc(E,E->mesh.levmax);

   temperatures_conform_bcs(E);

/*
   sprintf(output_file,"%s.XandT",E->control.data_file);
   if ( (fp = fopen(output_file,"w")) != NULL) {
      for (j=1;j<=E->mesh.nno;j++)
         fprintf(fp,"X[%05d] = %.6e Z[%05d] = %.6e T[%05d] = %.6e C[%05d] = %.6e\n",j,E->X[1][j],j,E->X[2][j],j,E->T[j],j,E->C[j]);
      }        

   fclose(fp);     
*/

    thermal_buoyancy(E);

    return; 
    }

 void convection_initial_markers(E)
   struct All_variables *E;
  {
    int el,i,j,k,p,node,ii,jj;
    float half_dist, init_height,dx,dr;
    char input_s[100],output_file[255];
    FILE *fp;
    void  process_restart_mk();
    void get_C_from_markers();

    half_dist = E->advection.marker_maxdist/2.0;
    init_height = E->viscosity.z410;

   if (E->control.restart==1)
      process_restart_mk(E);
   else {

    node = 0;
    p = pow((double)E->advection.markers_per_ele,(double)(1.0/E->mesh.dof));
    for (el=1;el<=E->mesh.nel;el++)  {
      dx = (E->X[1][E->ien[el].node[3]] - E->X[1][E->ien[el].node[1]])/p;
      dr = (E->X[2][E->ien[el].node[3]] - E->X[2][E->ien[el].node[1]])/p;
      for (i=1;i<=p;i++)
      for (j=1;j<=p;j++)  {
        node ++;
        E->XMC[1][node] = E->X[1][E->ien[el].node[1]] + dx*(i-0.5);
        E->XMC[2][node] = E->X[2][E->ien[el].node[1]] + dr*(j-0.5);
        E->CElement[node] = el;
        if (E->XMC[2][node]>E->viscosity.zcrust1)
              E->C12[node] = 0;
        else
              E->C12[node] = 1;
        }
      }

    E->advection.markers = node;

    get_C_from_markers(E,E->C,E->CElement);
    }
 return;
  }

/* ====================================================================== */

void process_restart_tc(E,lev)
    struct All_variables *E;
   int lev;
{
    int fileid[20];
    int i,j,k,ii,size2;
 char input_s[200],output_file[255],in_file[255];
 FILE *fp;
 float t1;
 float temp;
   sprintf(output_file,"%s/temp.%d",E->convection.old_T_file,E->monitor.solution_cycles);
   fp=fopen(output_file,"r");

   fgets(input_s,200,fp);
   sscanf(input_s,"%d %d %g",&i,&E->advection.timesteps,&E->monitor.elapsed_time);
   for (i=1;i<=E->mesh.NNO[lev];i++)   {
     fgets(input_s,200,fp);
     sscanf(input_s,"%g %g %g %g",&temp,&t1,&t1,&t1);
//     E->T[i]=(temp+0.2)>1.0 ? 1.0:(temp+0.2);
//     E->T[i]=(temp-0.2)<0.0 ? 0.0:(temp-0.2);
     E->T[i]=temp;
     }


// E->monitor.solution_cycles = E->advection.timesteps;

 fclose(fp);

  return;
  }


/* ====================================================================== */

void process_restart_mk(E)
  struct All_variables *E;
{
 int fileid[20];
 int i,j,k,ii,size2;
 char input_s[200],output_file[255],in_file[255];
 FILE *fp;
 float t1;
 void get_C_from_markers();

 fp = fopen(E->convection.old_T_file,"r");
   fgets(input_s,200,fp);
   sscanf(input_s,"%d %d %g",&E->advection.markers,&E->advection.timesteps,&E->monitor.elapsed_time);
   for (i=1;i<=E->advection.markers;i++)  {
      fgets(input_s,100,fp);
      sscanf(input_s,"%g %g %d",&E->XMC[1][i],&E->XMC[2][i],&E->CElement[i]);
      }

 E->monitor.solution_cycles = E->advection.timesteps;

 get_C_from_markers(E,E->XMC,E->C,E->CElement);

 fclose(fp);

  return;
  }

