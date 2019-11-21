/* Set up the finite element problem to suit: returns with all memory */
/* allocated, temperature, viscosity, node locations and how to use */
/* them all established. 8.29.92 or 29.8.92 depending on your nationality*/
#include <signal.h>
#include <math.h>
#include <malloc.h>
#include <sys/types.h>
#include <string.h>
#include "element_definitions.h"
#include "global_defs.h"

int Emergency_stop;

void read_instructions(E,argc,argv)
     struct All_variables *E;
     int argc;
     char **argv;
{
    int get_process_identifier();

    void allocate_common_vars();  
    void common_initial_fields();
    void read_initial_settings();
    void global_default_values();
    void global_derived_values();
    void construct_ien();
    void construct_masks();
    void construct_shape_functions();
    void construct_id();
    void construct_lm();
    void construct_sub_element();
    void mass_matrix();
    void construct_node_ks();
    void construct_node_maps();
    void construct_mat_group();
    void interuption();
    void set_up_nonmg_aliases();
    void check_bc_consistency();
    void node_locations();
    void allocate_velocity_vars();                 

    void setup_parser();

  
    double start_time, CPU_time0(),vmag;
    double vdot();

    int *temp, i;

 
    /* =====================================================
       Global interuption handling routine defined once here
       =====================================================  */

    start_time=CPU_time0();
    Emergency_stop = 0;
    signal(SIGINT,interuption);
    signal(SIGTERM,interuption);

    E->control.PID=get_process_identifier(); 

    /* ==================================================
       Initialize from the command line 
       from startup files. (See Parsing.c).
       ==================================================  */

    setup_parser(E,argc,argv);

    global_default_values(E); 
    read_initial_settings(E);  
    (E->problem_derived_values)(E);   /* call this before global_derived_  */
    global_derived_values(E);

    allocate_common_vars(E);                 
    (E->problem_allocate_vars)(E);
    (E->solver_allocate_vars)(E);
 
           /* logical domain */
    construct_ien(E);
    construct_sub_element(E);

           /* physical domain */
    node_locations (E);             
    allocate_velocity_vars(E);                 
    (E->problem_boundary_conds)(E);
    
    check_bc_consistency(E);
 
    construct_masks(E);		/* order is important here */
    construct_id(E);
    construct_lm(E);
    construct_mat_group(E);

    construct_shape_functions(E);
    mass_matrix(E);
printf("ok1\n");
    (E->problem_initial_fields)(E);   /* temperature/chemistry/melting etc */
printf("ok2\n");
    common_initial_fields(E);  /* velocity/pressure/viscosity (viscosity must be done LAST) */

printf("ok3\n");
    shutdown_parser(E);
/*
*/
 
    return;
}


/* ===================================
   Functions which set up details 
   common to all problems follow ...
   ===================================  */

void allocate_common_vars(E) 
     struct All_variables *E;

{ 
    double **dmatrix();
    float **fmatrix();
    void set_up_nonmg_aliases();
    int i,j,l,nno_l,npno_l,nozl,nnov_l,nxyz;

    const int vpts=vpoints[E->mesh.nsd];
 
    E->P = (double *) malloc ((E->mesh.npno+1)*sizeof(double));
    E->S = (double *) malloc ((E->mesh.npno+1)*sizeof(double));
 
    E->C        = (float *) malloc((E->mesh.nno+1)*sizeof(float));
    E->CE       = (float *) malloc((E->mesh.nel+1)*sizeof(float));

    E->T        = (float *) malloc((E->mesh.nno+1)*sizeof(float));
    E->buoyancy = (float *) malloc((E->mesh.nno+1)*sizeof(float));
    E->rho_r    = (float *) malloc((E->mesh.nno+1)*sizeof(float));
    E->NP       = (float *) malloc((E->mesh.nno+1)*sizeof(float));
    E->edot     = (float *) malloc((E->mesh.nno+1)*sizeof(float));
    E->heatflux    = (float *) malloc((E->mesh.nno+1)*sizeof(float));
    E->heatflux_adv    = (float *) malloc((E->mesh.nno+1)*sizeof(float));
    E->MassForHF   = (float *) malloc((E->mesh.nno+1)*sizeof(float));

    E->heating_latent    = (float *) malloc((E->mesh.nel+1)*sizeof(float));
//    E->heating_visc    = (float *) malloc((E->mesh.nel+1)*sizeof(float));
//    E->heating_adi    = (float *) malloc((E->mesh.nel+1)*sizeof(float));
    E->heating_adi    = (struct HEATINGNODE *) malloc((E->mesh.nel+1)*sizeof(struct HEATINGNODE));
    E->heating_visc    = (struct HEATINGNODE *) malloc((E->mesh.nel+1)*sizeof(struct HEATINGNODE));
    //E->strain_rate_2   = (float *) malloc((E->mesh.nno+1)*sizeof(float));
    E->strain_rate_2   = (double *) malloc((E->mesh.nno+1)*sizeof(double));

    E->Fas670    = (float *) malloc((E->mesh.nno+1)*sizeof(float));
    E->Fas670_b  = (float *) malloc((E->mesh.nsf+1)*sizeof(float));
    E->Fas410    = (float *) malloc((E->mesh.nno+1)*sizeof(float));
    E->Fas410_b  = (float *) malloc((E->mesh.nsf+1)*sizeof(float));

  for(i=1;i<=E->mesh.nsd;i++)  {
    E->TB[i] = (float *)  malloc((E->mesh.nno+1)*sizeof(float));
    E->CB[i] = (float *)  malloc((E->mesh.nno+1)*sizeof(float));
    }

    E->Have.T    = (float *) malloc((E->mesh.noz+1)*sizeof(float));
    E->Have.adiT    = (float *) malloc((E->mesh.noz+1)*sizeof(float));
    E->Have.C    = (float *) malloc((E->mesh.noz+1)*sizeof(float));
    E->Have.F   = (float *) malloc((E->mesh.noz+1)*sizeof(float));
    E->Have.f    = (float *) malloc((E->mesh.noz+1)*sizeof(float));
    E->Have.Vi   = (float *) malloc((E->mesh.noz+1)*sizeof(float));
    E->Have.Rho  = (float *) malloc((E->mesh.noz+1)*sizeof(float));
    E->Have.vrms = (float *) malloc((E->mesh.noz+1)*sizeof(float));
    E->Have.ViscH = (float *) malloc((E->mesh.noz+1)*sizeof(float));
    E->Have.AdiH = (float *) malloc((E->mesh.noz+1)*sizeof(float));
    
    E->stress    = (float *)malloc((2*6*E->mesh.nsf+12)*sizeof(float));
    E->slice.tpg      = (float *)malloc((E->mesh.nsf+2)*sizeof(float));
    E->slice.tpgb     = (float *)malloc((E->mesh.nsf+2)*sizeof(float));
    E->slice.vline     = (float *)malloc((E->mesh.nsf+2)*sizeof(float));
    E->slice.vlinek    = (float *)malloc((E->mesh.nsf+2)*sizeof(float));
    E->slice.shflux    = (float *)malloc((E->mesh.nsf+2)*sizeof(float));
    E->slice.bhflux    = (float *)malloc((E->mesh.nsf+2)*sizeof(float));
    E->slice.cen_mflux = (float *)malloc((E->mesh.nsf+2)*sizeof(float));
    E->slice.vxsurf[1] = (float *)malloc((E->mesh.nsf+2)*sizeof(float));
    E->slice.vxsurf[2] = (float *)malloc((E->mesh.nsf+2)*sizeof(float));
   
    E->mat = (int *) malloc((E->mesh.nel+2)*sizeof(int));

  E->XP[1] = (float *)  malloc((E->mesh.nox+1)*sizeof(float));
  E->XP[2] = (float *)  malloc((E->mesh.noz+1)*sizeof(float));

  for(i=E->mesh.levmin;i<=E->mesh.levmax;i++)   {
    for(j=1;j<=E->mesh.nsd;j++)   {
      E->XX[i][j] = (float *)  malloc((E->mesh.NNO[i]+1)*sizeof(float));
      E->Interp[i][j] = (float *)  malloc((E->mesh.NNO[i]+1)*sizeof(float));
      }
    E->MASS[i]     = (float *) malloc((E->mesh.NNO[i]+1)*sizeof(float));

    E->ECO[i] = (struct COORD *) malloc((E->mesh.NNO[i]+2)*sizeof(struct COORD));
    E->ID[i]  = (struct ID *)    malloc((E->mesh.NNO[i]+2)*sizeof(struct ID));
    E->IEN[i] = (struct IEN *)   malloc((E->mesh.NEL[i]+2)*sizeof(struct IEN));
    E->EL[i]  = (struct SUBEL *) malloc((E->mesh.NEL[i]+2)*sizeof(struct SUBEL));
    E->LMD[i] = (struct LM *)    malloc((E->mesh.NEL[i]+2)*sizeof(struct LM));
    E->elt_del[i]=(struct EG *)  malloc((E->mesh.NEL[i]+1)*sizeof(struct EG));
    E->elt_c[i] = (struct EC *)  malloc((E->mesh.NEL[i]+1)*sizeof(struct EC));
    E->BPI[i]    = (double *)    malloc((E->mesh.NPNO[i]+1)*sizeof(double));

    E->EVI[i] = (float *) malloc((E->mesh.NEL[i]+2)*vpoints[E->mesh.nsd]*sizeof(float));

    E->TW[i]  = (float  *)       malloc((E->mesh.NNO[i]+2)*sizeof(float));
    E->VI[i]  = (float *)        malloc((E->mesh.NNO[i]+2)*sizeof(float));
    E->NODE[i] = (unsigned int *)malloc((E->mesh.NNO[i]+2)*sizeof(unsigned int));
      E->TWW[i] =(struct FNODE *) malloc((E->mesh.NEL[i]+2)*sizeof(struct FNODE));

    E->NEI[i].nels     = (int *) malloc((E->mesh.NNO[i]+2)*sizeof(int));
    E->NEI[i].lnode    = (int *) malloc((E->mesh.NNO[i]+2)*enodes[E->mesh.nsd]*sizeof(int));
    E->NEI[i].element  = (int *) malloc((E->mesh.NNO[i]+2)*enodes[E->mesh.nsd]*sizeof(int));   
    }

  for(i=E->mesh.levmin;i<=E->mesh.levmax;i++)  {
    if (E->mesh.nsd==2)  {
        nxyz = max(E->mesh.nox,E->mesh.noz);
        }
    else if (E->mesh.nsd==3)  {
        nxyz = max(E->mesh.nox*E->mesh.noz,E->mesh.nox*E->mesh.noy);
        nxyz = 2*max(nxyz,E->mesh.noz*E->mesh.noy);
        }

    E->sien         = (struct SIEN *) malloc((nxyz+2)*sizeof(struct SIEN));
    E->surf_element = (int *) malloc((nxyz+2)*sizeof(int));
    E->surf_node    = (int *)malloc((E->mesh.nsf+2)*sizeof(int));
    }

  for(i=1;i<=E->mesh.nno;i++)     {
    E->T[i] = E->buoyancy[i] = E->rho_r[i] = 0.0;
    for(j=1;j<=E->mesh.nsd;j++)
      E->TB[j][i] = 0.0;
    }
  
  for(l=E->mesh.levmin;l<=E->mesh.levmax;l++)
    for(i=1;i<=E->mesh.NNO[l];i++)  {
      E->NODE[l][i] = (INTX | INTY | INTZ);  /* and any others ... */
      E->VI[l][i] = 1.0;
      E->TW[l][i] = 0.0;
      for(j=1;j<=E->mesh.nsd;j++)
         E->XX[l][j][i] = 0.0;
      }

  for(i=1;i<E->mesh.nel;i++)   {
     E->mat[i]=1;
     E->heating_latent[i] = 0.0;
//     E->heating_visc[i] = E->heating_adi[i] = 0.0;
     for(j=1;j<=vpts;j++) {
	     E->heating_visc[i].node[j] = 0;
	     E->heating_adi[i].node[j] = 0;
     }
  }

  for(i=1;i<=E->mesh.npno;i++)
      E->P[i] = 0.0;

  set_up_nonmg_aliases(E); 
  
  return; 
  }
  
/*  =========================================================  */

void allocate_velocity_vars(E) 
     struct All_variables *E;

{ 
    int i,j,l;
 
  E->mesh.nnov = E->mesh.nno;
  E->mesh.NEQ[E->mesh.levmax] = E->mesh.nnov * E->mesh.nsd;

  E->F = (double *) malloc((E->mesh.nsd*E->mesh.nnov+1)*sizeof(double));
  E->U = (double *) malloc((E->mesh.nsd*E->mesh.nnov+1)*sizeof(double));
 
  for(i=1;i<=E->mesh.nsd;i++)  {
      E->V[i]  = (float *)  malloc((E->mesh.nnov+1)*sizeof(float));
      E->VB[i] = (float *)  malloc((E->mesh.nnov+1)*sizeof(float));
      }

  for(i=E->mesh.levmin;i<=E->mesh.levmax;i++)  {
    E->BI[i] = (double *) malloc((E->mesh.NEQ[i]+2)*sizeof(double)); 
    E->EQN[i] = (unsigned int *) malloc((E->mesh.NEQ[i]+2)*sizeof(unsigned int)); 
    }


  for(l=E->mesh.levmin;l<=E->mesh.levmax;l++)
    for(i=0;i<E->mesh.NEQ[l];i++) {
      E->BI[l][i]=0.0;
      E->EQN[l][i]=0;
      }

  for(i=0;i<E->mesh.NEQ[E->mesh.levmax];i++)
    E->U[i]=0.0;

  for(i=1;i<=E->mesh.nnov;i++)
    for(j=1;j<=E->mesh.nsd;j++)
       E->V[j][i] =
         E->VB[j][i] = 0.0;
 
  return;  
 }
  

/*  =========================================================  */


void interuption()
    
{  if (Emergency_stop++) exit(0);
   fprintf(stderr,"Cleaning up before exit\n");
   return;  }


void global_default_values(E)
     struct All_variables *E;
{
    FILE *fp;

  /* FIRST: values which are not changed routinely by the user */
 		
  E->control.v_steps_low = 10;
  E->control.v_steps_upper = 1;
  E->control.max_res_red_each_p_mg = 1.0e-3;
  E->control.accuracy = 1.0e-6;
  E->control.vaccuracy = 1.0e-8;
  E->control.true_vcycle=0;
  E->control.depth_dominated=0;
  E->control.eqn_zigzag=0;
  E->control.verbose=0; /* debugging/profiles */

  /* SECOND: values for which an obvious default setting is useful */

  E->control.ORTHO = 1; /* for orthogonal meshes by default */
  E->control.ORTHOZ = 1; /* for orthogonal meshes by default */

 
    E->control.KERNEL = 0;
    E->control.stokes=0;
    E->control.CONVECTION = 0;
    E->control.CART2D = 0;
    E->control.CART3D = 0;
    E->control.CART2pt5D = 0;
    E->control.AXI = 0;
    E->control.CONJ_GRAD = 0;
    E->control.NMULTIGRID = 0;
    E->control.EMULTIGRID = 0;
    E->control.COMPRESS = 1;
    E->control.augmented_Lagr = 0;
    E->control.augmented = 0.0;

    /* Default: all optional modules set to `off' */
    E->control.MELTING_MODULE = 0;
    E->control.CHEMISTRY_MODULE = 0;

    E->control.composition = 0;
    E->control.comp_diff = 0.0;

    E->control.GRID_TYPE=1;
    E->mesh.hwidth[1]=E->mesh.hwidth[2]=E->mesh.hwidth[3]=1.0; /* divide by this one ! */
    E->mesh.magnitude[1]=E->mesh.magnitude[2]=E->mesh.magnitude[3]=0.0;
    E->mesh.offset[1]=E->mesh.offset[2]=E->mesh.offset[3]=0.0;

  E->mesh.levmax=0;
  E->mesh.levmin=0;
  E->mesh.noz = 1;    E->mesh.noz = 1;
  E->mesh.noy = 1;    E->mesh.noy = 1;  

  E->monitor.T_interior=1.0;

  E->viscosity.guess = 0;
  sprintf(E->viscosity.old_file,"initialize");
 
  E->control.precondition = 0;	/* for larger visc contrasts turn this back on  */
  E->control.vprecondition = 1;	
 
  E->mesh.toptbc = 1; /* fixed t */
  E->mesh.bottbc = 1;
  E->mesh.topvbc = 0; /* stress */
  E->mesh.botvbc = 0;
  E->mesh.sidevbc=0;
  E->mesh.periodic_x=0; /* reflection is default*/
  E->mesh.periodic_y=0;
  E->control.VBXtopval=0.0;
  E->control.VBYtopval=0.0;
  E->control.VBXbotval=0.0;
  E->control.VBYbotval=0.0;

  E->data.layer_km = 2800.0; /* Earth, whole mantle defaults */
  E->data.grav_acc = 9.81;
  E->data.therm_exp = 3.28e-5;
  E->data.Cp = 1200.0;
  E->data.therm_diff = 8.0e-7;
  E->data.therm_cond = 3.168;
  E->data.density = 3340.0;
  E->data.res_density = 3295.0;  /* density when X = ... */
  E->data.res_density_X = 0.3;
  E->data.melt_density = 2800.0;
  E->data.permeability = 3.0e-10;
  E->data.density_above = 1030.0;    /* sea water */
  E->data.gas_const = 8.3;
  E->data.surf_heat_flux = 4.4e-2;
  E->data.grav_const = 6.673e-11;
  E->data.surf_temp = 0.0;
  E->data.youngs_mod = 1.0e11;
  E->data.Te = 0.0;
  E->data.T_sol0 = 1373.0;	/* Dave's values 1991 (for the earth) */
  E->data.Tsurf = 273.0;
  E->data.dTsol_dz = 3.4e-3 ;
  E->data.dTsol_dF = 440.0;
  E->data.dT_dz = 0.48e-3;
  E->data.delta_S = 250.0;
  E->data.ref_temperature = 2 * 1350.0; /* fixed temperature ... delta T */
    
  /* THIRD: you forgot and then went home, let's see if we can help out */
 
    sprintf(E->control.data_file,"citcom.tmp.%d",getpid());
  
    E->control.NASSEMBLE = 0;
  
    E->mesh.layer[1] =  E->mesh.layer[2] =  E->mesh.layer[3] = 1.0;
    E->monitor.elapsed_time=0.0;
 
  return;  }


void global_derived_values(E)
     struct All_variables *E;

{
    int d,lx,lz,ly,i,nox,noz,noy;

 if (E->control.NMULTIGRID||E->control.EMULTIGRID)
    { E->mesh.levmax=E->mesh.levels-1;
      E->mesh.nox = E->mesh.mgunitx * (int) pow(2.0,((double)E->mesh.levmax)) + 1;
      E->mesh.noz = E->mesh.mgunitz *(int) pow(2.0,((double)E->mesh.levmax)) + 1; 
      if(E->mesh.nsd == 3 ) 
	  E->mesh.noy = E->mesh.mgunity * (int) pow(2.0,((double)E->mesh.levmax)) + 1; 
    }

 if(E->mesh.nsd != 3) 
   E->mesh.noy = 1;

  E->mesh.nnx[1] = E->mesh.nox;	
  E->mesh.nnx[2] = E->mesh.noz;	
  E->mesh.nnx[3] = E->mesh.noy;	
  E->mesh.elx = E->mesh.nox-1;	
  E->mesh.elz = E->mesh.noz-1;
  E->mesh.ely = max(E->mesh.noy-1,1);

  E->mesh.nel = E->mesh.elx*E->mesh.ely*E->mesh.elz;
  E->mesh.nno = E->mesh.nox*E->mesh.noy*E->mesh.noz;
  E->mesh.nnov = E->mesh.nno;
  E->mesh.neq = E->mesh.nnov*E->mesh.nsd;

  E->mesh.npno = E->mesh.nel;
  E->mesh.nsf = E->mesh.nox*E->mesh.noy;
  for(i=E->mesh.levmax;i>=E->mesh.levmin;i--)  /* set up dimensions for different grids  */
    { if (E->control.NMULTIGRID||E->control.EMULTIGRID)
	{ nox = E->mesh.mgunitx * (int) pow(2.0,(double)i) + 1;
	  noz = E->mesh.mgunitz * (int) pow(2.0,(double)i) + 1;
	  if(E->mesh.nsd==3)
	    noy = E->mesh.mgunity * (int) pow(2.0,(double)i) + 1;
	  else 
	    noy = 1;
	}
    else 
	{ nox = E->mesh.nox;
	  noz = E->mesh.noz;
	  noy = E->mesh.noy;
	}

      E->mesh.ELX[i] = nox-1;
      E->mesh.ELZ[i] = noz-1;
      E->mesh.ELY[i] = max(noy-1,1);
      E->mesh.NNO[i] = nox * noz * noy;
      E->mesh.NEL[i] = (nox-1) * (noz-1) * max((noy-1),1);;
      E->mesh.NPNO[i] = E->mesh.NEL[i] ;
      E->mesh.NOX[i] = nox;
      E->mesh.NOZ[i] = noz;
      E->mesh.NOY[i] = noy;
      E->mesh.NNX[i][1] = nox;	
      E->mesh.NNX[i][2] = noz;	
      E->mesh.NNX[i][3] = noy;	

      E->mesh.NNOV[i] = E->mesh.NNO[i];
      E->mesh.NEQ[i] = E->mesh.nsd * E->mesh.NNOV[i] ;  

    }

    if(E->control.print_convergence)
	fprintf(stderr,"Problem has %d x %d x %d nodes\n",E->mesh.nox,E->mesh.noz,E->mesh.noy);
  	


   return; } 


void read_initial_settings(E)
     struct All_variables *E;
  
{
    void set_convection_defaults();
    void set_2dc_defaults();
    void set_3dc_defaults();
    void set_cg_defaults();
    void set_mg_defaults();
    char logfile[100];
    FILE *fp;
 
  /* first the problem type (defines subsequent behaviour) */

    input_string("Problem",E->control.PROBLEM_TYPE,NULL);
    if ( strcmp(E->control.PROBLEM_TYPE,"convection") == 0)  {
	E->control.CONVECTION = 1; 
	set_convection_defaults(E);
    }

    else if ( strcmp(E->control.PROBLEM_TYPE,"convection-chemical") == 0) {
	E->control.CONVECTION = 1;
	E->control.CHEMISTRY_MODULE=1;
	set_convection_defaults(E);
    }
    
    else {
	fprintf(E->fp,"Unable to determine problem type, assuming convection ... \n");
	E->control.CONVECTION = 1;
	set_convection_defaults(E);
    }
      
  input_string("Geometry",E->control.GEOMETRY,NULL); 
  if ( strcmp(E->control.GEOMETRY,"cart2d") == 0)
    { E->control.CART2D = 1; 
      set_2dc_defaults(E);}
  else if ( strcmp(E->control.GEOMETRY,"axi") == 0)
    { E->control.AXI = 1; 
      }
  else if ( strcmp(E->control.GEOMETRY,"cart2pt5d") == 0)
    { E->control.CART2pt5D = 1; 
      set_2pt5dc_defaults(E);}
  else if ( strcmp(E->control.GEOMETRY,"cart3d") == 0)
    { E->control.CART3D = 1;
      set_3dc_defaults(E);}
  else
    { fprintf(E->fp,"Unable to determine geometry, assuming cartesian 2d ... \n");
      E->control.CART2D = 1; 
      set_2dc_defaults(E); }

  input_string("Solver",E->control.SOLVER_TYPE,NULL);
  if ( strcmp(E->control.SOLVER_TYPE,"cgrad") == 0)
    { E->control.CONJ_GRAD = 1;
      set_cg_defaults(E);}
  else if ( strcmp(E->control.SOLVER_TYPE,"multigrid") == 0)
    { E->control.NMULTIGRID = 1;
      set_mg_defaults(E);}
  else if ( strcmp(E->control.SOLVER_TYPE,"multigrid-el") == 0)
    { E->control.EMULTIGRID = 1;
      set_mg_defaults(E);}
  else
    { fprintf(stderr,"Unable to determine how to solve, specify Solver=VALID_OPTION \n");
      exit(0); 
    }

 
  /* admin */

  input_string("Spacing",E->control.NODE_SPACING,"regular");
  if ( strcmp(E->control.NODE_SPACING,"regular") == 0)
    E->control.GRID_TYPE = 1; 
  else if ( strcmp(E->control.NODE_SPACING,"bound_lyr") == 0)
    E->control.GRID_TYPE = 2;
  else if ( strcmp(E->control.NODE_SPACING,"region") == 0)
    E->control.GRID_TYPE = 3;
  else if ( strcmp(E->control.NODE_SPACING,"ortho_files") == 0)
    E->control.GRID_TYPE = 4;
  else
    {  E->control.GRID_TYPE = 1; }

    /* Information on which files to print, which variables of the flow to calculate and print.
       Default is no information recorded (apart from special things for given applications.
     */
    
    input_string("datatypes",E->control.which_data_files,"");
    input_string("averages",E->control.which_horiz_averages,"");
    input_string("timelog",E->control.which_running_data,"");
    input_string("observables",E->control.which_observable_data,"");

    input_string("datafile",E->control.data_file,"initialize");
    input_string("process_command",E->control.output_written_external_command,"");
    input_boolean("IBM_DX",&(E->control.DX),"off");
    input_boolean("CONMAN",&(E->control.CONMAN),"off");

    /* As early as possible, set up the log file to 
       record information about the progress of the 
       program as it runs 
       */

    sprintf(logfile,"%s/a.log",E->control.data_file, E->control.data_file);
    
    if((fp=fopen(logfile,"w")) == NULL)
	E->fp = stdout;
    else
	E->fp = fp;

    
    if (E->control.NMULTIGRID||E->control.EMULTIGRID) {
	input_int("mgunitx",&(E->mesh.mgunitx),"1");
	input_int("mgunitz",&(E->mesh.mgunitz),"1");
	input_int("mgunity",&(E->mesh.mgunity),"1");
	input_int("levels",&(E->mesh.levels),"0");
    }

    input_boolean("node_assemble",&(E->control.NASSEMBLE),"off");
                                    /* general mesh structure */

    input_boolean("verbose",&(E->control.verbose),"off");
    input_boolean("see_convergence",&(E->control.print_convergence),"off");
    input_boolean("COMPRESS",&(E->control.COMPRESS),"on");
    input_float("sobtol",&(E->control.sob_tolerance),"0.0001");

    input_int("obs_maxlongk",&(E->slice.maxlong),"100,1");
    input_int("obs_minlongk",&(E->slice.minlong),"1,1");

    input_int("stokes_flow_only",&(E->control.stokes),"0");

    input_int("slab_nz",&(E->control.SLAB),"1");
    input_int("plume_nz",&(E->control.PLUME),"1");

        /* for phase change    */

    input_float("Ra_670",&(E->control.Ra_670),"0.0");
    input_float("clapeyron670",&(E->control.clapeyron670),"0.0");
    input_float("transT670",&(E->control.transT670),"0.0");
    input_float("width670",&(E->control.width670),"0.0");

    input_float("Ra_410",&(E->control.Ra_410),"0.0");
    input_float("clapeyron410",&(E->control.clapeyron410),"0.0");
    input_float("transT410",&(E->control.transT410),"0.0");
    input_float("width410",&(E->control.width410),"0.0");


    input_int("restart",&(E->control.restart),"0");
    
    input_int("topvbc",&(E->mesh.topvbc),"0");
    input_int("botvbc",&(E->mesh.botvbc),"0");
    input_int("sidevbc",&(E->mesh.sidevbc),"0");

    input_boolean("periodicx",&(E->mesh.periodic_x),"off");
    input_boolean("periodicy",&(E->mesh.periodic_y),"off");
    input_boolean("depthdominated",&(E->control.depth_dominated),"off");
    input_boolean("eqnzigzag",&(E->control.eqn_zigzag),"off");
    input_boolean("eqnviscosity",&(E->control.eqn_viscosity),"off");

    input_float("topvbxval",&(E->control.VBXtopval),"0.0");
    input_float("botvbxval",&(E->control.VBXbotval),"0.0");
    input_float("topvbyval",&(E->control.VBYtopval),"0.0");
    input_float("botvbyval",&(E->control.VBYbotval),"0.0");
  
    input_int("toptbc",&(E->mesh.toptbc),"1");
    input_int("bottbc",&(E->mesh.bottbc),"1");
    input_float("toptbcval",&(E->control.TBCtopval),"0.0");
    input_float("bottbcval",&(E->control.TBCbotval),"1.0");
 
    input_float("blyr_hwx1",&(E->mesh.bl1width[1]),"nodefault");
    input_float("blyr_hwz1",&(E->mesh.bl1width[2]),"nodefault");
    input_float("blyr_hwy1",&(E->mesh.bl1width[3]),"nodefault");
    input_float("blyr_hwx2",&(E->mesh.bl2width[1]),"nodefault");
    input_float("blyr_hwz2",&(E->mesh.bl2width[2]),"nodefault");
    input_float("blyr_hwy2",&(E->mesh.bl2width[3]),"nodefault");
    input_float("blyr_mgx1",&(E->mesh.bl1mag[1]),"nodefault");
    input_float("blyr_mgz1",&(E->mesh.bl1mag[2]),"nodefault");
    input_float("blyr_mgy1",&(E->mesh.bl1mag[3]),"nodefault");
    input_float("blyr_mgx2",&(E->mesh.bl2mag[1]),"nodefault");
    input_float("blyr_mgz2",&(E->mesh.bl2mag[2]),"nodefault");
    input_float("blyr_mgy2",&(E->mesh.bl2mag[3]),"nodefault");
   
   
    input_float("region_wdx",&(E->mesh.width[1]),"nodefault");
    input_float("region_wdz",&(E->mesh.width[2]),"nodefault");
    input_float("region_wdy",&(E->mesh.width[3]),"nodefault");
    input_float("region_hwx",&(E->mesh.hwidth[1]),"nodefault");
    input_float("region_hwz",&(E->mesh.hwidth[2]),"nodefault");
    input_float("region_hwy",&(E->mesh.hwidth[3]),"nodefault");
    input_float("region_mgx",&(E->mesh.magnitude[1]),"nodefault");
    input_float("region_mgz",&(E->mesh.magnitude[2]),"nodefault");
    input_float("region_mgy",&(E->mesh.magnitude[3]),"nodefault");
    input_float("region_ofx",&(E->mesh.offset[1]),"nodefault");
    input_float("region_ofz",&(E->mesh.offset[2]),"nodefault");
    input_float("region_ofy",&(E->mesh.offset[3]),"nodefault");

    input_string("gridxfile",E->mesh.gridfile[1]," ");
    input_string("gridzfile",E->mesh.gridfile[2]," ");
    input_string("gridyfile",E->mesh.gridfile[3]," ");
    
    input_float("layerd",&(E->data.layer_km),"2800.0");

    E->sphere.ro_dim = E->data.layer_km*1000;
    E->sphere.ro = 1.0;
    E->sphere.ri = 0.0;
        /* for layers    */
    E->viscosity.zlm = 1.0;
    E->viscosity.zlith = 0.0;
    input_int("nz_lmantle",&(E->viscosity.nlm),"1");
    input_int("nz_410",&(E->viscosity.n410),"1");
    input_int("nz_lith",&(E->viscosity.nlith),"1");
    input_int("nz_mid_moho",&(E->viscosity.ncrust2),"1");
    input_float("z_lmantle",&(E->viscosity.zlm),"1.0");
    input_float("z_410",&(E->viscosity.z410),"1.0");
    input_float("z_lith",&(E->viscosity.zlith),"0.0");
    input_float("z_mid_moho",&(E->viscosity.zcrust2),"0.0");

    input_float("z_dd",&(E->viscosity.zdd),"1.0");
    input_int("nz_dd",&(E->viscosity.ndd),"1");
    input_int("nz_comp",&(E->viscosity.ncrust1),"1");
    input_float("z_comp",&(E->viscosity.zcrust1),"0.0");

    input_float("dimenx",&(E->mesh.layer[1]),"nodefault");
    input_float("dimenz",&(E->mesh.layer[2]),"nodefault");
    
    input_int("nodex",&(E->mesh.nox),"nodefault,1,nomax");
    input_int("nodez",&(E->mesh.noz),"nodefault,1,nomax");
    input_int("nodey",&(E->mesh.noy),"1,1,nomax");
    input_boolean("aug_lagr",&(E->control.augmented_Lagr),"off");
    input_double("aug_number",&(E->control.augmented),"0.0");

    input_float("tole_compressibility",&(E->control.tole_comp),"0.0");
    input_boolean("orthogonal",&(E->control.ORTHO),"on");

    input_int("storage_spacing",&(E->control.record_every),"10");
    input_int("storage_always_before",&(E->control.record_all_until),"5");
 
    input_boolean("precond",&(E->control.precondition),"off");
    input_boolean("vprecond",&(E->control.vprecondition),"on");
    input_int("mg_cycle",&(E->control.mg_cycle),"2,0,nomax");
    input_int("down_heavy",&(E->control.down_heavy),"1,0,nomax");
    input_int("up_heavy",&(E->control.up_heavy),"1,0,nomax");
    input_double("accuracy",&(E->control.accuracy),"1.0e-4,0.0,1.0");
    input_int("viterations",&(E->control.max_vel_iterations),"250,0,nomax");

 
    input_int("vhighstep",&(E->control.v_steps_high),"1,0,nomax");
    input_int("vlowstep",&(E->control.v_steps_low),"250,0,nomax");
    input_int("vupperstep",&(E->control.v_steps_upper),"1,0,nomax");
    input_int("piterations",&(E->control.p_iterations),"100,0,nomax");
    input_int("maxsamevisc",&(E->control.max_same_visc),"25,0,nomax");

    input_int("compress_iter_maxstep",&(E->control.compress_iter_maxstep),"100");

    input_float("dissipation_number",&(E->data.disptn_number),"0.0");

  /* data section */ 

  input_float("Ts",&(E->control.Ts),"0.0");
  input_float("ReferenceT",&(E->data.ref_temperature),"2600.0");

//  E->data.ref_temperature = E->data.ref_temperature - E->control.Ts;
  E->control.Ts = E->control.Ts/E->data.ref_temperature;

  E->rad_heat.num = 0;

  input_int("num_radioactives",&(E->rad_heat.num),"0");
  fprintf(stderr,"nn %d\n",E->rad_heat.num);
  fprintf(E->fp,"n_rad %d\n",E->rad_heat.num);

  if (E->rad_heat.num!=0)  {
     input_double("concen_u",&(E->rad_heat.concen_u),"0.0");
     input_double_vector("percent",E->rad_heat.num,(E->rad_heat.percent));
     input_double_vector("heat_g",E->rad_heat.num,(E->rad_heat.heat_g));
     input_double_vector("decay_time",E->rad_heat.num,(E->rad_heat.decay_t));
     E->rad_heat.concen[0]=E->rad_heat.concen_u;
     E->rad_heat.concen[1]=E->rad_heat.concen_u;
     E->rad_heat.concen[2]=E->rad_heat.concen_u*4;
     E->rad_heat.concen[3]=E->rad_heat.concen_u*10000;
     fprintf(E->fp,"Rad_heat %.4e %.4e %.4e %.4e\n",E->rad_heat.percent[0],E->rad_heat.heat_g[0],E->rad_heat.decay_t[0],E->rad_heat.concen[0]);
     fprintf(E->fp,"Rad_heat %.4e %.4e %.4e %.4e\n",E->rad_heat.percent[1],E->rad_heat.heat_g[1],E->rad_heat.decay_t[1],E->rad_heat.concen[1]);
     fprintf(E->fp,"Rad_heat %.4e %.4e %.4e %.4e\n",E->rad_heat.percent[2],E->rad_heat.heat_g[2],E->rad_heat.decay_t[2],E->rad_heat.concen[2]);
     fprintf(E->fp,"Rad_heat %.4e %.4e %.4e %.4e\n",E->rad_heat.percent[3],E->rad_heat.heat_g[3],E->rad_heat.decay_t[3],E->rad_heat.concen[3]);
     fflush(E->fp);
     }
  else
     input_float("Q0",&(E->control.Q0),"0.0");


  input_float("gravacc",&(E->data.grav_acc),"9.81");
  input_float("thermexp",&(E->data.therm_exp),"3.28e-5");
  input_float("thermexp_factor",&(E->data.therm_exp_factor),"3.28e-5");
  input_float("visc_factor",&(E->data.visc_factor),"3.28e-5");
  input_float("cp",&(E->data.Cp),"1200.0");
  input_float("thermdiff",&(E->data.therm_diff),"8.0e-7");
  input_float("thermcond",&(E->data.therm_cond),"3.168");
  input_float("density",&(E->data.density),"3340.0");
  input_float("mdensity",&(E->data.melt_density),"2800.0");
  input_float("wdensity",&(E->data.density_above),"1030.0");
  input_float("rdensity",&(E->data.res_density),"3295.0");
  input_float("heatflux",&(E->data.surf_heat_flux),"4.4e-2");
  input_float("refvisc",&(E->data.ref_viscosity),"nodefault");
  input_float("meltvisc",&(E->data.melt_viscosity),"nodefault");
  input_float("surftemp",&(E->data.surf_temp),"0.0");
  input_float("youngs",&(E->data.youngs_mod),"1.0e11");
  input_float("Te",&(E->data.Te),"0.0");
  input_float("Tsol0",&(E->data.T_sol0),"1373.0");
  input_float("dTsoldz",&(E->data.dTsol_dz),"3.4e-3");
  input_float("dTsoldF",&(E->data.dTsol_dF),"440.0");
  input_float("dTdz",&(E->data.dT_dz),"0.48e-3");
  input_float("deltaS",&(E->data.delta_S),"250.0");
  input_float("gasconst",&(E->data.gas_const),"8.3");     /* not much cause to change these ! */
  input_float("gravconst",&(E->data.grav_const),"6.673e-11");
  input_float("permeability",&(E->data.permeability),"3.0e-10");

  E->monitor.time_scale = E->sphere.ro_dim*E->sphere.ro_dim/(E->data.therm_diff*3600.0*24.0*365.25);   /* years*/

 (E->problem_settings)(E);


return; }

void check_bc_consistency(E)
     struct All_variables *E;

{ int i,lev;

  for(i=1;i<=E->mesh.nno;i++)
    { if ((E->node[i] & VBX) && (E->node[i] & SBX))
	printf("Inconsistent x velocity bc at %d\n",i);
      if ((E->node[i] & VBZ) && (E->node[i] & SBZ))
	printf("Inconsistent z velocity bc at %d\n",i);
      if ((E->node[i] & VBY) && (E->node[i] & SBY))
	printf("Inconsistent y velocity bc at %d\n",i);
      if ((E->node[i] & TBX) && (E->node[i] & FBX))
	printf("Inconsistent x temperature bc at %d\n",i);
      if ((E->node[i] & TBZ) && (E->node[i] & FBZ))
	printf("Inconsistent z temperature bc at %d\n",i);
      if ((E->node[i] & TBY) && (E->node[i] & FBY))
	printf("Inconsistent y temperature bc at %d\n",i); }

  for(lev=E->mesh.levmin;lev<=E->mesh.levmax;lev++)
    for(i=1;i<=E->mesh.NNO[lev];i++)
      { if ((E->NODE[lev][i] & VBX) && (E->NODE[lev][i]  & SBX))
	  printf("Inconsistent x velocity bc at %d,%d\n",lev,i);
	if ((E->NODE[lev][i]  & VBZ) && (E->NODE[lev][i]  & SBZ))
	  printf("Inconsistent z velocity bc at %d,%d\n",lev,i);
	if ((E->NODE[lev][i]  & VBY) && (E->NODE[lev][i]  & SBY))
	  printf("Inconsistent y velocity bc at %d,%d\n",lev,i);
	/* Tbc's not applicable below top level */ }

  return;

}

void set_up_nonmg_aliases(E)
     struct All_variables *E;
     
{ /* Aliases for functions only interested in the highest mg level */
 int i;

  E->eco = E->ECO[E->mesh.levmax]; 
  E->ien = E->IEN[E->mesh.levmax];
  E->id = E->ID[E->mesh.levmax];
  E->lm = E->LMD[E->mesh.levmax];
  E->Vi = E->VI[E->mesh.levmax];
  E->EVi = E->EVI[E->mesh.levmax];
  E->node = E->NODE[E->mesh.levmax];
  E->tw = E->TW[E->mesh.levmax];
  E->Mass = E->MASS[E->mesh.levmax];
  for (i=1;i<=E->mesh.nsd;i++)
    E->X[i] = E->XX[E->mesh.levmax][i];
 
  return; }

report(E,string)
     struct All_variables *E;
     char * string;
{ if(E->control.verbose)
    { fprintf(stderr,"%s\n",string);
      fflush(stderr);
    }
  return;
}

record(E,string)
     struct All_variables *E;
     char * string;
{ if(E->control.verbose)
    { fprintf(E->fp,"%s\n",string);
      fflush(E->fp);
    }

  return;
}



/* =============================================================
   Initialize values which are not problem dependent.
   NOTE: viscosity may be a function of all previous
   input fields (temperature, pressure, velocity, chemistry) and 
   so is always to be done last.
   ============================================================= */


void common_initial_fields(E)
    struct All_variables *E;
{
    void initial_pressure();
    void initial_velocity();
    void read_viscosity_option();

    report(E,"Initialize pressure field");
    initial_pressure(E);
    report(E,"Initialize velocity field");
    initial_velocity(E);
    report(E,"Initialize viscosity field");
    get_viscosity_option(E);

    return;

   }
/* ========================================== */

void initial_pressure(E)
     struct All_variables *E;
{
    int i,node,ii;

    for(i=1;i<=E->mesh.npno;i++)
	    E->P[i]=0.0;

  return; 
}

void initial_velocity(E)
     struct All_variables *E;
{
    int i,node,ii;

    for(i=1;i<=E->mesh.nnov;i++)   {
	E->V[1][i]=0.0;
	E->V[2][i]=0.0;
	}

    if(E->mesh.dof==3)   {
 	for(i=1;i<=E->mesh.nnov;i++)
	    E->V[3][i]=0.0;
        }

    return; 
}
