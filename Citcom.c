 /*CITCOM: A finite element convection program written at Caltech 1992 */
 /*Aims to include an iterative matrix solver based on Multigrid techniques */
 /*To do this requires the use of a mixed method and a conjugate-gradient */
 /*approach to determining the */

#include <math.h>
#include <malloc.h>
#include <sys/types.h>

#include "element_definitions.h"
#include "global_defs.h"

extern int Emergency_stop;

main(argc,argv)
     int argc;
     char **argv;
     
{	/* Functions called by main*/
  void general_stokes_solver();
  void read_instructions();
  void solve_constrained_flow();
  void solve_derived_velocities();
  void process_temp_field(); 
  void process_heating();
  void velocity_apply_periodic_bcs();
  void temperature_apply_periodic_bcs();
  

  float dot();

  int k, *temp;
  double CPU_time0(),time,initial_time,start_time;
 
  struct All_variables E;

 E.monitor.solution_cycles=0;

  start_time = time = CPU_time0();
 
  read_instructions(&E,argc,argv);

  E.control.keep_going=1;

     fprintf(stderr,"Input parameters taken from file '%s'\n",argv[1]);
     fprintf(stderr,"Initialization complete after %g seconds\n\n",CPU_time0()-time); fflush(E.fp);
     initial_time = CPU_time0()-time;
     E.monitor.cpu_time_on_vp_it = CPU_time0();

  general_stokes_solver(&E);
  process_new_velocity(&E,E.monitor.solution_cycles);
//  if(E.mesh.periodic_x || E.mesh.periodic_y)
//      velocity_apply_periodic_bcs(&E);

  if (E.control.stokes)  {
     E.control.keep_going=0;
     E.monitor.solution_cycles++; 
     }

  while ( E.control.keep_going   &&  (Emergency_stop == 0) )   {

      E.monitor.solution_cycles++; 
      if(E.monitor.solution_cycles>E.control.print_convergence)
         E.control.print_convergence=1;

      process_heating(&E);

      (E.next_buoyancy_field)(&E);

      process_temp_field(&E,E.monitor.solution_cycles); 

//      if(E.mesh.periodic_x || E.mesh.periodic_y)
//          temperature_apply_periodic_bcs(&E);

      E.monitor.elapsed_time += E.advection.timestep;

      general_stokes_solver(&E);
      process_new_velocity(&E,E.monitor.solution_cycles);
//      if(E.mesh.periodic_x || E.mesh.periodic_y)
//          velocity_apply_periodic_bcs(&E);

      if (E.control.composition && strcmp(E.control.comp_adv_method,"particle")==0) 
          (E.next_buoyancy_field)(&E);     /* correct with R-G */


        fprintf(E.fp,"CPU total = %g & CPU = %g for step %d time = %.4e dt = %.4e  maxT = %.4e sub_iteration%d\n",CPU_time0()-start_time,CPU_time0()-time,E.monitor.solution_cycles,E.monitor.elapsed_time,E.advection.timestep,E.monitor.T_interior,E.advection.last_sub_iterations);
        fflush(E.fp);

        time = CPU_time0();

      }
  
     E.monitor.cpu_time_on_vp_it=CPU_time0()-E.monitor.cpu_time_on_vp_it;
     fprintf(E.fp,"Initialization overhead = %f\n",initial_time);
     fprintf(E.fp,"Average cpu time taken for velocity step = %f\n",
	 E.monitor.cpu_time_on_vp_it/((float)(E.monitor.solution_cycles)));
     fprintf(stderr,"Average cpu time taken for velocity step = %f\n",
	 E.monitor.cpu_time_on_vp_it/((float)(E.monitor.solution_cycles)));

  fclose(E.fp);

  return;  

  } 
