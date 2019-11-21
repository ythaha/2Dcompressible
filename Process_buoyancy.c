/*  Here are the routines which process the results of each buoyancy solution, and call
    any relevant output routines. Much of the information has probably been output along
    with the velocity field. (So the velocity vectors and other data are fully in sync).
    However, heat fluxes and temperature averages are calculated here (even when they
    get output the next time around the velocity solver);
    */

#include <math.h>
#include <malloc.h>
#include <sys/types.h>
#include <stdlib.h> /* for "system" command */

#include "element_definitions.h"
#include "global_defs.h"

void process_temp_field(E,ii)
 struct All_variables *E;
    int ii;
{ 
    void heat_flux();
    void plume_heat_flux();
    void output_temp();

    double return_bulk_value();
    void return_horiz_ave();
    void heat_gint_to_nodes();
    void ele_to_nodes();
    static int been_here = 0;
    double *ValueAtNode,*ValueAtNode2;
    float *ValueAtNode_float;
    int i;
    FILE *fp;
    char filename[250];

    
    const int nno=E->mesh.nno;
    const int lev = E->mesh.levmax;

    if ( ((ii % E->control.record_every) == 0))    {
        heat_flux(E);
/*
      output_temp(E,ii);
*/
      }

    ValueAtNode = (double *) malloc((nno+1)*sizeof(double));
    ValueAtNode2 = (double *) malloc((nno+1)*sizeof(double));
    ValueAtNode_float = (float *) malloc((nno+1)*sizeof(float));

    for(i=1;i<=nno;i++) ValueAtNode[i] = E->T[i];
    E->monitor.T_BulkAve = return_bulk_value(E,ValueAtNode,1);
    heat_gint_to_nodes(E,E->heating_visc,ValueAtNode,lev);
//    ele_to_nodes(E,E->heating_visc,ValueAtNode,lev);
    E->monitor.Phi_BulkAve = return_bulk_value(E,ValueAtNode,1);
    heat_gint_to_nodes(E,E->heating_adi,ValueAtNode2,lev);
//    ele_to_nodes(E,E->heating_adi,ValueAtNode2,lev);
    for(i=1;i<=nno;i++) ValueAtNode2[i] *= E->rho_r[i];
    E->monitor.W_BulkAve = return_bulk_value(E,ValueAtNode2,1);
    if(E->monitor.Phi_BulkAve!=0)
	    E->monitor.Err_Phi_W = (E->monitor.Phi_BulkAve - E->monitor.W_BulkAve)/E->monitor.Phi_BulkAve;
    if(E->mesh.bottbc==0 && E->control.TBCbotval!=0)
	    E->monitor.Effi_of_Phi= E->monitor.Phi_BulkAve/E->control.TBCbotval;

    if ((ii%(10*E->control.record_every)) == 0) {
	sprintf(filename,"%s/heating.%d",E->control.data_file,E->monitor.solution_cycles);
	fp=fopen(filename,"w");
	fprintf(fp,"QQ %g %g %g\n",E->control.Ra_temp,E->data.disptn_number,E->rad_heat.total);
	for(i=1;i<=nno;i++)
		fprintf(fp,"%.4e %.4e %.4e\n",ValueAtNode[i],ValueAtNode2[i],E->strain_rate_2[i]);
	fclose(fp);

        //plume_heat_flux(E,ValueAtNode,ValueAtNode2);
    }

    for(i=1;i<=nno;i++) ValueAtNode_float[i] = ValueAtNode[i];
    return_horiz_ave(E,ValueAtNode_float,E->Have.ViscH);
    for(i=1;i<=nno;i++) ValueAtNode_float[i] = ValueAtNode2[i];
    return_horiz_ave(E,ValueAtNode_float,E->Have.AdiH);




    free((void *) ValueAtNode);
    free((void *) ValueAtNode2);
    free((void *) ValueAtNode_float);
    return;
}
/* ===================
    Surface heat flux  
   =================== */

void heat_flux(E)
    struct All_variables *E;
{
    int e,i,j,node,lnode;
    float *mass,*flux,*SU,*RU;
    float VZ[9],u[9],T[9],dTdz[9],area,uT,uT_adv,uT_adv_s,T1[9];

    double xk[3][5];
    struct Shape_function GN;
    struct Shape_function_dA dOmega;
    struct Shape_function_dx GNx;
    void get_global_shape_fn();
    void return_horiz_ave();
    void temperature_boundary_conditions();

    int startstep;
    static int been_here=0;
    static float old_time;

    const int dims=E->mesh.nsd,dofs=E->mesh.dof;
    const int vpts=vpoints[dims];
    const int ppts=ppoints[dims];
    const int ends=enodes[dims];
    const int nno=E->mesh.nno;
    const int lev = E->mesh.levmax;

    
    flux = (float *) malloc((1+nno)*sizeof(float));

    return_horiz_ave(E,E->T,E->Have.T);

    for(i=1;i<=nno;i++)   {
      flux[i] = 0.0;
      E->heatflux[i] = 0.0;
      E->heatflux_adv[i] = 0.0;
      }
    
    for(e=1;e<=E->mesh.nel;e++) {
      get_global_shape_fn(E,e,&GN,&GNx,&dOmega,2,E->mesh.levmax);

        for(j=1;j<=ends;j++)
          VZ[j] = E->V[2][E->ien[e].node[j]];

       

      uT = 0.0;
      uT_adv = 0.0;
      uT_adv_s = 0.0;
      area = 0.0;
      for(i=1;i<=ppts;i++)   {
        u[i] = 0.0;
        T[i] = 0.0;
        T1[i] = 0.0;
        dTdz[i] = 0.0;
        for(j=1;j<=ends;j++)  {
          lnode = (E->ien[e].node[j]-1)%E->mesh.noz+1;
	  if(E->control.compressibility == 0) {
              u[i] += VZ[j]*E->N.ppt[GNPINDEX(j,i)];
	  }
	  else {
	      u[i] += E->rho_r[E->ien[e].node[j]] * VZ[j] * E->N.ppt[GNPINDEX(j,i)];
	  }
          T[i] += E->T[E->ien[e].node[j]]*E->N.ppt[GNPINDEX(j,i)];
          T1[i] += (E->T[E->ien[e].node[j]]-E->Have.T[lnode])*E->N.ppt[GNPINDEX(j,i)];
          dTdz[i] += -E->T[E->ien[e].node[j]]*GNx.ppt[GNPXINDEX(1,j,i)];
          }
        uT = uT + (u[i]*T[i] + dTdz[i])*dOmega.ppt[i];
        uT_adv = uT_adv + u[i]*T1[i]*dOmega.ppt[i];
        uT_adv_s = uT_adv_s + u[i]*fabs(T1[i])*dOmega.ppt[i];
        area += dOmega.ppt[i];
        }

      uT /= area;
      uT_adv /= area;
      uT_adv_s /= area;

      for(j=1;j<=ends;j++)  {
        E->heatflux[E->ien[e].node[j]] += uT*E->TWW[E->mesh.levmax][e].node[j];
        E->heatflux_adv[E->ien[e].node[j]] += uT_adv*E->TWW[E->mesh.levmax][e].node[j];
        flux[E->ien[e].node[j]] += uT_adv_s*E->TWW[E->mesh.levmax][e].node[j];
        }
    }             /* end of e */

    for(i=1;i<=nno;i++)   {
      flux[i] = flux[i]*E->Mass[i];
      E->heatflux[i] = E->heatflux[i]*E->Mass[i];
      E->heatflux_adv[i] = E->heatflux_adv[i]*E->Mass[i];
      }

    for(i=1;i<=E->mesh.nsf;i++)   {
      E->slice.shflux[i] = 2*E->heatflux[E->surf_node[i]]
                           - E->heatflux[E->surf_node[i]-1];

      E->slice.bhflux[i] = 2*E->heatflux[E->surf_node[i]-E->mesh.noz+1]
                           - E->heatflux[E->surf_node[i]-E->mesh.noz+2];
      }

 return_horiz_ave(E,E->heatflux,E->Have.Rho);
 return_horiz_ave(E,E->heatflux_adv,E->Have.F);
 return_horiz_ave(E,flux,E->Have.f);

    for(i=1;i<=nno;i++)  
       E->heatflux[i] = flux[i];

////////for CMB cooling
    /*
    startstep=30000;
    if(been_here==0 && E->monitor.solution_cycles>=startstep) {
	    old_time = E->monitor.elapsed_time;
	    been_here = 1;
    }
    //in Honda_Iwase_1996_EPSL, K=2.0
    if(E->monitor.solution_cycles>=startstep) {
    E->control.TBCbotval -= 2.0 * E->Have.Rho[1] * (E->monitor.elapsed_time - old_time);
    temperature_boundary_conditions(E);
    old_time=E->monitor.elapsed_time;
    }
    */
////////

    
   free((void *)flux);

  return;  
  }
  

void heat_flux_CBF(E,T_old,T,time_interval,diff)
    struct All_variables *E;
    float *T_old;
    float *T;
    float time_interval;
    float diff;
{
    const int dims=E->mesh.nsd;
    const int ends=enodes[dims];
    const int nel=E->mesh.nel;
    const int nno=E->mesh.nno;
    const int elz=E->mesh.elz;

    void get_global_shape_fn();
    void pg_shape_fn();
    void element_residual();

    float *Tdot1;
    int el,node,i,j;
    double Eres[9];   

    struct Shape_function PG;
    struct Shape_function GN;
    struct Shape_function_dA dOmega;
    struct Shape_function_dx GNx;

    Tdot1= (float *)malloc((nno+1)*sizeof(float));
    for(i=1;i<=nno;i++) {
        E->heatflux[i] = 0.0;
        Tdot1[i] = (T[i]-T_old[i])/time_interval;
    }
    for(el=1;el<=E->mesh.nel;el++)    {
        get_global_shape_fn(E,el,&GN,&GNx,&dOmega,0,E->mesh.levmax);
        pg_shape_fn(E,el,&PG,&GNx,E->V,diff);
        element_residual(E,el,PG,GNx,dOmega,E->V,T,Tdot1,Eres,diff);

        for(j=1;j<=ends;j++) {
            node = E->ien[el].node[j];
	    if( (0 == el % elz) || (1 == el % elz) ) {
                E->heatflux[node] += Eres[j];
            }
        }

    } //end el
    for(i=1;i<=nno;i++) {
        E->heatflux[i] *= E->MassForHF[i];
    }

    free((void *) Tdot1 );
	
    for(i=1;i<=E->mesh.nsf;i++) {
      E->slice.shflux[i] = E->heatflux[E->surf_node[i]];
      E->slice.bhflux[i] = E->heatflux[E->surf_node[i]-E->mesh.noz+1];
    }

}


void plume_heat_flux(E,viscH,adiH)
    struct All_variables *E;
    double *viscH,*adiH;
{
    void get_global_shape_fn();

    const int nno=E->mesh.nno;
    const int noz=E->mesh.noz;
    const int nox=E->mesh.nox;

    const int nel=E->mesh.nel;
    const int elz=E->mesh.elz;
    const int elx=E->mesh.elx;

    const int vpts = vpoints[E->mesh.nsd];
    const int ends = enodes[E->mesh.nsd];

    struct Shape_function GN;
    struct Shape_function_dx GNx;
    struct Shape_function_dA dOmega;
    
    int i,j;
    int node,node1,node2,el,element;
    int *tagPlume,*tagPlumeNode;
    float *maxT,*plumeHF,*plumeArea,*maxTadi;
    float *minT;
    float aveT,aveV,aveRho,aveArea;
    float Di,tThreshold;
    float plumeAdiH,downwellingAdiH;
    char filename[250];
    FILE *fp;

    Di=0.5;
    tThreshold = 0.1;


    tagPlume = (int *) malloc((nel+1)*sizeof(int));
    tagPlumeNode = (int *) malloc((nno+1)*sizeof(int));

    maxT = (float *) malloc((elz+1)*sizeof(float));
    maxTadi = (float *) malloc((elz+1)*sizeof(float));
    plumeHF = (float *) malloc((elz+1)*sizeof(float));
    plumeArea = (float *) malloc((elz+1)*sizeof(float));

    minT = (float *) malloc((elz+1)*sizeof(float));

    for(i=1;i<=nel;i++)
	tagPlume[i] = 0;
    for(i=1;i<=nno;i++)
	tagPlumeNode[i] = 0;
    

    //for layer j
    for(j=5;j<=noz-5;j++) {
	maxT[j] = 0.0;
	maxTadi[j] = 0.0;
	plumeHF[j] = 0.0;
	plumeArea[j] = 0.0;
	minT[j] = 1.0;
	//find maxT and minT
	for(i=1;i<=nox;i++) {
	    node = j + (i-1)*noz;
	    if(E->T[node] > maxT[j]) {
		maxT[j]=E->T[node];
	    }
	    if(E->T[node] < minT[j]) {
		minT[j]=E->T[node];
	    }
        }
	//detect plume area,compute plume heat flux
	for(i=1;i<=elx;i++) {
	    element = j + (i-1)*elz;
	    node1 = j + (i-1)*noz;
	    node2 = j + i*noz;
	    aveT = (E->T[node1]+E->T[node2])/2.0;
	    aveV = (E->V[2][node1]+E->V[2][node2])/2.0;
	    aveRho = (E->rho_r[node1]+E->rho_r[node2])/2.0;
	    aveArea = fabs(E->X[1][node2]-E->X[1][node1]);
	    if( (aveT > E->Have.T[j]+tThreshold*(maxT[j]-E->Have.T[j])) 
	        && (aveV > 0) ){
		tagPlumeNode[node1] = tagPlumeNode[node2] = 1;
		tagPlume[element] = 1;
		plumeHF[j] += aveRho*(aveT-E->Have.T[j])*aveV*aveArea;
		plumeArea[j] += aveArea;
	    }
	    else if (aveT < E->Have.T[j] && aveV < 0) {
		tagPlumeNode[node1] = tagPlumeNode[node2] = -1;
		tagPlume[element] = -1;
	    }
	}
    }

    //get adiabatic cooling for plume area
    plumeAdiH=0.0;
    downwellingAdiH=0.0;
    for(el=1;el<nel;el++) {
	get_global_shape_fn(E,el,&GN,&GNx,&dOmega,0,E->mesh.levmax);
	if(1==tagPlume[el]) {
	    for(j=1;j<=vpts;j++)
            for(i=1;i<=ends;i++) {
                node = E->ien[el].node[i];
                plumeAdiH += adiH[node] * E->N.vpt[GNVINDEX(i,j)] * dOmega.vpt[j];
	    }
	}
	if(-1==tagPlume[el]) {
	    for(j=1;j<=vpts;j++)
            for(i=1;i<=ends;i++) {
                node = E->ien[el].node[i];
                downwellingAdiH += adiH[node] * E->N.vpt[GNVINDEX(i,j)] * dOmega.vpt[j];
	    }
	}
    }


    maxTadi[11]=maxT[11];
    for(j=12;j<=noz-11;j++) {
	maxTadi[j] = maxTadi[11]/exp(Di*(E->X[2][j]-E->X[2][11]));
    }
    
    //output
	sprintf(filename,"%s/plumeheatflux.%d",E->control.data_file,E->monitor.solution_cycles);
	fp=fopen(filename,"w");
	for(j=11;j<=noz-11;j++)
		fprintf(fp,"%.4e %.4e %.4e %.4e %.4e %.4e\n",E->X[2][j],plumeHF[j],plumeArea[j],maxT[j],maxTadi[j],minT[j]);
	fclose(fp);
	sprintf(filename,"%s/plumeArea.%d",E->control.data_file,E->monitor.solution_cycles);
	fp=fopen(filename,"w");
	for(i=1;i<=nno;i++)
		fprintf(fp,"%.4e %.4e %d \n",E->X[1][i],E->X[2][i],tagPlumeNode[i]);
	fprintf(fp,"plumeAdiH %g downwellingAdiH %g\n",plumeAdiH,downwellingAdiH);
	fclose(fp);

    free((void *)tagPlume);
    free((void *)tagPlumeNode);
    free((void *)maxT);
    free((void *)maxTadi);
    free((void *)plumeHF);
    free((void *)plumeArea);
    free((void *)minT);
}
