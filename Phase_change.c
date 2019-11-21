#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include <sys/types.h>
#include "element_definitions.h"
#include "global_defs.h"

void phase_change(E,Bb,Bb_b,Bt,Bt_b)
  struct All_variables *E;
  float *Bb,*Bb_b,*Bt,*Bt_b;
{
  static int been_here=0;

  FILE *fp1,*fp2;
  char output_file[255];

  int i,j,k,n,ns;
  float e_pressure,temp;

  static float pt5=0.5;
  static float one=1.0;

  if (been_here++==0)    {

      E->control.Ra_670 = E->control.Ra_670*E->control.Ra_temp
              /(E->data.density*E->data.therm_exp*E->data.ref_temperature);
      E->control.Ra_410 = E->control.Ra_410*E->control.Ra_temp
              /(E->data.density*E->data.therm_exp*E->data.ref_temperature);

      E->control.clapeyron670 = E->control.clapeyron670*E->data.ref_temperature/
                          (E->data.density*E->data.grav_acc*E->sphere.ro_dim);
      E->control.clapeyron410 = E->control.clapeyron410*E->data.ref_temperature/
                          (E->data.density*E->data.grav_acc*E->sphere.ro_dim);

      E->control.width670 = E->sphere.ro_dim/E->control.width670;
      E->control.width410 = E->sphere.ro_dim/E->control.width410;

      E->control.transT670 = E->control.transT670/E->data.ref_temperature;
      E->control.transT410 = E->control.transT410/E->data.ref_temperature;

fprintf(E->fp,"Rab410 670=%g %g Clap410 670=%g %g Di=%g %g \n",E->control.Ra_410,E->control.Ra_670,E->control.clapeyron410,E->control.clapeyron670,E->data.disptn_number,E->data.ref_viscosity);

      }

  for(i=1;i<=E->mesh.nno;i++)  {
    e_pressure = E->viscosity.zlm - E->X[2][i] -
            E->control.clapeyron670*(E->T[i]-E->control.transT670);
    Bb[i] = pt5*(one+tanh(E->control.width670*e_pressure));
    }

  for(i=1;i<=E->mesh.nno;i++)  {
    e_pressure = E->viscosity.z410 - E->X[2][i] -
            E->control.clapeyron410*(E->T[i]-E->control.transT410);
    Bt[i] = pt5*(one+tanh(E->control.width410*e_pressure));
    }

if (E->advection.timesteps%E->control.record_every == 0)   {

    for (j=1;j<=E->mesh.nox;j++)  {
      Bb_b[j]=0.0;
      for (i=1;i<E->mesh.noz;i++)   {
        n = (j-1)*E->mesh.noz + i;
        if (Bb[n]>=pt5&&Bb[n+1]<=pt5)   {
          Bb_b[j]=(E->X[2][n+1]-E->X[2][n])*(pt5-Bb[n])/(Bb[n+1]-Bb[n])+E->X[2][n];
          break;
          }
        }
      }

   for (j=1;j<=E->mesh.nox;j++)  {
      Bt_b[j]=0.0;
      for (i=1;i<E->mesh.noz;i++)   {
        n = (j-1)*E->mesh.noz + i;
        if (Bt[n]>=pt5&&Bt[n+1]<=pt5)  {
          Bt_b[j]=(E->X[2][n+1]-E->X[2][n])*(pt5-Bt[n])/(Bt[n+1]-Bt[n])+E->X[2][n];
          break;
          }
        }
      }

  sprintf(output_file,"%s/fas.%d",E->control.data_file,E->advection.timesteps);
  fp1=fopen(output_file,"w");
  for (j=1;j<=E->mesh.nox;j++)
    fprintf(fp1,"%.4e %.5e %.5e\n",E->X[1][j*E->mesh.noz],Bt_b[j],Bb_b[j]);
  fclose(fp1);

   }


  return;
  }
