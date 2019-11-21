
/* ================================================ */

#include <malloc.h>
#include <sys/types.h>
#include <math.h>
#include "element_definitions.h"
#include "global_defs.h"

static float xxsh[5][3] =
        {       { 0.0, 0.0, 0.0 },
                { 0.0,-1.0, -1.0 },
                { 0.0, 1.0, -1.0 },
                { 0.0, 1.0, 1.0 },
                { 0.0,-1.0, 1.0 } };

void Runge_Kutta(E,XMC,XMCpred,C,V,VO)
  struct All_variables *E;
  float *XMC[4],*XMCpred[4],*C,*V[4],*VO[4];
  {

  void get_C_from_markers();
  void velocity_markers();
  void element_markers();

  int i,j;
  float *Vpred[4];

  for(j=1;j<=E->mesh.nsd;j++)   {
    Vpred[j] = (float *) malloc ((E->advection.markers+1)*sizeof(float));
    }

        /*   predicted velocity Vpred at predicted marker positions at t+dt  */
  velocity_markers(E,Vpred,V,XMCpred,E->CElement);
 
         /*   final marker positions at t+dt from modified Euler */
  for (i=1;i<=E->advection.markers;i++)   {
    XMC[1][i] = XMC[1][i] + 0.5*E->advection.timestep*(VO[1][i]+Vpred[1][i]);
    XMC[2][i] = XMC[2][i] + 0.5*E->advection.timestep*(VO[2][i]+Vpred[2][i]);
    XMC[1][i] = min(XMC[1][i],E->XP[1][E->mesh.nox]);
    XMC[1][i] = max(XMC[1][i],E->XP[1][1]);
    XMC[2][i] = min(XMC[2][i],E->XP[2][E->mesh.noz]);
    XMC[2][i] = max(XMC[2][i],E->XP[2][1]);
    }
  
  element_markers(E,XMC,E->CElement);
  get_C_from_markers(E,C,E->CElement);

  for(j=1;j<=E->mesh.nsd;j++)  {  
    free ((void *)Vpred[j]);
    }

  return;
  }

/* ================================================ */

void Euler(E,XMC,XMCpred,C,V,VO)
  struct All_variables *E;
  float *XMC[4],*XMCpred[4],*C,*V[4],*VO[4];
  {

  void get_C_from_markers();
  void velocity_markers();
  void element_markers();

  int i,j;

                 /*   velocity VO at t and x=XMC  */
  velocity_markers(E,VO,V,XMC,E->CElement);
  
                 /*   predicted marker positions at t+dt  */
  for (i=1;i<=E->advection.markers;i++)  {
    XMCpred[1][i] = XMC[1][i] + E->advection.timestep*VO[1][i];
    XMCpred[2][i] = XMC[2][i] + E->advection.timestep*VO[2][i];
    XMCpred[1][i] = min(XMCpred[1][i],E->XP[1][E->mesh.nox]);
    XMCpred[1][i] = max(XMCpred[1][i],E->XP[1][1]);
    XMCpred[2][i] = min(XMCpred[2][i],E->XP[2][E->mesh.noz]);
    XMCpred[2][i] = max(XMCpred[2][i],E->XP[2][1]);
    }

                 /*   predicted compositional field at t+dt  */
  element_markers(E,XMCpred,E->CElement);
  get_C_from_markers(E,C,E->CElement);

  return;
  }

/* ================================================ 
 ================================================  */

 void  get_C_from_markers(E,C,Element)
  struct All_variables *E;
  float *C;
  int *Element;
  {

  int el,i,imark,j,node;
  float C1,temp3,temp1,temp2,temp0;
  static int been_here=0;
  static int *element[3];

  const int elx=E->mesh.elx;
  const int elz=E->mesh.elz;
  const int ely=E->mesh.ely;
  const int nox=E->mesh.nox;
  const int noz=E->mesh.noz;
  const int nno=E->mesh.nno;
  const int nel=E->mesh.nel;
  const int dims=E->mesh.nsd;
  const int ends=enodes[dims];
  const int lev=E->mesh.levmax;

  if (been_here==0)  {
     been_here++;
     element[0] =(int *)malloc((nel+1)*sizeof(int));
     element[1] =(int *)malloc((nel+1)*sizeof(int));
     }

  for (el=1;el<=nel;el++)   {
    element[0][el] = 0;
    element[1][el] = 0;
    }

  
  for (i=1;i<=nno;i++)   {
    C[i] = 0.0;
    }

       /* for each element, count dense and regular marks  */ 
  for (imark=1;imark<=E->advection.markers;imark++)   {
    element[E->C12[imark]][Element[imark]] ++; 
    }

  for (el=1;el<=nel;el++)   {
    temp0 = element[0][el];
    temp1 = element[1][el];
    if (element[0][el] || element[1][el] )
       temp3 = temp1/(temp0+temp1);    /* elemental C */
    else
       temp3 = E->CE[el];    /* elemental C */
    for(j=1;j<=ends;j++) {
       node = E->ien[el].node[j];
       C[node] += E->TWW[lev][el].node[j] * temp3;
       }
    E->CE[el] = temp3;
    }

  for(node=1;node<=nno;node++)  {
     C[node] = C[node]*E->Mass[node];
     }

  return;
  }

/* ================================================ */
 void  element_markers(E,XMC,Element)
  struct All_variables *E;
  float *XMC[4];
  int *Element;
  {
  FILE *fp0;
  char filename1[100];
  int eln,elo,i,j,el,n1,n2,n3,n4;
  float area,XMCold[4],dX[4],weigh1,weigh2,weigh3,weigh4;
  int get_element();

  E->advection.markerIX=1;
  E->advection.markerIZ=1;

  for (i=1;i<=E->advection.markers;i++)  {

    el = get_element(E,XMC[1][i],XMC[2][i],dX);

    Element[i] = el;

    }

  return;
  }

/* ================================================ */
 void  velocity_markers(E,V1,V,XMC,Element)
  struct All_variables *E;
  float *XMC[4],*V1[4],*V[4];
  int *Element;
  {
  FILE *fp0;
  char filename1[100];
  int eln,elo,i,j,el,n1,n2,n3,n4;
  float area,XMCold[4],dX[4],weigh1,weigh2,weigh3,weigh4;
  int get_element();
  static int onf=0;

/*
  sprintf(filename1,"markers%d.%d",E->advection.timesteps,onf);
  fp0=fopen(filename1,"w"); 
  onf=(onf==0)?1:0;
*/

/*  el can also be obtained from CElement[i] and dX can then be
 dX[1]=XMC[1][i]-X[1][ien[el].node[1]]
and  dX[2]=XMC[2][i]-X[2][ien[el].node[1]]. So no element number
is needed to be sought. But since it is so easy to get anyway
for 2D, we do not want to implement this yet
*/

  for (i=1;i<=E->advection.markers;i++)  {

    el = get_element(E,XMC[1][i],XMC[2][i],dX);

    weigh1 = (E->eco[el].size[1]-dX[1])*(E->eco[el].size[2]-dX[2]);
    weigh4 = dX[1]*(E->eco[el].size[2]-dX[2]);
    weigh3 = dX[1]*dX[2];
    weigh2 = (E->eco[el].size[1]-dX[1])*dX[2];
    area = E->eco[el].size[1]*E->eco[el].size[2];

    V1[1][i] = (weigh1*V[1][E->ien[el].node[1]] + weigh2*V[1][E->ien[el].node[2]] 
              + weigh3*V[1][E->ien[el].node[3]] + weigh4*V[1][E->ien[el].node[4]])
               /area;
    V1[2][i] = (weigh1*V[2][E->ien[el].node[1]] + weigh2*V[2][E->ien[el].node[2]] 
              + weigh3*V[2][E->ien[el].node[3]] + weigh4*V[2][E->ien[el].node[4]])
               /area;

    Element[i] = el;

/*
if (i<1200)fprintf(E->fp,"%d %d %g %g %g %g %g %g %g %g\n",i,el,E->eco[el].size[1],E->eco[el].size[2],dX[1],dX[2],XMC[1][i],XMC[2][i],V1[1][i],V1[2][i]);
if (i<1200)fprintf(E->fp,"V1 %d %d %g %g %g %g\n",i,el,V[1][E->ien[el].node[1]],V[1][E->ien[el].node[2]],V[1][E->ien[el].node[3]],V[1][E->ien[el].node[4]]);
if (i<1200)fprintf(E->fp,"V2 %d %d %g %g %g %g\n",i,el,V[2][E->ien[el].node[1]],V[2][E->ien[el].node[2]],V[2][E->ien[el].node[3]],V[2][E->ien[el].node[4]]);
*/

    }

  return;
  }

/* ================================================ */

  int get_element(E,XMC1,XMC2,dX)
  struct All_variables *E;
  float XMC1,XMC2,dX[4];
  {
  int done,i,i1,i2,ii,j,j1,j2,jj,el;
  const int nox = E->mesh.nox;
  const int noz = E->mesh.noz;
  const int elx = E->mesh.elx;
  const int elz = E->mesh.elz;
  static int been_here=0;
  static float dx;

  if (been_here++==0)  {
    dx = E->XP[1][nox]/elx;
    for (i=1;i<=noz;i++)
      fprintf(E->fp,"%g\n",E->XP[2][i]);
    for (i=1;i<=nox;i++)
      fprintf(E->fp,"%g\n",E->XP[1][i]);
    }
  jj=noz+4;

  E->advection.markerIX = min(XMC1/dx+1,elx);
  dX[1] = XMC1-E->XP[1][E->advection.markerIX];

  done = 0;
  i = 1;
  j=1;
  do {
    
    if (XMC2>=E->XP[2][i])  {
       if (XMC2<=E->XP[2][i+1])  {
          E->advection.markerIZ = i;
          dX[2] = XMC2-E->XP[2][i];
          done = 1;
          }
       else  {
          i = min(i+1,elz);
          }
       }
    j++;
    if (j>=jj)  {
       done=1;
       }
    } while (done==0);


  el = E->advection.markerIZ + (E->advection.markerIX-1)*elz;
    if (j>=jj)  {
       fprintf(E->fp,"!!!overflow %g %g %d\n",XMC1,XMC2,el); fflush(E->fp);
       exit(11);
       }

  return (el);
  }

/* ====================================================  */

/* ============================================= */

 float area_of_4node1(x1,y1,x2,y2,x3,y3,x4,y4)
 float x1,y1,x2,y2,x3,y3,x4,y4;

 {
 float temp1,temp2,area;

 temp1 = 0.5*(x1+x2);
 temp2 = 0.5*(y1+y2);

 if (fabs(x1-x2)==2.0)
   area = 2.0*fabs(temp2-y3);
 else if (fabs(y1-y2)==2.0)
   area = 2.0*fabs(temp1-x3);

 return area;
 }

/* ============================================= */

 float area_of_3node(x1,y1,x2,y2,x3,y3)
 float x1,y1,x2,y2,x3,y3;

 {
 float area;

 area = 0.5*max(fabs(y3-y1),fabs(y3-y2))*max(fabs(x3-x1),fabs(x3-x2));

 return area;
 }


/* ============================================= */

 float mean_of_3node(a,x1,y1,x2,y2,x3,y3)
 float x1,y1,x2,y2,x3,y3;
 int a;
 {
 float mean,xm,ym;

 xm = (x1 + x2 + x3)/3.0;
 ym = (y1 + y2 + y3)/3.0;

 mean = 0.25*(1.0+xxsh[a][1]*xm)*(1.0+xxsh[a][2]*ym);

 return mean;
 }

/* ============================================= */

 float mean_of_4node(a,x1,y1,x2,y2,x3,y3,x4,y4)
 float x1,y1,x2,y2,x3,y3,x4,y4;
 int a;
 {
 float mean,xm,ym;

 xm = (x1 + x2 + x3 + x4)*0.25;
 ym = (y1 + y2 + y3 + y4)*0.25;

 mean = 0.25*(1.0+xxsh[a][1]*xm)*(1.0+xxsh[a][2]*ym);

 return mean;
 }

/* ============================================= */

 float mean_of_5node(a,x1,y1,x2,y2,x3,y3,x4,y4,x5,y5)
 float x1,y1,x2,y2,x3,y3,x4,y4,x5,y5;
 int a;
 {

 float mean,xm,ym;

 xm = (x1 + x2 + x3 + x4 + x5)*0.2;
 ym = (y1 + y2 + y3 + y4 + y5)*0.2;

 mean = 0.25*(1.0+xxsh[a][1]*xm)*(1.0+xxsh[a][2]*ym);


 return mean;
 }

/* ================================================ */

 float dist1(XO,XN)
  float XO[4],XN[4];
  {

  float dist2;

  dist2 = sqrt( (XO[1]-XN[1])*(XO[1]-XN[1])
              + (XO[2]-XN[2])*(XO[2]-XN[2]) );

  return (dist2);
  }
