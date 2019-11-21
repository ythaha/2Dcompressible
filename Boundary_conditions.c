#include "element_definitions.h"
#include "global_defs.h"
#include <math.h>

 
/* ========================================== */

void velocity_boundary_conditions(E)
     struct All_variables *E;
{
  void velocity_refl_vert_bc();
  void velocity_imp_vert_bc();
  void horizontal_bc();
  void velocity_apply_periodic_bcs();
  int lv; 
  int node,d;
 
  for(lv=E->mesh.levmax;lv>=E->mesh.levmin;lv--)  {
    if(E->mesh.botvbc != 1) {
	  horizontal_bc(E,E->VB,1,1,0.0,VBX,0,lv);	 
	  horizontal_bc(E,E->VB,1,2,0.0,VBZ,1,lv);
	  horizontal_bc(E,E->VB,1,1,E->control.VBXbotval,SBX,1,lv);	 
	  horizontal_bc(E,E->VB,1,2,0.0,SBZ,0,lv);
	  if(E->mesh.nsd==3)     {
	    horizontal_bc(E,E->VB,1,3,E->control.VBYbotval,SBY,1,lv);	
	    horizontal_bc(E,E->VB,1,3,0.0,VBY,0,lv);	 
	    }
	  }
    if(E->mesh.topvbc != 1) {
	  horizontal_bc(E,E->VB,E->mesh.NOZ[lv],1,0.0,VBX,0,lv);	 
	  horizontal_bc(E,E->VB,E->mesh.NOZ[lv],2,0.0,VBZ,1,lv);
	  horizontal_bc(E,E->VB,E->mesh.NOZ[lv],1,E->control.VBXtopval,SBX,1,lv); 
	  horizontal_bc(E,E->VB,E->mesh.NOZ[lv],2,0.0,SBZ,0,lv);  
	  if(E->mesh.nsd==3)     {
	    horizontal_bc(E,E->VB,E->mesh.NOZ[lv],3,E->control.VBYtopval,SBY,1,lv); 
	    horizontal_bc(E,E->VB,E->mesh.NOZ[lv],3,0.0,VBY,0,lv);	 
	    } 
	  } 
    }
 
//  velocity_refl_vert_bc(E);				 /* default */

 if(E->mesh.periodic_x || E->mesh.periodic_y)
    velocity_apply_periodic_bcs(E);
 else
    velocity_refl_vert_bc(E);                              /* default */

 for(lv=E->mesh.levmax;lv>=E->mesh.levmin;lv--) {
     if(E->mesh.botvbc == 1) {
	 horizontal_bc(E,E->VB,1,1,E->control.VBXbotval,VBX,1,lv);	 
	 horizontal_bc(E,E->VB,1,2,0.0,VBZ,1,lv);
	 horizontal_bc(E,E->VB,1,1,0.0,SBX,0,lv);	 
	 horizontal_bc(E,E->VB,1,2,0.0,SBZ,0,lv);
	 if(E->mesh.nsd==3) {
	     horizontal_bc(E,E->VB,1,3,E->control.VBYbotval,VBY,1,lv);	 
	     horizontal_bc(E,E->VB,1,3,0.0,SBY,0,lv);	 
	 }
     }
      if(E->mesh.topvbc == 1) {
	  horizontal_bc(E,E->VB,E->mesh.NOZ[lv],1,E->control.VBXtopval,VBX,1,lv);
	  horizontal_bc(E,E->VB,E->mesh.NOZ[lv],2,0.0,VBZ,1,lv);
	  horizontal_bc(E,E->VB,E->mesh.NOZ[lv],1,0.0,SBX,0,lv);	 
	  horizontal_bc(E,E->VB,E->mesh.NOZ[lv],2,0.0,SBZ,0,lv); 
	  if(E->mesh.nsd==3) {
	      horizontal_bc(E,E->VB,E->mesh.NOZ[lv],3,E->control.VBYtopval,VBY,1,lv);
	      horizontal_bc(E,E->VB,E->mesh.NOZ[lv],3,0.0,SBY,0,lv);	 
	  }
      }
 }


if (E->control.verbose)     {
  if(E->mesh.nsd==3)
    for (node=1;node<=E->mesh.nnov;node++)
      fprintf(E->fp,"VB== %d %g %g %g\n",node,E->VB[1][node],E->VB[2][node],E->VB[3][node]);
  else
    for (node=1;node<=E->mesh.nnov;node++)
      fprintf(E->fp,"VB== %d %g %g \n",node,E->VB[1][node],E->VB[2][node]);

  for(lv=E->mesh.levmax;lv>=E->mesh.levmin;lv--)    {
    fprintf(E->fp,"VBB level=%d %d\n",lv,E->mesh.NNO[lv]);
    for (node=1;node<=E->mesh.NNO[lv];node++)
      fprintf(E->fp,"VB== %d %u %u %u\n",node,E->NODE[lv][node]&VBX,E->NODE[lv][node]&VBZ,E->NODE[lv][node]&VBY);
    }
  }

   return;
   }
/* ========================================== */

void composition_boundary_conditions(E)
     struct All_variables *E;
{
  void composition_refl_vert_bc();
  void compositions_conform_bcs();
  void horizontal_bc();
  
  if(E->mesh.botcbc == 1)
    { horizontal_bc(E,E->CB,1,2,E->control.CBCbotval,CBZ,1,E->mesh.levmax);	
      horizontal_bc(E,E->CB,1,2,E->control.CBCbotval,HBZ,0,E->mesh.levmax); }
  else
    { horizontal_bc(E,E->CB,1,2,E->control.CBCbotval,CBZ,0,E->mesh.levmax);	
      horizontal_bc(E,E->CB,1,2,E->control.CBCbotval,HBZ,1,E->mesh.levmax); }
 
  if(E->mesh.topcbc == 1)
    { horizontal_bc(E,E->CB,E->mesh.noz,2,E->control.CBCtopval,CBZ,1,E->mesh.levmax);	
      horizontal_bc(E,E->CB,E->mesh.noz,2,E->control.CBCtopval,HBZ,0,E->mesh.levmax); }
  else
    { horizontal_bc(E,E->CB,E->mesh.noz,2,E->control.CBCtopval,CBZ,0,E->mesh.levmax);	
      horizontal_bc(E,E->CB,E->mesh.noz,2,E->control.CBCtopval,HBZ,1,E->mesh.levmax); }
 
 
  composition_refl_vert_bc(E);				/* default */
 
/* zero flux is used for all the compositional field */

/*
  compositions_conform_bcs(E);  
*/

   return; }

/* ========================================== */

void temperature_boundary_conditions(E)
     struct All_variables *E;
{
  void temperature_refl_vert_bc();
  void temperatures_conform_bcs();
  void horizontal_bc();
  void temperature_apply_periodic_bcs();
  void temperature_imposed_vert_bcs();
  
  if(E->mesh.bottbc == 1)
    { horizontal_bc(E,E->TB,1,2,E->control.TBCbotval,TBZ,1,E->mesh.levmax);	
      horizontal_bc(E,E->TB,1,2,E->control.TBCbotval,FBZ,0,E->mesh.levmax); }
  else
    { horizontal_bc(E,E->TB,1,2,E->control.TBCbotval,TBZ,0,E->mesh.levmax);	
      horizontal_bc(E,E->TB,1,2,E->control.TBCbotval,FBZ,1,E->mesh.levmax); }
 
  if(E->mesh.toptbc == 1)
    { horizontal_bc(E,E->TB,E->mesh.noz,2,E->control.TBCtopval,TBZ,1,E->mesh.levmax);	
      horizontal_bc(E,E->TB,E->mesh.noz,2,E->control.TBCtopval,FBZ,0,E->mesh.levmax); }
  else
    { horizontal_bc(E,E->TB,E->mesh.noz,2,E->control.TBCtopval,TBZ,0,E->mesh.levmax);	
      horizontal_bc(E,E->TB,E->mesh.noz,2,E->control.TBCtopval,FBZ,1,E->mesh.levmax); }
 
 
//  temperature_refl_vert_bc(E);				/* default */
 

  if(E->mesh.periodic_x || E->mesh.periodic_y)
    temperature_apply_periodic_bcs(E);
  else
    temperature_refl_vert_bc(E);                          /* default */

  temperatures_conform_bcs(E);

   return; }

/* ========================================== */

void velocity_refl_vert_bc(E)
     struct All_variables *E;
{
  int i,j,ii,jj;
  int node1,node2;
  int level,nox,noy,noz;
  const int dims=E->mesh.nsd;

  /* for two YOZ planes if 3-D, or two OZ side walls for 2-D */

    for(j=1;j<=E->mesh.noy;j++)
      for(i=1;i<=E->mesh.noz;i++)  {
        node1 = i + (j-1)*E->mesh.noz*E->mesh.nox;
        node2 = node1 + (E->mesh.nox-1)*E->mesh.noz;

           E->VB[1][node1] = 0.0;
           E->VB[1][node2] = 0.0;
           if((i != 1) && (i != E->mesh.noz)) {
              E->VB[2][node1] = 0.0;  
              E->VB[2][node2] = 0.0;
              }
        }      /* end loop for i and j */

  /* for two XOZ planes if 3-D */
	
  if (E->mesh.nsd == 3) {
      for(j=1;j<=E->mesh.nox;j++)
        for(i=1;i<=E->mesh.noz;i++)       {
          node1 = i + (j-1)*E->mesh.noz;
          node2 = node1 +  (E->mesh.noy-1)*E->mesh.noz*E->mesh.nox;
          if (E->mesh.nsd==3) {
              E->VB[3][node2] = 0.0;
              E->VB[3][node1] = 0.0;
          }

          if((i != 1) && (i != E->mesh.noz)){
              E->VB[2][node2] = 0.0;
              E->VB[2][node1] = 0.0;
          }

          }    /* end of loop i & j */

    }           /* end of if */

 
  /* all vbc's apply at all levels  */
  for(level=E->mesh.levmax;level>=E->mesh.levmin;level--) {
    noz = E->mesh.NOZ[level] ;
    noy = E->mesh.NOY[level] ;
    nox = E->mesh.NOX[level] ;

      for(j=1;j<=noy;j++)
        for(i=1;i<=noz;i++) {
          node1 = i + (j-1)*noz*nox;
          node2 = node1 + (nox-1)*noz;
            E->NODE[level][node1] = E->NODE[level][node1] | VBX;
            E->NODE[level][node1] = E->NODE[level][node1] & (~SBX);
	    //
            if((i!=1) && (i!=noz)) {
               E->NODE[level][node1] = E->NODE[level][node1] & (~VBY);
               E->NODE[level][node1] = E->NODE[level][node1] | SBY;
               E->NODE[level][node1] = E->NODE[level][node1] & (~ VBZ);
               E->NODE[level][node1] = E->NODE[level][node1] | SBZ;    
               }
            E->NODE[level][node2] = E->NODE[level][node2] | VBX;
            E->NODE[level][node2] = E->NODE[level][node2] & (~SBX);
            if((i!=1) && (i!=noz)) {
              E->NODE[level][node2] = E->NODE[level][node2] & (~VBY);
              E->NODE[level][node2] = E->NODE[level][node2] | SBY;
              E->NODE[level][node2] = E->NODE[level][node2] & (~ VBZ);
              E->NODE[level][node2] = E->NODE[level][node2] | SBZ;
           	  }
	    //
	    
	    //for no-slip lateral boundary condition
	    /*
            if((i!=1) && (i!=noz)) {
               E->NODE[level][node1] = E->NODE[level][node1] & (~SBY);
               E->NODE[level][node1] = E->NODE[level][node1] | VBY;
               E->NODE[level][node1] = E->NODE[level][node1] & (~ SBZ);
               E->NODE[level][node1] = E->NODE[level][node1] | VBZ;    
               }
            E->NODE[level][node2] = E->NODE[level][node2] | VBX;
            E->NODE[level][node2] = E->NODE[level][node2] & (~SBX);
            if((i!=1) && (i!=noz)) {
              E->NODE[level][node2] = E->NODE[level][node2] & (~SBY);
              E->NODE[level][node2] = E->NODE[level][node2] | VBY;
              E->NODE[level][node2] = E->NODE[level][node2] & (~ SBZ);
              E->NODE[level][node2] = E->NODE[level][node2] | VBZ;
           	  }
	    */
	  }   /* end for loop i & j */
	   

    if (E->mesh.nsd == 3)  {
        for(j=1;j<=nox;j++) 
          for(i=1;i<=noz;i++) {
            node1 = i + (j-1)*noz;

            E->NODE[level][node1] = E->NODE[level][node1] | VBY;
            E->NODE[level][node1] = E->NODE[level][node1] & (~SBY);
            if((i!= 1) && (i != noz))  {
                E->NODE[level][node1] = E->NODE[level][node1] & (~VBZ);
                E->NODE[level][node1] = E->NODE[level][node1] | SBZ;
                } 
            if((j!=1) && (j!=nox) && (i!=1) && (i!=noz)){
                E->NODE[level][node1] = E->NODE[level][node1] & (~VBX);
                E->NODE[level][node1] = E->NODE[level][node1] | SBX;
                }
	        }    /* end for loop i & j  */

        for(j=1;j<=nox;j++)
          for(i=1;i<=noz;i++)       {
            node2 = (noy-1)*noz*nox + i + (j-1)*noz;
            E->NODE[level][node2] = E->NODE[level][node2] | VBY;
            E->NODE[level][node2] = E->NODE[level][node2] & (~SBY);
            if((i!= 1) && (i != noz))  {
                E->NODE[level][node2] = E->NODE[level][node2] & (~VBZ);
                E->NODE[level][node2] = E->NODE[level][node2] | SBZ;
                } 
            if((j!=1) && (j!=nox) && (i!=1) && (i!=noz)){
                E->NODE[level][node2] = E->NODE[level][node2] & (~VBX);
                E->NODE[level][node2] = E->NODE[level][node2] | SBX;
                }
            }
      }               /* end for if dims=3 */
  }                   /* end for loop level */

 
  return;
}

/* =========================================== */
void temperature_refl_vert_bc(E)
     struct All_variables *E;
{
  int i,j;
  int node1,node2;
  const int dims=E->mesh.nsd;

 /* Temps and bc-values  at top level only */

    for(j=1;j<=E->mesh.noy;j++)
      for(i=1;i<=E->mesh.noz;i++) {
        node1 = i + (j-1)*E->mesh.noz*E->mesh.nox;
	node2 = node1 + (E->mesh.nox-1)*E->mesh.noz;
          E->node[node1] = E->node[node1] & (~TBX);
          E->node[node1] = E->node[node1] | FBX;   
          E->TB[1][node1] = 0.0;
          E->node[node2] = E->node[node2] & (~TBX);
          E->node[node2] = E->node[node2] | FBX;
          E->TB[1][node2] = 0.0;
        }       /* end for loop i & j */
	
  if (E->mesh.nsd == 3)  {  
      for(j=1;j<=E->mesh.nox;j++)
        for(i=1;i<=E->mesh.noz;i++) {
          node1 = i + (j-1)*E->mesh.noz;
          E->node[node1] = E->node[node1] & (~TBY);
	      E->node[node1] = E->node[node1] | FBY;
	      E->TB[3][node1] = 0.0;
	      }
      for(j=1;j<=E->mesh.nox;j++)
        for(i=1;i<=E->mesh.noz;i++) {
          node2 = i +(j-1)*E->mesh.noz + (E->mesh.noy-1)*E->mesh.noz*E->mesh.nox;
	  E->node[node2] = E->node[node2] & (~TBY);
	  E->node[node2] = E->node[node2] | FBY;
	  E->TB[3][node2] = 0.0;
          }    /* end loop for i and j */
    }        /* end for if ==3    */


  return;
}

/* =========================================== */
void composition_refl_vert_bc(E)
     struct All_variables *E;
{
  int i,j;
  int node1,node2;
  const int dims=E->mesh.nsd;

 /* Temps and bc-values  at top level only */

    for(j=1;j<=E->mesh.noy;j++)
      for(i=1;i<=E->mesh.noz;i++) {
        node1 = i + (j-1)*E->mesh.noz*E->mesh.nox;
	node2 = node1 + (E->mesh.nox-1)*E->mesh.noz;
          E->node[node1] = E->node[node1] & (~CBX);
          E->node[node1] = E->node[node1] | HBX;   
          E->CB[1][node1] = 0.0;
          E->node[node2] = E->node[node2] & (~CBX);
          E->node[node2] = E->node[node2] | HBX;
          E->CB[1][node2] = 0.0;
        }       /* end for loop i & j */
	
  if (E->mesh.nsd == 3)  {  
      for(j=1;j<=E->mesh.nox;j++)
        for(i=1;i<=E->mesh.noz;i++) {
          node1 = i + (j-1)*E->mesh.noz;
          E->node[node1] = E->node[node1] & (~CBY);
	      E->node[node1] = E->node[node1] | HBY;
	      E->CB[3][node1] = 0.0;
	      }
      for(j=1;j<=E->mesh.nox;j++)
        for(i=1;i<=E->mesh.noz;i++) {
          node2 = i +(j-1)*E->mesh.noz + (E->mesh.noy-1)*E->mesh.noz*E->mesh.nox;
	  E->node[node2] = E->node[node2] & (~CBY);
	  E->node[node2] = E->node[node2] | HBY;
	  E->CB[3][node2] = 0.0;
          }    /* end loop for i and j */
    }        /* end for if ==3    */


  return;
}


/*  =========================================================  */
    

void horizontal_bc(E,BC,ROW,dirn,value,mask,onoff,level)
     struct All_variables *E;
     float *BC[];
     int ROW;
     int dirn;
     float value;
     unsigned int mask;
     char onoff;
     int level;

{
  int i,j,node,rowl;
  const int dims=E->mesh.nsd;

    /* safety feature */
  if(dirn > E->mesh.nsd) 
     return;

  if (ROW==1) 
      rowl = 1;
  else 
      rowl = E->mesh.NOZ[level];
   
  if ( ROW==1 ||
       ROW==E->mesh.NOZ[level] ) {

    /* turn bc marker to zero */
    if (onoff == 0)          {
      for(j=1;j<=E->mesh.NOY[level];j++)
    	for(i=1;i<=E->mesh.NOX[level];i++)     {
    	  node = rowl+(i-1)*E->mesh.NOZ[level]+(j-1)*E->mesh.NOX[level]*E->mesh.NOZ[level];
    	  E->NODE[level][node] = E->NODE[level][node] & (~ mask);
    	  }        /* end for loop i & j */
      }

    /* turn bc marker to one */    
    else        {
      for(j=1;j<=E->mesh.NOY[level];j++)
        for(i=1;i<=E->mesh.NOX[level];i++)       {
    	  node = rowl+(i-1)*E->mesh.NOZ[level]+(j-1)*E->mesh.NOX[level]*E->mesh.NOZ[level];
    	  E->NODE[level][node] = E->NODE[level][node] | (mask);
    	  if(level==E->mesh.levmax)   /* NB */
    	    BC[dirn][node] = value;   
    	  }     /* end for loop i & j */
      }

    }             /* end for if ROW */
    
  return;
}


void velocity_apply_periodic_bcs(E)
    struct All_variables *E;
{
  int i,j,ii,jj;
  int node1,node2;
  int node3,node4;
  int level,nox,noy,noz;
  const int dims=E->mesh.nsd;
  static int been_here = 0;

  if(been_here == 0) {
  
  /* for two YOZ planes if 3-D, or two OZ side walls for 2-D */

    for(j=1;j<=E->mesh.noy;j++)
      for(i=1;i<=E->mesh.noz;i++)  {
        node1 = i + (j-1)*E->mesh.noz*E->mesh.nox;
        node2 = node1 + (E->mesh.nox-1)*E->mesh.noz;

           E->VB[1][node1] = 0.0;
           E->VB[1][node2] = 0.0;
           if((i != 1) && (i != E->mesh.noz)) {
              E->VB[2][node1] = 0.0;  
              E->VB[2][node2] = 0.0;
              }
        }      /* end loop for i and j */

  /* for two XOZ planes if 3-D */
	
  if (E->mesh.nsd == 3) {
      for(j=1;j<=E->mesh.nox;j++)
        for(i=1;i<=E->mesh.noz;i++)       {
          node1 = i + (j-1)*E->mesh.noz;
          node2 = node1 +  (E->mesh.noy-1)*E->mesh.noz*E->mesh.nox;
          if (E->mesh.nsd==3) {
              E->VB[3][node2] = 0.0;
              E->VB[3][node1] = 0.0;
          }

          if((i != 1) && (i != E->mesh.noz)){
              E->VB[2][node2] = 0.0;
              E->VB[2][node1] = 0.0;
          }

          }    /* end of loop i & j */

    }           /* end of if */

 
  /* all vbc's apply at all levels  */
  for(level=E->mesh.levmax;level>=E->mesh.levmin;level--) {
    noz = E->mesh.NOZ[level] ;
    noy = E->mesh.NOY[level] ;
    nox = E->mesh.NOX[level] ;

      for(j=1;j<=noy;j++)
        for(i=1;i<=noz;i++) {
          node1 = i + (j-1)*noz*nox;
          node2 = node1 + (nox-1)*noz;
            E->NODE[level][node1] = E->NODE[level][node1] | VBX;
            E->NODE[level][node1] = E->NODE[level][node1] & (~SBX);
            if((i!=1) && (i!=noz)) {
               E->NODE[level][node1] = E->NODE[level][node1] & (~SBY);
               E->NODE[level][node1] = E->NODE[level][node1] | VBY;
               E->NODE[level][node1] = E->NODE[level][node1] & (~SBZ);
               E->NODE[level][node1] = E->NODE[level][node1] | VBZ;    
               }
            E->NODE[level][node2] = E->NODE[level][node2] | VBX;
            E->NODE[level][node2] = E->NODE[level][node2] & (~SBX);
            if((i!=1) && (i!=noz)) {
              E->NODE[level][node2] = E->NODE[level][node2] & (~SBY);
              E->NODE[level][node2] = E->NODE[level][node2] | VBY;
              E->NODE[level][node2] = E->NODE[level][node2] & (~SBZ);
              E->NODE[level][node2] = E->NODE[level][node2] | VBZ;
           	  }
	  }   /* end for loop i & j */
	   

    if (E->mesh.nsd == 3)  {
        for(j=1;j<=nox;j++) 
          for(i=1;i<=noz;i++) {
            node1 = i + (j-1)*noz;

            E->NODE[level][node1] = E->NODE[level][node1] | VBY;
            E->NODE[level][node1] = E->NODE[level][node1] & (~SBY);
            if((i!= 1) && (i != noz))  {
                E->NODE[level][node1] = E->NODE[level][node1] & (~SBZ);
                E->NODE[level][node1] = E->NODE[level][node1] | VBZ;
                } 
            if((j!=1) && (j!=nox) && (i!=1) && (i!=noz)){
                E->NODE[level][node1] = E->NODE[level][node1] & (~SBX);
                E->NODE[level][node1] = E->NODE[level][node1] | VBX;
                }
	        }    /* end for loop i & j  */

        for(j=1;j<=nox;j++)
          for(i=1;i<=noz;i++)       {
            node2 = (noy-1)*noz*nox + i + (j-1)*noz;
            E->NODE[level][node2] = E->NODE[level][node2] | VBY;
            E->NODE[level][node2] = E->NODE[level][node2] & (~SBY);
            if((i!= 1) && (i != noz))  {
                E->NODE[level][node2] = E->NODE[level][node2] & (~SBZ);
                E->NODE[level][node2] = E->NODE[level][node2] | VBZ;
                } 
            if((j!=1) && (j!=nox) && (i!=1) && (i!=noz)){
                E->NODE[level][node2] = E->NODE[level][node2] & (~SBX);
                E->NODE[level][node2] = E->NODE[level][node2] | VBX;
                }
            }
      }               /* end for if dims=3 */
  }                   /* end for loop level */

  fprintf(E->fp,"Periodic boundary conditions\n");
  been_here ++;
  }

  else {
  /* for two YOZ planes if 3-D, or two OZ side walls for 2-D */
    for(j=1;j<=E->mesh.noy;j++)
      for(i=1;i<=E->mesh.noz;i++)  {
        node1 = i + (j-1)*E->mesh.noz*E->mesh.nox;
        node2 = node1 + (E->mesh.nox-1)*E->mesh.noz;

	node3 = node1 + E->mesh.noz;
	node4 = node2 - E->mesh.noz;

           E->VB[1][node1] = E->V[1][node4];
           E->VB[1][node2] = E->V[1][node3];
           if((i != 1) && (i != E->mesh.noz)) {
              E->VB[2][node1] = E->V[2][node4];  
              E->VB[2][node2] = E->V[2][node3];
              }
        }      /* end loop for i and j */

  /* for two XOZ planes if 3-D */
	
  if (E->mesh.nsd == 3) {
      for(j=1;j<=E->mesh.nox;j++)
        for(i=1;i<=E->mesh.noz;i++)       {
          node1 = i + (j-1)*E->mesh.noz;
          node2 = node1 +  (E->mesh.noy-1)*E->mesh.noz*E->mesh.nox;

	  node3 = node1 + E->mesh.noz*E->mesh.nox;
	  node4 = node2 - E->mesh.noz*E->mesh.nox;
          if (E->mesh.nsd==3) {
              E->VB[3][node2] = E->V[3][node3];
              E->VB[3][node1] = E->V[3][node4];
          }

          if((i != 1) && (i != E->mesh.noz)){
              E->VB[2][node2] = E->V[2][node3];
              E->VB[2][node1] = E->V[2][node4];
          }

          }    /* end of loop i & j */

    }           /* end of if */

  } /* end of else */
  return;
  }

void temperature_apply_periodic_bcs(E)
    struct All_variables *E;
{
  int i,j;
  int node1,node2;
  int node3,node4;
  const int dims=E->mesh.nsd;
  static int been_here = 0;

  if(been_here ==0) {

 /* Temps and bc-values  at top level only */

    for(j=1;j<=E->mesh.noy;j++)
      for(i=1;i<=E->mesh.noz;i++) {
        node1 = i + (j-1)*E->mesh.noz*E->mesh.nox;
	node2 = node1 + (E->mesh.nox-1)*E->mesh.noz;
          E->node[node1] = E->node[node1] & (~FBX);
          E->node[node1] = E->node[node1] | TBX;   
          E->node[node2] = E->node[node2] & (~FBX);
          E->node[node2] = E->node[node2] | TBX;
        }       /* end for loop i & j */
	
  if (E->mesh.nsd == 3)  {  
/*
      for(j=1;j<=E->mesh.nox;j++)
        for(i=1;i<=E->mesh.noz;i++) {
          node1 = i + (j-1)*E->mesh.noz;
          E->node[node1] = E->node[node1] & (~FBY);
	      E->node[node1] = E->node[node1] | TBY;
	      E->TB[3][node1] = 0.0;
	      }
      for(j=1;j<=E->mesh.nox;j++)
        for(i=1;i<=E->mesh.noz;i++) {
          node2 = i +(j-1)*E->mesh.noz + (E->mesh.noy-1)*E->mesh.noz*E->mesh.nox;
	  E->node[node2] = E->node[node2] & (~FBY);
	  E->node[node2] = E->node[node2] | TBY;
	  E->TB[3][node2] = 0.0;
          }   
*/
    }        /* end for if ==3    */

  fprintf(E->fp,"Periodic temperature boundary conditions\n");
  been_here ++;
  }
  else {
 /* Temps and bc-values  at top level only */

    for(j=1;j<=E->mesh.noy;j++)
      for(i=1;i<=E->mesh.noz;i++) {
        node1 = i + (j-1)*E->mesh.noz*E->mesh.nox;
	node2 = node1 + (E->mesh.nox-1)*E->mesh.noz;

	node3 = node1 + E->mesh.noz;
	node4 = node2 - E->mesh.noz;
          E->TB[1][node1] = E->T[node4];
          E->TB[1][node2] = E->T[node3];
        }       /* end for loop i & j */
	
  if (E->mesh.nsd == 3)  {  
    }        /* end for if == 3    */

  } /* end of else */
   
  return;
  }



void strip_bcs_from_residual(E,Res,level)
    struct All_variables *E;
    double *Res;
    int level;
{
    int i;
    const int dims=E->mesh.nsd;
     
    for(i=1;i<=E->mesh.NNO[level];i++) {
        if(E->NODE[level][i] & OFFSIDE)
            continue;
        if ( (E->NODE[level][i] & VBX) != 0 )
            Res[ E->ID[level][i].doff[1] ] = 0.0;
        if ( (E->NODE[level][i] & VBZ) != 0 )
            Res[ E->ID[level][i].doff[2] ] = 0.0;
        if (3==dims && ((E->NODE[level][i] & VBY) != 0))
            Res[ E->ID[level][i].doff[3] ] = 0.0;
    }

    return;
    }


void temperatures_conform_bcs(E)
     struct All_variables *E;
{
    int node;
    unsigned int type;

    for(node=1;node<=E->mesh.nno;node++)  {
	if(E->node[node] & OFFSIDE)
	    continue;
      
	type = (E->node[node] & (TBX | TBZ | TBY));
	
	switch (type) {
	case 0:  /* no match, next node */
	    break;
	case TBX: 
	    E->T[node] = E->TB[1][node];
	    break;
	case TBZ:
	    E->T[node] = E->TB[2][node];
	    break; 
	case TBY:
	    E->T[node] = E->TB[3][node];
	    break; 
	case (TBX | TBZ):     /* clashes ! */
	    E->T[node] = 0.5 * (E->TB[1][node] + E->TB[2][node]);
	    break;
	case (TBX | TBY):     /* clashes ! */
	    E->T[node] = 0.5 * (E->TB[1][node] + E->TB[3][node]);
	    break;
	case (TBZ | TBY):     /* clashes ! */
	    E->T[node] = 0.5 * (E->TB[2][node] + E->TB[3][node]);
	    break;
	case (TBZ | TBY | TBX):     /* clashes ! */
	    E->T[node] = 0.3333333 * (E->TB[1][node] + E->TB[2][node] + E->TB[3][node]);
	    break;
	}
	
	/* next node */ 
    }

return;

 }


void velocities_conform_bcs(E,U)
    struct All_variables *E;
    double *U;
{ 
    int node,d;

    const unsigned int typex = VBX;
    const unsigned int typez = VBZ;
    const unsigned int typey = VBY;
    const int addi_dof = additional_dof[E->mesh.nsd];

    const int dofs = E->mesh.dof;
    const int nno = E->mesh.nno;

    for(node=1;node<=nno;node++) {
	if(E->node[node] & OFFSIDE)
	    continue;

        if (E->node[node] & typex)  
	      U[E->id[node].doff[1]] = E->VB[1][node]; 
	if (E->node[node] & typez)  
	      U[E->id[node].doff[2]] = E->VB[2][node]; 
 	if (3==dofs && E->node[node] & typey)  
	      U[E->id[node].doff[3]] = E->VB[3][node]; 

    } 

    return;
}
