/*================================================================================*/
/*
   Initilisation file for a radiativly driven stellar wind with a non-rigid dipole  
   configuration magnetic field.

   The boundary and initial conditions are taken from Runacres and Owocki (2002)    

   The method for calculating the radiative acceleration comes from CAK (1975)      

   The model only works with polar corrdinates in 2D, with the 
   MHD module. 1D, 3D, HD, RHD, RMHD and other geometries do not 
   work at the moment.

*/
/*================================================================================*/
#include "pluto.h"                                                                  
/*================================================================================*/
/*
  Initialise the Grid acouding to C.A.K. steady state Hydro model.          
*/                                                                                  
void Init (double *v, double x1, double x2, double x3){
/*================================================================================*/
 double Mratio    = g_inputParam[M_RATIO],                                          
        Lratio    = g_inputParam[L_RATIO],                                          
        Bcgs      = g_inputParam[B_CGS],                                            
        T         = g_inputParam[TT],                                               
        mu        = g_inputParam[MU],                                               
        a         = g_inputParam[AA],                                              
        b         = g_inputParam[Bb],                                               
        Q         = g_inputParam[QQ],                                               
        a_eff     = g_inputParam[aa_eff],                                           
        M_star    = (Mratio*CONST_Msun/UNIT_MASS),                                  
        Edd       = (2.6e-5*(Lratio)*(1.0/Mratio)),                                 
        L         = (Lratio*L_sun/UNIT_L),                                          
        c         = 3.0e+5,                                                         
        M_dot     = pow(1.0+a_eff,-(1.0/a_eff)) * a_eff * pow(1.0-a_eff,-1)*
                    (L/(c*c))*pow(((Q*Edd)*pow(1.0-Edd,-1)),pow(a_eff,-1)-1.0),   
        cs        = sqrt(UNIT_kB*T/(mu*(CONST_AH/UNIT_MASS)*CONST_amu)),            
        Bq        = Bcgs/UNIT_B,                                                    
        v_esc     = sqrt(2.0*UNIT_G*M_star*(1.0-Edd)),                              
        v_inf     = v_esc * sqrt((a/(1.0-a))),                                      
        vv        = v_inf*pow(1.0-(1.0/x1),b),                                      
        beta      = g_inputParam[BB],
        x         = 0.0,
        y         = 0.0,
        z         = 0.0,
        xp        = 0.0,
        yp        = 0.0,
        zp        = 0.0,
        r         = 0.0,
        theta     = 0.0;
/*================================================================================*/
 g_smallPressure = (v[RHO])*T/(KELVIN*mu); /**< Small value for pressure fix. */
 #if EOS == IDEAL
 g_gamma = 1.05;
 #endif
 #if EOS == ISOTHERMAL                                                  
 g_isoSoundSpeed = sqrt(UNIT_kB*T/(mu*(CONST_AH/UNIT_MASS)*CONST_amu));
 #endif
 if(isnan(-vv) || isnan(vv)){vv = 0.00518*v_inf;}                                   
 #if ROTATING_FRAME == YES                                                          
 g_OmegaZ = (0.5*sqrt((8.0*UNIT_G*M_star)/27.0));                                  
 #endif                                                                             
 #if EOS == ISOTHERMAL                                                              
 v[RHO] = (M_dot/(4.0*CONST_PI*vv*x1*x1));                                         
 #endif                                                                             
 #if EOS == IDEAL                                                              
 v[RHO] = (M_dot/(4.0*CONST_PI*vv*x1*x1));                                         
 v[PRS] = (v[RHO]*T/(KELVIN*mu));                                                  
 #endif                                                               
 EXPAND(v[VX1] = vv;,                                                 
        v[VX2] = 0.0;,                               
        v[VX3] = 0.0;)                 
 v[TRC] = 0.0;                                                                     
 #if PHYSICS == MHD                                   
 beta *= 0.0174532925;
 x = x1*sin(x2)*cos(x3);
 y = x1*sin(x2)*sin(x3);
 z = x1*cos(x2);
 xp = x*cos(beta) - z*sin(beta);
 yp = y;
 zp = x*sin(beta) + z*cos(beta);
 r = sqrt(xp*xp + yp*yp + zp*zp);
 theta = acos(zp/r);
 EXPAND(v[BX1] = Bq*pow(r,-3)*cos(theta);,            
        v[BX2] = (Bq/2.)*pow(r,-3)*sin(theta);,                 
        v[BX3] = 0.0;)                                                   
 #endif                                                                    
}                                                                          
/*================================================================================*/
void Analysis (const Data *d, Grid *grid){}                               
/*================================================================================*/
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid) {        
/*================================================================================*/
  int i, j, k, ip, kp, jp;

  double Cs_p   = g_inputParam[Cs_P],
         Mratio = g_inputParam[M_RATIO],
         Lratio = g_inputParam[L_RATIO],
         T      = g_inputParam[TT],
         mu     = g_inputParam[MU],
         a      = g_inputParam[AA],
         b      = g_inputParam[Bb],
         Q      = g_inputParam[QQ],
         a_eff  = g_inputParam[aa_eff],
         M_star = (Mratio*CONST_Msun/UNIT_MASS),
         Edd    = (2.6e-5*(Lratio)*(1.0/Mratio)),
         L      = (Lratio*L_sun/UNIT_L),
         c      = 3.0e+5,
         M_dot  = pow(1.0+a_eff,-(1.0/a_eff)) * a_eff * pow(1.0-a_eff,-1)*
                  (L*pow(c,-2))*pow(((Q*Edd)*pow(1.0-Edd,-1)),pow(a_eff,-1)-1.0),
         ke     = ((4.0*CONST_PI*UNIT_G*M_star*c*Edd)/L),
         Omega2 = pow(0.5,2)*(8.0/27.0)*UNIT_G*M_star,
         A      = ((1.0/(1.0-a))*((ke*L*Q)/(4.0*CONST_PI*c))),
         *x1    = grid[IDIR].x,                                                  
         *x2    = grid[JDIR].x,                                                  
         *x3    = grid[KDIR].x,
         *dx1   = grid[IDIR].dx,
         *dx2   = grid[JDIR].dx,
         *dx3   = grid[KDIR].dx,
         ***vx1 = d->Vc[VX1],
         ***vx2 = d->Vc[VX2],                                                  
         ***vx3 = d->Vc[VX3],
         ***rho = d->Vc[RHO],
          ***prs = d->Vc[PRS],
         #if PHYSICS == MHD
         ***bx1 = d->Vc[BX1],
         ***bx2 = d->Vc[BX2],
         ***bx3 = d->Vc[BX3],
         #endif
         Bcgs   = g_inputParam[B_CGS],                                        
         cs     = sqrt(UNIT_kB*T/(mu*(CONST_AH/UNIT_MASS)*CONST_amu)),       
         Bq     = Bcgs/UNIT_B,                                                    
         #if EOS == ISOTHERMAL                                                  
         g_isoSoundSpeed = sqrt(UNIT_kB*T/(mu*(CONST_AH/UNIT_MASS)*CONST_amu)),
          #endif       
         dvdx1 = 0.0,
         dvdx12  = 0.0,
         nu2_c  = 0.0,
         B      = 0.0,
         sigma  = 0.0,
         f      = 0.0,
         gLx1   = 0.0,
         gLx2   = 0.0,
         gcx1   = 0.0,
         gcx2   = 0.0,
         gg     = 0.0,
         beta   = g_inputParam[BB],
         x      = 0.0,
         y      = 0.0,
         z      = 0.0,
         xp     = 0.0,
         yp     = 0.0,
         zp     = 0.0,
         r      = 0.0,
         theta  = 0.0;
/*================================================================================*/
 
#if EOS == IDEAL
  g_gamma = 1.05;
#endif
#if PHYSICS == MHD   
  beta *= 0.0174532925;
#endif

  if(side == X1_BEG){BOX_LOOP(box,k,j,i){                           

    rho[k][j][i] = (M_dot/(4.0*CONST_PI*(cs/Cs_p)));          
    prs[k][j][i] = ((rho[k][j][i])*T/(KELVIN*mu));          

    EXPAND(vx1[k][j][i] = cs/Cs_p;,
           vx2[k][j][i] = 0.0;,
           vx3[k][j][i] = 0.0;)

    //printf("rho=%e, prs=%e \n", rho[k][j][i], prs[k][j][i]);

#if PHYSICS == MHD   
    x = x1[i]*sin(x2[j])*cos(x3[k]);
    y = x1[i]*sin(x2[j])*sin(x3[k]);
    z = x1[i]*cos(x2[j]);
    xp = x*cos(beta) - z*sin(beta);
    yp = y;
    zp = x*sin(beta) + z*cos(beta);
    r = sqrt(xp*xp + yp*yp + zp*zp);
    theta = acos(zp/r);
    EXPAND(bx1[k][j][i] = Bq*pow(r,-3)*cos(theta);,            
           bx2[k][j][i] = (Bq/2.)*pow(r,-3)*sin(theta);,                 
           bx3[k][j][i] = 0.0;)                                                   
#endif                                                            
  }}
                                                                   
/*
  if(side == X2_BEG){X2_BEG_LOOP(k,j,i){                           

   if (k > NX3/2){
     kp = k - NX3_TOT/2 + (NX2_TOT - NX2)/2;
   } 
   else if (k < NX3/2){
     kp = k + NX3_TOT/2 + (NX2_TOT - NX2)/2;
   }
   jp = (NX2_TOT - NX2)/2 + 2 - j;
   ip = i;

   //printf("KEND=%i, NX3_TOT=%i, NX3=%i, k=%i, kp=%i, G_zomes=%i, j=%i, jp=%i \n",KEND,NX3_TOT,NX3,k,kp,NX2_TOT-NX2,j,jp);

 #if EOS == IDEAL                                                    
   rho[k][j][i] = rho[kp][jp][ip];          
   prs[k][j][i] = prs[kp][jp][ip];          
 #endif                                                               
 #if EOS == ISOTHERMAL                                                      
   rho[k][j][i] = rho[kp][jp][ip];           
 #endif                                                                         
   vx1[k][j][i] = vx1[kp][jp][ip];     
   vx2[k][j][i] = vx2[kp][jp][ip];     
   vx3[k][j][i] = vx3[kp][jp][ip];     
 #if PHYSICS == MHD   
   bx1[k][j][i] = bx1[kp][jp][ip];    
   bx2[k][j][i] = bx2[kp][jp][ip];    
   bx3[k][j][i] = bx3[kp][jp][ip];  
 #endif                                                            
 }} 

 if(side == X2_END){X2_END_LOOP(k,j,i){   
   if (k > NX3/2){
     kp = k - NX3_TOT/2 + (NX2_TOT - NX2)/2;
   } 
   else if (k < NX3/2){
     kp = k + NX3_TOT/2 + (NX2_TOT - NX2)/2;
   }
   jp = j - (NX2_TOT - NX2)/2;
   ip = i;

   //printf("KEND=%i, NX3_TOT=%i, NX3=%i, k=%i, kp=%i, G_zomes=%i, j=%i, jp=%i \n",KEND,NX3_TOT,NX3,k,kp,NX2_TOT-NX2,j,jp);

 #if EOS == IDEAL                                                    
   rho[k][j][i] = rho[kp][jp][ip];          
   prs[k][j][i] = prs[kp][jp][ip];          
 #endif                                                               
 #if EOS == ISOTHERMAL                                                      
   rho[k][j][i] = rho[kp][jp][ip];           
 #endif                                                                         
   vx1[k][j][i] = vx1[kp][jp][ip];     
   vx2[k][j][i] = vx2[kp][jp][ip];     
   vx3[k][j][i] = vx3[kp][jp][ip];     
 #if PHYSICS == MHD   
   bx1[k][j][i] = bx1[kp][jp][ip];    
   bx2[k][j][i] = bx2[kp][jp][ip];    
   bx3[k][j][i] = bx3[kp][jp][ip];  
 #endif    
 }} 
*/
 if(side == 0){DOM_LOOP(k,j,i){
   if (d->Vc[PRS][k][j][i] < (rho[k][j][i])*T/(KELVIN*mu)){
/*
    if (isnan(rho[k][j][i])){
      printf("rho[k][j][i] = %e \n", rho[k][j][i]);
    }
*/
    d->Vc[PRS][k][j][i] = (rho[k][j][i])*T/(KELVIN*mu);
   }
 }}

}                                                                          
/*================================================================================*/
#if BODY_FORCE != NO
#if CAK == YES
void BodyForceVector(double v1, double *v, double v3, double *g, 
                     double r1, double x1, double r3, double x2, double x3)
{
  double L, A, ke, a, M_star, gg, Edd, Mratio, Lratio, Q;
  double h, c, dvdx1, nu2_c, B, sigma, f, gLx1, mu, temp;

  c      = 3.0e+5;
  a      = g_inputParam[AA];
  Q      = g_inputParam[QQ];               
  mu     = g_inputParam[MU];
  Mratio = g_inputParam[M_RATIO];
  M_star = (Mratio*CONST_Msun/UNIT_MASS);
  Lratio = g_inputParam[L_RATIO];
  L      = (Lratio*L_sun/UNIT_L),
  Edd    = (2.6e-5*(Lratio)*(1.0/Mratio));

  gg = -UNIT_G*M_star*(1.0-Edd)/x1/x1;

  h = r3 - x1;
  dvdx1 = fabs((v3-v[VX1])/h);
  nu2_c = (1.0-(1.0/(x1*x1)));
  ke     = ((4.0*CONST_PI*UNIT_G*M_star*c*Edd)/L);
  B     = ((v[RHO])*Q*c*ke);
  sigma = (x1/fabs(v[VX1]))*(dvdx1)-1.0; 
  f    = ((pow(1.0+sigma,1.0+a)-pow(1.0+sigma*nu2_c,1.0+a))/((1.0+a)*(1.0-nu2_c)*sigma*pow(1.0+sigma,a)));  
  A      = ((1.0/(1.0-a))*((ke*L*Q)/(4.0*CONST_PI*c)));
  gLx1 = (f*A*pow(x1,-2)*pow(dvdx1/B,a));
  if (dvdx1 < 1.0e-4){gLx1 = 0.0;}

  temp = v[PRS]*KELVIN*mu/v[RHO];

  if (temp > 1.0e+6) {
    gLx1 = 0.0;
  }

  //printf("gLx1 = %e \n", gLx1);

  g[IDIR] = gg + gLx1;
  g[JDIR] = 0.0;
  g[KDIR] = 0.0;

}
#endif
#if CAK == NO
void BodyForceVector(double *v, double *g, double x1, double x2, double x3)
{
  double M_star, gg, Edd, Mratio, Lratio;

  Mratio  = g_inputParam[M_RATIO];
  M_star  = (Mratio*CONST_Msun/UNIT_MASS);
  Lratio  = g_inputParam[L_RATIO];
  Edd     = (2.6e-5*(Lratio)*(1.0/Mratio));

  gg = -UNIT_G*M_star*(1.0-Edd)/x1/x1;

  g[IDIR] = gg;
  g[JDIR] = 0.0;
  g[KDIR] = 0.0;
}
#endif
#endif
/*================================================================================*/

