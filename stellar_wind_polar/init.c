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
  double Mratio, Lratio, Bcgs, T, mu, a, b, Q, a_eff, M_star, Edd, eta, Rratio;
  double L, c, M_dot, cs, Bq, v_esc, v_inf, vv, beta, M_dot_cgs, v_inf_cgs;
  double x, y, z, xp, yp, zp, r, theta, Rcgs, omega;

  eta = g_inputParam[Eta];
  Rratio = g_inputParam[R_RATIO];
  Mratio = g_inputParam[M_RATIO];
  Lratio = g_inputParam[L_RATIO];
  omega = g_inputParam[OMEGA];
  T = g_inputParam[TT];
  mu = g_inputParam[MU];
  a = g_inputParam[AA];
  b = g_inputParam[Bb];
  Q = g_inputParam[QQ];
  a_eff = g_inputParam[aa_eff];
  M_star = (Mratio*CONST_Msun/UNIT_MASS);
  Edd = (2.6e-5*(Lratio)*(1.0/Mratio));
  L = (Lratio*L_sun/UNIT_L);
  c = 3.0e+5;
/*  M_dot = pow(1.0+a_eff,-(1.0/a_eff)) * a_eff * pow(1.0-a_eff,-1)*
         (L/(c*c))*pow(((Q*Edd)*pow(1.0-Edd,-1)),pow(a_eff,-1)-1.0);
  M_dot = (1.0e-7*CONST_Msun/tyear)*UNIT_TIME/UNIT_MASS;
*/
  M_dot= (L/(c*c))*(a/(1.-a))*pow(Q*Edd/(1.-Edd),((1.-a)/a));
  M_dot=M_dot*pow(1+a,-1./a);
  M_dot_cgs = M_dot*UNIT_MASS/UNIT_TIME;
  cs = sqrt(UNIT_kB*T/(mu*(CONST_AH/UNIT_MASS)*CONST_amu));
  v_esc = sqrt(2.0*UNIT_G*M_star*(1.0-Edd));                              
  v_inf = v_esc * 1.77582;//sqrt((a/(1.0-a)));                                      
  v_inf_cgs = v_inf*UNIT_VELOCITY;
  vv = v_inf*pow(1.0-(1.0/x1),b);                                      
  beta = g_inputParam[BB];

  Bcgs = sqrt(eta*M_dot_cgs*v_inf_cgs/pow(UNIT_LENGTH, 2));
  Bq = Bcgs/UNIT_B;

  g_smallPressure = (v[RHO])*T/(KELVIN*mu); /**< Small value for pressure fix. */
#if EOS == IDEAL
  g_gamma = 1.05;
#endif
#if EOS == ISOTHERMAL                                                  
  g_isoSoundSpeed = sqrt(UNIT_kB*T/(mu*(CONST_AH/UNIT_MASS)*CONST_amu));
#endif
  if(isnan(-vv) || isnan(vv)){vv = 0.00518*v_inf;}                                   
#if ROTATING_FRAME == YES                                                          
  g_OmegaZ = (omega*sqrt((8.0*UNIT_G*M_star)/27.0));                                  
#endif                                                                             
  if(x1 < 1.02 && x2 < CONST_PI/100. && x3 < CONST_PI/100.){
    printf("Bcgs=%e, M_dotcgs=%e, Edd_gam=%e \n",Bcgs,M_dot_cgs/6.35e25,Edd);
  }
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
#if BACKGROUND_FIELD == NO
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
#if BACKGROUND_FIELD == YES
  v[BX1] = v[BX2] = v[BX3] =
  v[AX1] = v[AX2] = v[AX3] = 0.0;
#endif
#endif

}                                                                          

/*================================================================================*/
void Analysis (const Data *d, Grid *grid)
{
}
/*================================================================================*/

/*================================================================================*/
#if BACKGROUND_FIELD == YES
void BackgroundField (double x1, double x2, double x3, double *B0)                                                    
{                                                                                                                    
  double Rratio, Lratio, Mratio;
  double M_dot, v_inf, a, Q, Edd, a_ff, L, c;
  double eta, M_dot_cgs, v_inf_cgs, Rcgs;
  double Bq, Bcgs, beta, r, v_esc, a_eff, M_star;
  double x, y, z;
  double xp, yp, zp;
  double theta;

  beta = g_inputParam[BB];
  eta = g_inputParam[Eta];
  Rratio = g_inputParam[R_RATIO];
  Mratio = g_inputParam[M_RATIO];
  Lratio = g_inputParam[L_RATIO];
  a = g_inputParam[AA];
  Q = g_inputParam[QQ];
  a_eff = g_inputParam[aa_eff];
  Rcgs = Rratio*UNIT_LENGTH;
  Edd = (2.6e-5*(Lratio)*(1.0/Mratio));
  L = (Lratio*L_sun/UNIT_L);
  M_star = (Mratio*CONST_Msun/UNIT_MASS);
  c = 3.0e+5;
  v_esc = sqrt(2.0*UNIT_G*M_star*(1.0-Edd));                              
  v_inf = v_esc * 1.77582;//sqrt((a/(1.0-a)));                                      
  v_inf_cgs = v_inf*UNIT_VELOCITY;

  M_dot= (L/(c*c))*(a/(1.-a))*pow(Q*Edd/(1.-Edd),((1.-a)/a));
  M_dot=M_dot*pow(1+a,-1./a);
  M_dot_cgs = M_dot*UNIT_MASS/UNIT_TIME;

  Bcgs = sqrt(eta*M_dot_cgs*v_inf_cgs/pow(UNIT_LENGTH, 2));
  Bq = Bcgs/UNIT_B;
  //printf("Bcgs=%e, Bq=%e \n",Bcgs, Bq);

  beta *= 0.0174532925;
  x = x1*sin(x2)*cos(x3);
  y = x1*sin(x2)*sin(x3);
  z = x1*cos(x2);
  xp = x*cos(beta) - z*sin(beta);
  yp = y;
  zp = x*sin(beta) + z*cos(beta);
  r = sqrt(xp*xp + yp*yp + zp*zp);
  theta = acos(zp/r);
  EXPAND(B0[0] = Bq*pow(r,-3)*cos(theta);,
         B0[1] = (Bq/2.)*pow(r,-3)*sin(theta);,
         B0[2] = 0.0;)

  //printf("B0=%e, B1=%e, B2=%e \n",B0[0],B0[1],B0[2]);

}
#endif
/*================================================================================*/

void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid) {        
/*================================================================================*/

  int i, j, k, ip, kp, jp, ghost;

  double Cs_p, Mratio, Lratio, T, mu, a, b, Q, a_eff, M_star, Edd, Rratio;
  double L, c, M_dot, ke, Omega2, A, Bcgs, cs, eta, Bq;
  double nu2_c, B, sigma, f, gLx1, gLx2, gcx1, gcx2, gg, beta, Rcgs, vv;
  double x, y, z, xp, yp, zp, r, theta, v_inf, v_esc, v_inf_cgs, M_dot_cgs;
  double vradial, vtheta, vphi;
  double dvdx1, dvdx2, dvdx3;

  double *x1 = grid[IDIR].x;                                                  
  double *x2 = grid[JDIR].x;                                                  
  double *x3 = grid[KDIR].x;
  double *dx1 = grid[IDIR].dx;
  double *dx2 = grid[JDIR].dx;
  double *dx3 = grid[KDIR].dx;
  double ***vx1 = d->Vc[VX1];
  double ***vx2 = d->Vc[VX2];                                                  
  double ***vx3 = d->Vc[VX3];
  double ***rho = d->Vc[RHO];
#if EOS == IDEAL
  double ***prs = d->Vc[PRS];
#endif
#if PHYSICS == MHD
  double ***bx1 = d->Vc[BX1];
  double ***bx2 = d->Vc[BX2];
  double ***bx3 = d->Vc[BX3];
#endif

  eta = g_inputParam[Eta];
  Rratio = g_inputParam[R_RATIO];
  Cs_p = g_inputParam[Cs_P];
  Mratio = g_inputParam[M_RATIO];
  Lratio = g_inputParam[L_RATIO];
  T = g_inputParam[TT];
  mu = g_inputParam[MU];
  a = g_inputParam[AA];
  b = g_inputParam[Bb];
  Q = g_inputParam[QQ];
  a_eff = g_inputParam[aa_eff];
  M_star = (Mratio*CONST_Msun/UNIT_MASS);
  Edd = (2.6e-5*(Lratio)*(1.0/Mratio));
  L = (Lratio*L_sun/UNIT_L);
  c = 3.0e+5;


  M_dot= (L/(c*c))*(a/(1.-a))*pow(Q*Edd/(1.-Edd),((1.-a)/a));
  M_dot=M_dot*pow(1+a,-1./a);
  M_dot_cgs = M_dot*UNIT_MASS/UNIT_TIME;

  ke = ((4.0*CONST_PI*UNIT_G*M_star*c*Edd)/L);
  Omega2 = pow(0.5,2)*(8.0/27.0)*UNIT_G*M_star;
  A = ((1.0/(1.0-a))*((ke*L*Q)/(4.0*CONST_PI*c)));
  cs = sqrt(UNIT_kB*T/(mu*(CONST_AH/UNIT_MASS)*CONST_amu));
  Bq = Bcgs/UNIT_B;
  v_esc = sqrt(2.0*UNIT_G*M_star*(1.0-Edd));                              
  v_inf = v_esc * 1.77582;//sqrt((a/(1.0-a)));                                      
  v_inf_cgs = v_inf*UNIT_VELOCITY;
#if EOS == ISOTHERMAL                                                  
  g_isoSoundSpeed = sqrt(UNIT_kB*T/(mu*(CONST_AH/UNIT_MASS)*CONST_amu));
#endif       
  beta   = 0.0174532925*g_inputParam[BB];
#if EOS == IDEAL
  g_gamma = 1.05;
#endif

  Bcgs = sqrt(eta*M_dot_cgs*v_inf_cgs/pow(UNIT_LENGTH, 2));
  Bq = Bcgs/UNIT_B;
  //printf("eta = %e, M_dot_cgs = %e, v_inf_cgs = %e, Rcgs = %e, Bcgs=%e, Bq=%e \n",eta,M_dot_cgs,v_inf_cgs,Rcgs,Bcgs,Bq);

  ghost = (NX1_TOT - NX1)/2;

  if(side == X1_BEG){                          
    if(box->vpos == CENTER){
      BOX_LOOP(box,k,j,i){ 
  
        rho[k][j][i] = (M_dot/(4.0*CONST_PI*(cs/Cs_p)));          
#if EOS == IDEAL
        prs[k][j][i] = ((rho[k][j][i])*T/(KELVIN*mu));          
#endif
        if (vx1[k][j][ghost] > cs){
          EXPAND(vradial = cs;,
                 vtheta = 0.0;,
                 vphi = 0.0;)
        } else if (vx1[k][j][ghost] < -cs){
          EXPAND(vradial = -cs;,
                 vtheta = 0.0;,
                 vphi = 0.0;)
        } else if (vx1[k][j][ghost] < fabs(cs) || eta < 5.0) {
          EXPAND(vradial = 2.0*vx1[k][j][ghost] - vx1[k][j][ghost+1];,
                 vtheta = 0.0 ;, // 2.0*vx2[k][j][ghost] - vx2[k][j][ghost+1];,
                 vphi = 0.0 ;) //2.0*vx3[k][j][ghost] - vx3[k][j][ghost+1];)
         } else if (vx1[k][j][ghost] < fabs(cs) || eta > 5.0) {
          EXPAND(vradial = 2.0*vx1[k][j][ghost] - vx1[k][j][ghost+1];,
                 vtheta = 2.0*vx2[k][j][ghost] - vx2[k][j][ghost+1];,
                 vphi = 0.0 ;) //2.0*vx3[k][j][ghost] - vx3[k][j][ghost+1];)
        }

        EXPAND(vx1[k][j][i] = vradial;,
               vx2[k][j][i] = vtheta;,
               vx3[k][j][i] = vphi;)


        rho[k][j][i] = (M_dot/(4.0*CONST_PI*(cs/Cs_p)));          
#if EOS == IDEAL
        prs[k][j][i] = ((rho[k][j][i])*T/(KELVIN*mu));          
#endif
/*        EXPAND(vx1[k][j][i] = cs/Cs_p;,
               vx2[k][j][i] = 0.0;,
               vx3[k][j][i] = 0.0;)
*/
        //printf("eta=%e, vradial=%e \n", eta, vradial);

#if PHYSICS == MHD   
#if BACKGROUND_FIELD == NO
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
#if BACKGROUND_FIELD == YES
        EXPAND(bx1[k][j][i] = 0.0;,
               bx1[k][j][i] = 0.0;,
               bx1[k][j][i] = 0.0;)
#endif
#endif                                                            
  }}}
                                                                   
#if EOS == IDEAL
  if(side == 0){DOM_LOOP(k,j,i){
    if (d->Vc[PRS][k][j][i] < (rho[k][j][i])*T/(KELVIN*mu)){
      d->Vc[PRS][k][j][i] = (rho[k][j][i])*T/(KELVIN*mu);
    }
  }}
#endif
}                                                                          
/*================================================================================*/
#if BODY_FORCE != NO
#if CAK == YES
void BodyForceVector(double vm1, double *v, double vp1, double *g, 
                     double xm1, double x1, double xp1, double x2, double x3)
{
/*
  double L, A, ke, a, M_star, gg, Edd, Mratio, Lratio, Q, T;
  double h, c, dvdx1, nu2_c, B, sigma, f, gLx1, mu, temp;
  double Idr, Iv1, Iv2,dxi,dxim1;
  c      = 3.0e+5;
  a      = g_inputParam[AA];
  Q      = g_inputParam[QQ];               
  mu     = g_inputParam[MU];
  Mratio = g_inputParam[M_RATIO];
  M_star = (Mratio*CONST_Msun/UNIT_MASS);
  Lratio = g_inputParam[L_RATIO];
  L      = (Lratio*L_sun/UNIT_L);
  Edd    = (2.6e-5*(Lratio)*(1.0/Mratio));
  T      = g_inputParam[TT];

  gg     = -UNIT_G*M_star*(1.0-Edd)/x1/x1;
  dxi=xp1-x1 ;
  dxim1=x1-xm1;
  dvdx1=-dxi*vm1/(dxim1*(dxi+dxim1)) + (dxi-dxim1)*v[VX1]/(dxi*dxim1) + dxim1*vp1/(dxi*(dxi+dxim1));

  nu2_c = (1.0-(1.0/(x1*x1)));
  ke     = ((4.0*CONST_PI*UNIT_G*M_star*c*Edd)/L);
  B     = ((v[RHO])*Q*c*ke);
  sigma = (x1/fabs(v[VX1]))*(dvdx1)-1.0; 
  f    = ((pow(1.0+sigma,1.0+a)-pow(1.0+sigma*nu2_c,1.0+a))/((1.0+a)*(1.0-nu2_c)*sigma*pow(1.0+sigma,a)));  
  A      = ((1.0/(1.0-a))*((ke*L*Q)/(4.0*CONST_PI*c)));
  gLx1 = (f*A*pow(x1,-2)*pow(dvdx1/B,a));
 // if (dvdx1 < 1.0e-4){gLx1 = 0.0;}

#if EOS == IDEAL
  temp = v[PRS]*KELVIN*mu/v[RHO];
    gLx1 = gLx1*exp(-4.*log(2.)*pow((2.-temp/T-T/temp),2));
#endif

//  if (v[VX1] < 0.0){
    //printf("x1 = %e, dvdx1 = %e, gLx1 = %e, gg = %e, vx1 = %e \n", x1, dvdx1, gLx1, gg, v[VX1]);
//  }

  g[IDIR] = gg + gLx1;
  g[JDIR] = 0.0;
  g[KDIR] = 0.0;
  
*/
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

