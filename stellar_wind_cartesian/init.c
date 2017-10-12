/*================================================================================*/
/*
   Initilisation file for a radiativly driven stellar wind with a non-rigid dipole  
   configuration magnetic field.

   The boundary and initial conditions are taken from Runacres and Owocki (2002)    

   The method for calculating the radiative acceleration comes from CAK (1975)      

   The model only works with polar corrdinates in 2D & 3D, with the 
   MHD module. 1D, HD, RHD, RMHD and other geometries do not 
   work at the moment.
*/
#include "pluto.h"                                                                  
/*================================================================================*/

void Init (double *v, double x1, double x2, double x3){

  double 
  
  Mratio    = g_inputParam[M_RATIO],                                          
  Lratio    = g_inputParam[L_RATIO],                                          
  Bcgs      = g_inputParam[B_CGS],                                            
  T         = g_inputParam[TT],                                               
  mu        = g_inputParam[MU],                                               
  a         = g_inputParam[AA],                                              
  b         = g_inputParam[b_law],                                               
  Q         = g_inputParam[QQ],                                               
  a_eff     = g_inputParam[aa_eff],                                           
  beta      = g_inputParam[BB],
  M_star    = (Mratio*CONST_Msun/UNIT_MASS),                                  
  Edd       = (2.6e-5*(Lratio)*(1.0/Mratio)),                                 
  L         = (Lratio*L_SUN/UNIT_L),                                          
  c         = 3.0e+5,                                                         
  M_dot     = pow(1.0+a_eff,-(1.0/a_eff)) * a_eff * pow(1.0-a_eff,-1)*
              (L/(c*c))*pow(((Q*Edd)*pow(1.0-Edd,-1)),pow(a_eff,-1)-1.0),     
  cs        = sqrt(UNIT_kB*T/(mu*(CONST_AH/UNIT_MASS)*CONST_amu)),            
  Cs_p      = g_inputParam[Cs_P],
  Bq        = Bcgs/UNIT_B,                                                    
  v_esc     = sqrt(2.0*UNIT_G*M_star*(1.0-Edd)),                              
  v_inf     = v_esc * sqrt((a/(1.0-a))),                                      
  vv        = 0.0,
  Omega     = 0.0,
  x         = 0.0,
  y         = 0.0,
  z         = 0.0,
  xp        = 0.0,
  yp        = 0.0,
  zp        = 0.0,
  bx        = 0.0,
  by        = 0.0,
  bz        = 0.0,
  bxp       = 0.0,
  byp       = 0.0,
  bzp       = 0.0,
  r         = 0.0,
  r2        = 0.0,
  rp        = 0.0,
  rp2       = 0.0,
  rb        = 0.0,
  rb2       = 0.0,
  theta     = 0.0;

  double vel_x1, vel_x2, vel_x3, vel_mag, vel_mag2;

  #if EOS == IDEAL
  g_gamma = 1.05;
  #endif

  #if EOS == ISOTHERMAL                                                  
  g_isoSoundSpeed = sqrt(UNIT_kB*T/(mu*(CONST_AH/UNIT_MASS)*CONST_amu));
  #endif

  r2 = EXPAND(x1*x1,+x2*x2,+x3*x3);
  r  = sqrt(r2);
  beta *= 0.0174532925;

  #if DIMENSIONS == 2
  xp = x1*cos(beta) - x2*sin(beta);
  yp = x1*sin(beta) + x2*cos(beta);
  #endif

  #if DIMENSIONS == 3
  xp = x1*cos(beta) - x3*sin(beta);
  yp = x2;
  zp = x1*sin(beta) + x3*cos(beta);
  #endif

  rp2 = EXPAND(xp*xp, + yp*yp, + zp*zp);
  rp = sqrt(rp2);

  Omega = (0.5*sqrt((8.0*UNIT_G*M_star)/27.0));

  if (r > 1.0) {

    #if DIMENSIONS == 2 
    theta = atan2(xp,yp);
    #endif

    #if DIMENSIONS == 3            
    theta = acos(zp/rp);
    #endif

    vv = v_inf*pow(1.0 - (1.0/r),b);
    v[RHO] = (M_dot/(4.0*CONST_PI*vv*r2));

    #if EOS == IDEAL                                                              
    v[PRS] = (v[RHO]*T/(KELVIN*mu));
    #endif

    D_EXPAND(v[VX1] = vv*x1/r;,                                                 
             v[VX2] = vv*x2/r;,                               
             v[VX3] = vv*x3/r;)

    #if PHYSICS == MHD 
    #if BACKGROUND_FIELD == NO
    #if DIMENSIONS == 2

    bx = 3.0*xp*yp*Bq*pow(rp,-5);
    by = (3.0*pow(yp,2)-pow(rp,2))*Bq*pow(rp,-5);
    bxp = bx*cos(beta) + by*sin(beta);;
    byp = -bx*sin(beta) + by*cos(beta);
    EXPAND(v[BX1] = bxp;,
           v[BX2] = byp;,
           v[BX3] = 0.0;)
    #endif
    #if DIMENSIONS == 3
    bx = 3.0*xp*zp*Bq*pow(rp,-5);
    by = 3.0*yp*zp*Bq*pow(rp,-5);
    bz = (3.0*pow(zp,2)-pow(rp,2))*Bq*pow(rp,-5);
    bxp = bx*cos(beta) + bz*sin(beta);
    byp = by;
    bzp = -bx*sin(beta) + bz*cos(beta);
    EXPAND(v[BX1] = bxp;,
           v[BX2] = byp;,
           v[BX3] = bzp;)
    #endif
    #endif
    #endif

  } else if (r > 0.5 && r <= 1.0) {

    v[RHO] = (M_dot/(4.0*CONST_PI*(cs/Cs_p)));

    #if EOS == IDEAL                                                              
    v[PRS] = (v[RHO]*T/(KELVIN*mu));
    #endif

    vel_mag = v_inf*pow(1.0 - (1.0/1.001),b);

    D_EXPAND(vel_x1 = 0.0;,
             vel_x2 = 0.0;,
             vel_x3 = 0.0;)

    D_EXPAND(v[VX1] = vel_mag*x1/r;,                                                 
             v[VX2] = vel_mag*x2/r;,                               
             v[VX3] = vel_mag*x3/r;)

    #if PHYSICS == MHD 
    #if BACKGROUND_FIELD == NO
    #if DIMENSIONS == 2	

    bx = 3.0*xp*yp*Bq*pow(rp,-5);
    by = (3.0*pow(yp,2)-pow(rp,2))*Bq*pow(rp,-5);
    bxp = bx*cos(beta) + by*sin(beta);;
    byp = -bx*sin(beta) + by*cos(beta);
    EXPAND(v[BX1] = bxp;,
           v[BX2] = byp;,
           v[BX3] = 0.0;)
    #endif
    #if DIMENSIONS == 3
    bx = 3.0*xp*zp*Bq*pow(rp,-5);
    by = 3.0*yp*zp*Bq*pow(rp,-5);
    bz = (3.0*pow(zp,2)-pow(rp,2))*Bq*pow(rp,-5);
    bxp = bx*cos(beta) + bz*sin(beta);
    byp = by;
    bzp = -bx*sin(beta) + bz*cos(beta);
    EXPAND(v[BX1] = bxp;,
           v[BX2] = byp;,
           v[BX3] = bzp;)
    #endif
    #endif
    #endif

  } else if (r <= 0.5 ) {

    v[RHO] = (M_dot/(4.0*CONST_PI*(cs/Cs_p)));

    #if EOS == IDEAL                                                              
    v[PRS] = (v[RHO]*T/(KELVIN*mu));
    #endif

    D_EXPAND(v[VX1] = 0.0;,                                                 
             v[VX2] = 0.0;,                               
             v[VX3] = 0.0;)

    #if PHYSICS == MHD 
    #if BACKGROUND_FIELD == NO
    #if DIMENSIONS == 2
    bx = 0.0;
    by = Bq*16.0;
    bxp = bx*cos(beta) + by*sin(beta);;
    byp = -bx*sin(beta) + by*cos(beta);
    EXPAND(v[BX1] = bxp;,
           v[BX2] = byp;,
           v[BX3] = 0.0;)
    #endif
    #if DIMENSIONS == 3
    bx = 0.0;
    by = 0.0;
    bz = 16.0*Bq;
    bxp = bx*cos(beta) + bz*sin(beta);
    byp = by;
    bzp = -bx*sin(beta) + bz*cos(beta);
    EXPAND(v[BX1] = bxp;,
           v[BX2] = byp;,
           v[BX3] = bzp;)
    #endif 
    #endif
    #endif

  }     

}                                                                          

/*================================================================================*/


/*================================================================================*/
#if BACKGROUND_FIELD == YES
void BackgroundField (double x1, double x2, double x3, double *B0)
{ 

  double x, y, z;
  double xp, yp, zp;
  double bx, by, bz;
  double bxp, byp, bzp;
  double r, r2, rp, rp2, rb, rb2;
  double theta, beta, Bq, Bcgs;

  Bcgs = g_inputParam[B_CGS];                                            
  Bq = Bcgs/UNIT_B;
  beta = g_inputParam[BB];

  r2 = EXPAND(x1*x1,+x2*x2,+x3*x3);
  r  = sqrt(r2);
  beta *= 0.0174532925;

  #if DIMENSIONS == 2
  xp = x1*cos(beta) - x2*sin(beta);
  yp = x1*sin(beta) + x2*cos(beta);
  #endif

  #if DIMENSIONS == 3
  xp = x1*cos(beta) - x3*sin(beta);
  yp = x2;
  zp = x1*sin(beta) + x3*cos(beta);
  #endif

  rp2 = EXPAND(xp*xp, + yp*yp, + zp*zp);
  rp = sqrt(rp2);

  if (r <= 0.5 ) {

    #if DIMENSIONS == 2
    bx = 0.0;
    by = Bq*16.0;
    bxp = bx*cos(beta) + by*sin(beta);;
    byp = -bx*sin(beta) + by*cos(beta);
    B0[0] = bxp;
    B0[1] = byp;
    B0[2] = 0.0;
    #endif
    #if DIMENSIONS == 3
    bx = 0.0;
    by = 0.0;
    bz = 16.0*Bq;
    bxp = bx*cos(beta) + bz*sin(beta);
    byp = by;
    bzp = -bx*sin(beta) + bz*cos(beta);
    B0[0] = bxp;
    B0[1] = byp;
    B0[2] = bzp;
    #endif 

  } else if (r > 0.5){ 

    #if DIMENSIONS == 2	
    bx = 3.0*xp*yp*Bq*pow(rp,-5);
    by = (3.0*pow(yp,2)-pow(rp,2))*Bq*pow(rp,-5);
    bxp = bx*cos(beta) + by*sin(beta);;
    byp = -bx*sin(beta) + by*cos(beta);
    B0[0] = bxp;
    B0[1] = byp;
    B0[2] = 0.0;
    #endif
    #if DIMENSIONS == 3
    bx = 3.0*xp*zp*Bq*pow(rp,-5);
    by = 3.0*yp*zp*Bq*pow(rp,-5);
    bz = (3.0*pow(zp,2)-pow(rp,2))*Bq*pow(rp,-5);
    bxp = bx*cos(beta) + bz*sin(beta);
    byp = by;
    bzp = -bx*sin(beta) + bz*cos(beta);
    B0[0] = bxp;
    B0[1] = byp;
    B0[2] = bzp;
    #endif

  }

}
#endif
/*================================================================================*/




void Analysis (const Data *d, Grid *grid){}                               

/*================================================================================*/

void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid) {        

  int i, j, k, p, o, l, m, n, ghost, phase, kprime, jprime;                                                               

  double vel_x1, vel_x2, vel_x3, vel_mag, vel_mag2;

  double Cs_p, Mratio, Lratio, T, mu, a, b, Q, a_eff, beta, Omega, shell, M_star;
  double Edd, L, c, M_dot, ke, Omega2, A, Bcgs, cs, Bq, v_esc, v_inf;

  double       bx = 0.0,                by = 0.0,                bz = 0.0;
  double      bxp = 0.0,               byp = 0.0,               bzp = 0.0;
  double        x = 0.0,                 y = 0.0,                 z = 0.0;
  double       xp = 0.0,                yp = 0.0,                zp = 0.0;
  double       xa = 0.0,                ya = 0.0,                za = 0.0;
  double       xb = 0.0,                yb = 0.0,                zb = 0.0;
  double      *x1 = grid[IDIR].x,      *x2 = grid[JDIR].x,      *x3 = grid[KDIR].x;
  double     *dx1 = grid[IDIR].dx,    *dx2 = grid[JDIR].dx,    *dx3 = grid[KDIR].dx;
  double   ***vx1 = d->Vc[VX1],     ***vx2 = d->Vc[VX2],     ***vx3 = d->Vc[VX3];
  double      Ix1 = 0.0,               Ix2 = 0.0,               Ix3 = 0.0;
  double      vxI = 0.0,               vyI = 0.0,               vzI = 0.0;

  double V0x = 0.0, V1x = 0.0, V2x = 0.0, V3x = 0.0, V4x = 0.0, V5x = 0.0, V6x = 0.0, V7x = 0.0;
  double V0y = 0.0, V1y = 0.0, V2y = 0.0, V3y = 0.0, V4y = 0.0, V5y = 0.0, V6y = 0.0, V7y = 0.0;
  double V0z = 0.0, V1z = 0.0, V2z = 0.0, V3z = 0.0, V4z = 0.0, V5z = 0.0, V6z = 0.0, V7z = 0.0;
  double N0 = 0.0, N1 = 0.0, N2 = 0.0, N3 = 0.0, N4 = 0.0, N5 = 0.0, N6 = 0.0, N7 = 0.0, Ntot = 0.0;
  double P00 = 0.0, P01 = 0.0, P11 = 0.0, P12 = 0.0, P22 = 0.0; 
  double dvdr = 0.0, dvdr2 = 0.0, Beq = 0.0;
  double vr = 0.0, vr2 = 0.0, r = 0.0, r2 = 0.0, theta = 0.0, phi = 0.0;
  double dr2 = 0.0, dr = 0.0, ddr = 0.0, rb = 0.0, rb2 = 0.0, rp = 0.0, rp2 = 0.0, dI2 = 0.0, dI = 0.0, vrI[2];
  double gL = 0.0, gg = 0.0, gcx1 = 0.0, gcx2 = 0.0, gb = 0.0, nu2_c = 0.0, B = 0.0, sigma = 0.0, f = 0.0, vv = 0.0;

  Cs_p    = g_inputParam[Cs_P];
  Mratio  = g_inputParam[M_RATIO];
  Lratio  = g_inputParam[L_RATIO];
  T       = g_inputParam[TT];
  mu      = g_inputParam[MU];
  a       = g_inputParam[AA];
  b       = g_inputParam[b_law];
  Q       = g_inputParam[QQ];
  a_eff   = g_inputParam[aa_eff];
  beta    = g_inputParam[BB];
  Omega   = g_inputParam[OMEGA];
  shell   = g_inputParam[SHELL];
  M_star  = (Mratio*CONST_Msun/UNIT_MASS);
  Edd     = (2.6e-5*(Lratio)*(1.0/Mratio));
  L       = (Lratio*L_SUN/UNIT_L);
  c       = 3.0e+5;
  M_dot   = pow(1.0+a_eff,-(1.0/a_eff)) * a_eff * pow(1.0-a_eff,-1)*
            (L*pow(c,-2))*pow(((Q*Edd)*pow(1.0-Edd,-1)),pow(a_eff,-1)-1.0);
  ke      = ((4.0*CONST_PI*UNIT_G*M_star*c*Edd)/L);
  Omega2  = pow(Omega,2)*(8.0/27.0)*UNIT_G*M_star;
  A       = ((1.0/(1.0-a))*((ke*L*Q)/(4.0*CONST_PI*c)));
  Bcgs    = g_inputParam[B_CGS];        
  cs      = sqrt(UNIT_kB*T/(mu*(CONST_AH/UNIT_MASS)*CONST_amu));
  Bq      = Bcgs/UNIT_B;
  v_esc   = sqrt(2.0*UNIT_G*M_star*(1.0-Edd));
  v_inf   = v_esc * sqrt((a/(1.0-a)));

  #if EOS == IDEAL
  g_gamma = 1.05;
  #endif
  #if EOS == ISOTHERMAL                                                  
  g_isoSoundSpeed = sqrt(UNIT_kB*T/(mu*(CONST_AH/UNIT_MASS)*CONST_amu));
  #endif
  beta *= 0.0174532925;

  if(side == 0){DOM_LOOP(k,j,i){
    r2  = EXPAND(x1[i]*x1[i],+x2[j]*x2[j],+x3[k]*x3[k]);
    r   = sqrt(r2);

    #if DIMENSIONS == 2
    xp = x1[i]*cos(beta) - x2[j]*sin(beta);
    yp = x1[i]*sin(beta) + x2[j]*cos(beta);
    #endif

    #if DIMENSIONS == 3
    xp = x1[i]*cos(beta) - x3[k]*sin(beta);
    yp = x2[j];
    zp = x1[i]*sin(beta) + x3[k]*cos(beta);
    #endif

    rp2 = EXPAND(xp*xp, + yp*yp, + zp*zp);
    rp = sqrt(xp*xp + yp*yp + zp*zp);

    if (r <= 0.5 ) {

      d->Vc[RHO][k][j][i] = (M_dot/(4.0*CONST_PI*(cs/Cs_p)));

      #if EOS == IDEAL                                                              
      d->Vc[PRS][k][j][i] = (d->Vc[RHO][k][j][i]*T/(KELVIN*mu));
      #endif

      D_EXPAND(d->Vc[VX1][k][j][i] = 0.0;,                                                 
               d->Vc[VX2][k][j][i] = 0.0;,                               
               d->Vc[VX3][k][j][i] = 0.0;)

      #if PHYSICS == MHD 
      #if BACKGROUND_FIELD == NO
      #if DIMENSIONS == 2
      bx = 0.0;
      by = 16.0*Bq;
      bxp = bx*cos(beta) + by*sin(beta);
      byp = -bx*sin(beta) + by*cos(beta);
      EXPAND(d->Vc[BX1][k][j][i] = bxp;,
             d->Vc[BX2][k][j][i] = byp;,
             d->Vc[BX3][k][j][i] = 0.0;)
      #endif
      #if DIMENSIONS == 3
      bx = 0.0;
      by = 0.0;
      bz = 16.0*Bq;
      bxp = bx*cos(beta) + bz*sin(beta);
      byp = by;
      bzp = -bx*sin(beta) + bz*cos(beta);
      EXPAND(d->Vc[BX1][k][j][i] = bxp;,
             d->Vc[BX2][k][j][i] = byp;,
             d->Vc[BX3][k][j][i] = bzp;)
      #endif 
      #endif
      #endif

      d->flag[k][j][i] |= FLAG_INTERNAL_BOUNDARY;

    } else if (r > 0.5 && r <= 1.0) {
      
      d->Vc[RHO][k][j][i] = (M_dot/(4.0*CONST_PI*(cs/Cs_p)));

      #if EOS == IDEAL                                                              
      d->Vc[PRS][k][j][i] = (d->Vc[RHO][k][j][i]*T/(KELVIN*mu));
      #endif

      if(fabs(x1[i]) < 1.0 && fabs(x1[i+1]) > 1.0){
      
      }

      if(fabs(x2[j]) < 1.0 && fabs(x2[j+1]) < 1.0){
        
      }

      vel_mag2 = D_EXPAND(d->Vc[VX1][k][j][i]*d->Vc[VX1][k][j][i], + 
                          d->Vc[VX2][k][j][i]*d->Vc[VX2][k][j][i], +
                          d->Vc[VX3][k][j][i]*d->Vc[VX3][k][j][i]);
      vel_mag = sqrt(vel_mag2);

      D_EXPAND(vel_x1 = 0.0;,
               vel_x2 = 0.0;,
               vel_x3 = 0.0;)

      D_EXPAND(d->Vc[VX1][k][j][i] = vel_mag*x1[i]/r;,                                                 
               d->Vc[VX2][k][j][i] = vel_mag*x2[j]/r;,                               
               d->Vc[VX3][k][j][i] = vel_mag*x3[k]/r;)

      #if PHYSICS == MHD 
      #if BACKGROUND_FIELD == NO
      #if DIMENSIONS == 2	
      bx = 3.0*xp*yp*Bq*pow(rp,-5);
      by = (3.0*pow(yp,2)-pow(rp,2))*Bq*pow(rp,-5);
      bxp = bx*cos(beta) + by*sin(beta);
      byp = -bx*sin(beta) + by*cos(beta);
      EXPAND(d->Vc[BX1][k][j][i] = bxp;,
             d->Vc[BX2][k][j][i] = byp;,
             d->Vc[BX3][k][j][i] = 0.0;)
      #endif
      #if DIMENSIONS == 3
      bx = 3.0*xp*zp*Bq*pow(rp,-5);
      by = 3.0*yp*zp*Bq*pow(rp,-5);
      bz = (3.0*pow(zp,2)-pow(rp,2))*Bq*pow(rp,-5);
      bxp = bx*cos(beta) + bz*sin(beta);
      byp = by;
      bzp = -bx*sin(beta) + bz*cos(beta);
      EXPAND(d->Vc[BX1][k][j][i] = bxp;,
             d->Vc[BX2][k][j][i] = byp;,
             d->Vc[BX3][k][j][i] = bzp;)
      #endif
      #endif
      #endif

      d->flag[k][j][i] |= FLAG_INTERNAL_BOUNDARY;

    }

    /* - Radiative diriving calculation. - */
    if (r > 1.0){

      /* - Determine infinitesimal radial distance dr - */
      Ntot = EXPAND(dx1[i],*dx2[j],*dx3[k]);
      dr2  = EXPAND(dx1[i]*dx1[i],+dx2[j]*dx2[j],+dx3[k]*dx3[k]);
      dr   = sqrt(dr2);

      for (p = 0; p < 2; p++) {

        if (p == 0) {
          ddr = -0.4*dr;
        } else {
          ddr = 0.4*dr;
        }

        /* - Find nearest nighbours and record coordinates - */
        #if DIMENSIONS == 2
        theta = atan2(x1[i],x2[j]);
        Ix1 = (r+ddr)*sin(theta);
        Ix2 = (r+ddr)*cos(theta);

        if (Ix1>x1[i] && Ix2>x2[j]){
          N0 = (x1[i+1]-Ix1)*(x2[j+1]-Ix2);///Ntot;
          N1 = (Ix1-x1[i]  )*(x2[j+1]-Ix2);///Ntot;
          N2 = (x1[i+1]-Ix1)*(Ix2-x2[j]  );///Ntot;
          N3 = (Ix1-x1[i]  )*(Ix2-x2[j]  );///Ntot;
          Ntot = N0 + N1 + N2 + N3;
          N0 /= Ntot;
          N1 /= Ntot;
          N2 /= Ntot;
          N3 /= Ntot;
          vxI = (vx1[k][j][i]*N0+vx1[k][j][i+1]*N1+vx1[k][j+1][i]*N2+vx1[k][j+1][i+1]*N3);
          vyI = (vx2[k][j][i]*N0+vx2[k][j][i+1]*N1+vx2[k][j+1][i]*N2+vx2[k][j+1][i+1]*N3);
          vrI[p] = vxI*x1[i]/r + vyI*x2[j]/r;
        } else if (Ix1<x1[i] && Ix2>x2[j]){
          N0 = (x1[i]-Ix1)*(x2[j+1]-Ix2);///Ntot;
          N1 = (Ix1-x1[i-1])*(x2[j+1]-Ix2);///Ntot;
          N2 = (x1[i]-Ix1)*(Ix2-x2[j]);///Ntot;
          N3 = (Ix1-x1[i-1])*(Ix2-x2[j]);///Ntot;
          Ntot = N0 + N1 + N2 + N3;
          N0 /= Ntot;
          N1 /= Ntot;
          N2 /= Ntot;
          N3 /= Ntot;
          vxI = (vx1[k][j][i-1]*N0+vx1[k][j][i]*N1+vx1[k][j+1][i-1]*N2+vx1[k][j+1][i]*N3);
          vyI = (vx2[k][j][i-1]*N0+vx2[k][j][i]*N1+vx2[k][j+1][i-1]*N2+vx2[k][j+1][i]*N3);
          vrI[p] = vxI*x1[i]/r + vyI*x2[j]/r;
        } else if (Ix1<x1[i] && Ix2<x2[j]){
          N0 = (x1[i]-Ix1)*(x2[j]-Ix2);///Ntot;
          N1 = (Ix1-x1[i-1])*(x2[j]-Ix2);///Ntot;
          N2 = (x1[i]-Ix1)*(Ix2-x2[j-1]);///Ntot;
          N3 = (Ix1-x1[i-1])*(Ix2-x2[j-1]);///Ntot;
          Ntot = N0 + N1 + N2 + N3;
          N0 /= Ntot;
          N1 /= Ntot;
          N2 /= Ntot;
          N3 /= Ntot;
          vxI = (vx1[k][j-1][i-1]*N0+vx1[k][j-1][i]*N1+vx1[k][j][i-1]*N2+vx1[k][j][i]*N3);
          vyI = (vx2[k][j-1][i-1]*N0+vx2[k][j-1][i]*N1+vx2[k][j][i-1]*N2+vx2[k][j][i]*N3);
          vrI[p] = vxI*x1[i]/r + vyI*x2[j]/r;
        } else if (Ix1>x1[i] && Ix2<x2[j]){
          N0 = (x1[i+1]-Ix1)*(x2[j]-Ix2);///Ntot;
          N1 = (Ix1-x1[i])*(x2[j]-Ix2);///Ntot;
          N2 = (x1[i+1]-Ix1)*(Ix2-x2[j-1]);///Ntot;
          N3 = (Ix1-x1[i])*(Ix2-x2[j-1]);///Ntot;
          Ntot = N0 + N1 + N2 + N3;
          N0 /= Ntot;
          N1 /= Ntot;
          N2 /= Ntot;
          N3 /= Ntot;
          vxI = (vx1[k][j-1][i]*N0+vx1[k][j-1][i+1]*N1+vx1[k][j][i]*N2+vx1[k][j][i+1]*N3);
          vyI = (vx2[k][j-1][i]*N0+vx2[k][j-1][i+1]*N1+vx2[k][j][i]*N2+vx2[k][j][i+1]*N3);
          vrI[p] = vxI*x1[i]/r + vyI*x2[j]/r;
        }
        #endif

        #if DIMENSIONS == 3
        phi   = atan2(x2[j],x1[i]);
        theta = acos(x3[k]/r);
        Ix1   = (r+ddr)*sin(theta)*cos(phi);
        Ix2   = (r+ddr)*sin(theta)*sin(phi);
        Ix3   = (r+ddr)*cos(theta);

        if (Ix1>x1[i] && Ix2>x2[j] && Ix3>x3[k]){
          xa  = x1[i];              ya  = x2[j];              za  = x3[k];
          xb  = x1[i+1];            yb  = x2[j+1];            zb  = x3[k+1];
          V0x = vx1[k][j][i];       V0y = vx2[k][j][i];       V0z = vx3[k][j][i];
          V1x = vx1[k][j][i+1];     V1y = vx2[k][j][i+1];     V1z = vx3[k][j][i+1];
          V2x = vx1[k][j+1][i];     V2y = vx2[k][j+1][i];     V2z = vx3[k][j+1][i];
          V3x = vx1[k][j+1][i+1];   V3y = vx2[k][j+1][i+1];   V3z = vx3[k][j+1][i+1];
          V4x = vx1[k+1][j][i];     V4y = vx2[k+1][j][i];     V4z = vx3[k+1][j][i];
          V5x = vx1[k+1][j][i+1];   V5y = vx2[k+1][j][i+1];   V5z = vx3[k+1][j][i+1];
          V6x = vx1[k+1][j+1][i];   V6y = vx2[k+1][j+1][i];   V6z = vx3[k+1][j+1][i];
          V7x = vx1[k+1][j+1][i+1]; V7y = vx2[k+1][j+1][i+1]; V7z = vx3[k+1][j+1][i+1];
        }
        if (Ix1<x1[i] && Ix2>x2[j] && Ix3>x3[k]){
          xa  = x1[i-1];            ya  = x2[j];              za  = x3[k];
          xb  = x1[i];              yb  = x2[j+1];            zb  = x3[k+1];
          V0x = vx1[k][j][i-1];     V0y = vx2[k][j][i-1];     V0z = vx3[k][j][i-1];
          V1x = vx1[k][j][i];       V1y = vx2[k][j][i];       V1z = vx3[k][j][i];
          V2x = vx1[k][j+1][i-1];   V2y = vx2[k][j+1][i-1];   V2z = vx3[k][j+1][i-1];
          V3x = vx1[k][j+1][i];     V3y = vx2[k][j+1][i];     V3z = vx3[k][j+1][i];
          V4x = vx1[k+1][j][i-1];   V4y = vx2[k+1][j][i-1];   V4z = vx3[k+1][j][i-1];
          V5x = vx1[k+1][j][i];     V5y = vx2[k+1][j][i];     V5z = vx3[k+1][j][i];
          V6x = vx1[k+1][j+1][i-1]; V6y = vx2[k+1][j+1][i-1]; V6z = vx3[k+1][j+1][i-1];
          V7x = vx1[k+1][j+1][i];   V7y = vx2[k+1][j+1][i];   V7z = vx3[k+1][j+1][i];
        }
        if (Ix1>x1[i] && Ix2<x2[j] && Ix3>x3[k]){
          xa  = x1[i];              ya  = x2[j-1];            za  = x3[k];
          xb  = x1[i+1];            yb  = x2[j];              zb  = x3[k+1];
          V0x = vx1[k][j-1][i];     V0y = vx2[k][j-1][i];     V0z = vx3[k][j-1][i];
          V1x = vx1[k][j-1][i+1];   V1y = vx2[k][j-1][i+1];   V1z = vx3[k][j-1][i+1];
          V2x = vx1[k][j][i];       V2y = vx2[k][j][i];       V2z = vx3[k][j][i];
          V3x = vx1[k][j][i+1];     V3y = vx2[k][j][i+1];     V3z = vx3[k][j][i+1];
          V4x = vx1[k+1][j-1][i];   V4y = vx2[k+1][j-1][i];   V4z = vx3[k+1][j-1][i];
          V5x = vx1[k+1][j-1][i+1]; V5y = vx2[k+1][j-1][i+1]; V5z = vx3[k+1][j-1][i+1];
          V6x = vx1[k+1][j][i];     V6y = vx2[k+1][j][i];     V6z = vx3[k+1][j][i];
          V7x = vx1[k+1][j][i+1];   V7y = vx2[k+1][j][i+1];   V7z = vx3[k+1][j][i+1];
        }
        if (Ix1<x1[i] && Ix2<x2[j] && Ix3>x3[k]){
          xa  = x1[i-1];            ya  = x2[j-1];            za  = x3[k];
          xb  = x1[i];              yb  = x2[j];              zb  = x3[k+1];
          V0x = vx1[k][j-1][i-1];   V0y = vx2[k][j-1][i-1];   V0z = vx3[k][j-1][i-1];
          V1x = vx1[k][j-1][i];     V1y = vx2[k][j-1][i];     V1z = vx3[k][j-1][i];
          V2x = vx1[k][j][i-1];     V2y = vx2[k][j][i-1];     V2z = vx3[k][j][i-1];
          V3x = vx1[k][j][i];       V3y = vx2[k][j][i];       V3z = vx3[k][j][i];
          V4x = vx1[k+1][j-1][i-1]; V4y = vx2[k+1][j-1][i-1]; V4z = vx3[k+1][j-1][i-1];
          V5x = vx1[k+1][j-1][i];   V5y = vx2[k+1][j-1][i];   V5z = vx3[k+1][j-1][i];
          V6x = vx1[k+1][j][i-1];   V6y = vx2[k+1][j][i-1];   V6z = vx3[k+1][j][i-1];
          V7x = vx1[k+1][j][i];     V7y = vx2[k+1][j][i];     V7z = vx3[k+1][j][i];
        }

        // zI < zP 
        if (Ix1>x1[i] && Ix2>x2[j] && Ix3<x3[k]){
          xa  = x1[i];              ya  = x2[j];              za  = x3[k-1];
          xb  = x1[i+1];            yb  = x2[j+1];            zb  = x3[k];
          V0x = vx1[k-1][j][i];     V0y = vx2[k-1][j][i];     V0z = vx3[k-1][j][i];
          V1x = vx1[k-1][j][i+1];   V1y = vx2[k-1][j][i+1];   V1z = vx3[k-1][j][i+1];
          V2x = vx1[k-1][j+1][i];   V2y = vx2[k-1][j+1][i];   V2z = vx3[k-1][j+1][i];
          V3x = vx1[k-1][j+1][i+1]; V3y = vx2[k-1][j+1][i+1]; V3z = vx3[k-1][j+1][i+1];
          V4x = vx1[k][j][i];       V4y = vx2[k][j][i];       V4z = vx3[k][j][i];
          V5x = vx1[k][j][i+1];     V5y = vx2[k][j][i+1];     V5z = vx3[k][j][i+1];
          V6x = vx1[k][j+1][i];     V6y = vx2[k][j+1][i];     V6z = vx3[k][j+1][i];
          V7x = vx1[k][j+1][i+1];   V7y = vx2[k][j+1][i+1];   V7z = vx3[k][j+1][i+1];
        }
        if (Ix1<x1[i] && Ix2>x2[j] && Ix3<x3[k]){
          xa  = x1[i-1];            ya  = x2[j];              za  = x3[k-1];
          xb  = x1[i];              yb  = x2[j+1];            zb  = x3[k];
          V0x = vx1[k-1][j][i-1];   V0y = vx2[k-1][j][i-1];   V0z = vx3[k-1][j][i-1];
          V1x = vx1[k-1][j][i];     V1y = vx2[k-1][j][i];     V1z = vx3[k-1][j][i];
          V2x = vx1[k-1][j+1][i-1]; V2y = vx2[k-1][j+1][i-1]; V2z = vx3[k-1][j+1][i-1];
          V3x = vx1[k-1][j+1][i];   V3y = vx2[k-1][j+1][i];   V3z = vx3[k-1][j+1][i];
          V4x = vx1[k][j][i-1];     V4y = vx2[k][j][i-1];     V4z = vx3[k][j][i-1];
          V5x = vx1[k][j][i];       V5y = vx2[k][j][i];       V5z = vx3[k][j][i];
          V6x = vx1[k][j+1][i-1];   V6y = vx2[k][j+1][i-1];   V6z = vx3[k][j+1][i-1];
          V7x = vx1[k][j+1][i];     V7y = vx2[k][j+1][i];     V7z = vx3[k][j+1][i];
        }
        if (Ix1>x1[i] && Ix2<x2[j] && Ix3<x3[k]){
          xa  = x1[i];              ya  = x2[j-1];            za  = x3[k-1];
          xb  = x1[i+1];            yb  = x2[j];              zb  = x3[k];
          V0x = vx1[k-1][j-1][i];   V0y = vx2[k-1][j-1][i];   V0z = vx3[k-1][j-1][i];
          V1x = vx1[k-1][j-1][i+1]; V1y = vx2[k-1][j-1][i+1]; V1z = vx3[k-1][j-1][i+1];
          V2x = vx1[k-1][j][i];     V2y = vx2[k-1][j][i];     V2z = vx3[k-1][j][i];
          V3x = vx1[k-1][j][i+1];   V3y = vx2[k-1][j][i+1];   V3z = vx3[k-1][j][i+1];
          V4x = vx1[k][j-1][i];     V4y = vx2[k][j-1][i];     V4z = vx3[k][j-1][i];
          V5x = vx1[k][j-1][i+1];   V5y = vx2[k][j-1][i+1];   V5z = vx3[k][j-1][i+1];
          V6x = vx1[k][j][i];       V6y = vx2[k][j][i];       V6z = vx3[k][j][i];
          V7x = vx1[k][j][i+1];     V7y = vx2[k][j][i+1];     V7z = vx3[k][j][i+1];
        }
        if (Ix1<x1[i] && Ix2<x2[j] && Ix3<x3[k]){
          xa  = x1[i-1];            ya  = x2[j-1];            za  = x3[k-1];
          xb  = x1[i];              yb  = x2[j];              zb  = x3[k];
          V0x = vx1[k-1][j-1][i-1]; V0y = vx2[k-1][j-1][i-1]; V0z = vx3[k-1][j-1][i-1];
          V1x = vx1[k-1][j-1][i];   V1y = vx2[k-1][j-1][i];   V1z = vx3[k-1][j-1][i];
          V2x = vx1[k-1][j][i-1];   V2y = vx2[k-1][j][i-1];   V2z = vx3[k-1][j][i-1];
          V3x = vx1[k-1][j][i];     V3y = vx2[k-1][j][i];     V3z = vx3[k-1][j][i];
          V4x = vx1[k][j-1][i-1];   V4y = vx2[k][j-1][i-1];   V4z = vx3[k][j-1][i-1];
          V5x = vx1[k][j-1][i];     V5y = vx2[k][j-1][i];     V5z = vx3[k][j-1][i];
          V6x = vx1[k][j][i-1];     V6y = vx2[k][j][i-1];     V6z = vx3[k][j][i-1];
          V7x = vx1[k][j][i];       V7y = vx2[k][j][i];       V7z = vx3[k][j][i];
        }

        N0     = ((xb-Ix1)*(yb-Ix2)*(zb-Ix3))/Ntot;
        N1     = ((Ix1-xa)*(yb-Ix2)*(zb-Ix3))/Ntot;
        N2     = ((xb-Ix1)*(Ix2-ya)*(zb-Ix3))/Ntot;
        N3     = ((Ix1-xa)*(Ix2-ya)*(zb-Ix3))/Ntot;
        N4     = ((xb-Ix1)*(yb-Ix2)*(Ix3-za))/Ntot;
        N5     = ((Ix1-xa)*(yb-Ix2)*(Ix3-za))/Ntot;
        N6     = ((xb-Ix1)*(Ix2-ya)*(Ix3-za))/Ntot;
        N7     = ((Ix1-xa)*(Ix2-ya)*(Ix3-za))/Ntot;				
        vxI    = (V0x*N0+V1x*N1+V2x*N2+V3x*N3+V4x*N4+V5x*N5+V6x*N6+V7x*N7);
        vyI    = (V0y*N0+V1y*N1+V2y*N2+V3y*N3+V4y*N4+V5y*N5+V6y*N6+V7y*N7);
        vzI    = (V0z*N0+V1z*N1+V2z*N2+V3z*N3+V4z*N4+V5z*N5+V6z*N6+V7z*N7);
        vrI[p] = (vxI*x1[i]+vyI*x2[j]+vzI*x3[k])/r;			
        #endif

      }

      vr    = EXPAND(vx1[k][j][i]*x1[i]/r, + vx2[k][j][i]*x2[j]/r, + vx3[k][j][i]*x3[k]/r);
      P00   = vrI[0];
      P11   = vr;
      P22   = vrI[1];
      P01   = ((0.5*fabs(ddr)*P00)+(0.5*fabs(ddr)*P11))/fabs(ddr);
      P12   = ((0.5*fabs(ddr)*P11)+(0.5*fabs(ddr)*P22))/fabs(ddr);
      dvdr  = fabs(((1.0/12.0)*P00)-((2.0/3.0)*P01)+((2.0/3.0)*P12) -((1.0/12.0)*P22))/fabs(0.5*ddr);
      nu2_c = (1.0-(1.0/(r*r)));
      B     = ((d->Vc[RHO][k][j][i])*Q*c*ke);
      sigma = (r/fabs(vr))*(dvdr)-1.0;
      f     = ((pow(1.0+sigma,1.0+a)-pow(1.0+sigma*nu2_c,1.0+a))/((1.0+a)*(1.0-nu2_c)*sigma*pow(1.0+sigma,a)));

      #if AMR == YES
      #if DIMENSIONS == 3
      if (i < 4 || j < 4 || k < 4) {
        gL = 0.0;
      } else if (i > IEND-3 || j > JEND-3 || k > KEND-3){
        gL = 0.0;
      } else {
        gL = (f*A*pow(r,-2)*pow(dvdr/B,a));
      }
      #endif
      #if DIMENSIONS == 2
      if (i < 4 || j < 4) {
        gL = 0.0;
      } else if (i > IEND-3 || j > JEND-3){
        gL = 0.0;
      } else {
        gL = (f*A*pow(r,-2)*pow(dvdr/B,a));
        //printf("gl = %f, f = %f, A = %f, r = %f, dvdr = %f, B = %f, a = %f \n", gL, f, A, r, dvdr, B, a);
      }
      #endif
      #endif

      #if AMR == NO
      gL = (f*A*pow(r,-2)*pow(dvdr/B,a));
      #endif
      D_EXPAND(vx1[k][j][i] += (gL*x1[i]/r)*g_dt;,
               vx2[k][j][i] += (gL*x2[j]/r)*g_dt;,
               vx3[k][j][i] += (gL*x3[k]/r)*g_dt;)

/*
      #if DIMENSIONS == 2
      if (d->Vc[PRS][k][j][i] < 0.0){
        d->Vc[PRS][k][j][i] = (d->Vc[PRS][k][j][i+1] + 
                               d->Vc[PRS][k][j][i-1] +
                               d->Vc[PRS][k][j+1][i] +
                               d->Vc[PRS][k][j-1][i])/4.0;
      }
      if (d->Vc[RHO][k][j][i] < 0.0){
        d->Vc[RHO][k][j][i] = (d->Vc[RHO][k][j][i+1] + 
                               d->Vc[RHO][k][j][i-1] +
                               d->Vc[RHO][k][j+1][i] +
                               d->Vc[RHO][k][j-1][i])/4.0;
      }
      #endif

      #if DIMENSIONS == 3
        if (d->Vc[PRS][k][j][i] < 0.0){
          d->Vc[PRS][k][j][i] = (d->Vc[PRS][k][j][i+1] + 
                                 d->Vc[PRS][k][j][i-1] +
                                 d->Vc[PRS][k][j+1][i] +
                                 d->Vc[PRS][k][j-1][i] +
                                 d->Vc[PRS][k+1][j][i] + 
                                 d->Vc[PRS][k-1][j][i]  )/6.0;
      }
      if (d->Vc[RHO][k][j][i] < 0.0){
        d->Vc[RHO][k][j][i] = (d->Vc[RHO][k][j][i+1] + 
                               d->Vc[RHO][k][j][i-1] +
                               d->Vc[RHO][k][j+1][i] +
                               d->Vc[RHO][k][j-1][i] +
                               d->Vc[RHO][k+1][j][i] + 
                               d->Vc[RHO][k-1][j][i]  )/6.0;
      }
      #endif
*/

    } // end of if r > shell.


    #if EOS == IDEAL
    if (d->Vc[PRS][k][j][i] < d->Vc[RHO][k][j][i]*T/(KELVIN*mu)) {
      d->Vc[PRS][k][j][i] = (d->Vc[RHO][k][j][i]*T/(KELVIN*mu));
    }
    #endif

  }} // end of DOM loop.
} // End of function.



#if BODY_FORCE != NO
#if CAK == YES
void BodyForceVector(double *rubrik, double *v, double *g, double *x_box)
{
  double x1, x2, x3;
  /*
  x1 = *x_box[0][1][1][1];
  x2 = *x_box[1][1][1][1];
  x3 = *x_box[2][1][1][1];
  */ 
#endif
#if CAK == NO
void BodyForceVector(double *v, double *g, double x1, double x2, double x3)
{
#endif
  
  double M_star, Edd, Omega, Mratio, Lratio;
  double r2, r, omega_fr, Omega_star;
  double Fin_x1, Fin_x2, gg, g_in;

  Mratio  = g_inputParam[M_RATIO];
  Lratio  = g_inputParam[L_RATIO];
  Omega   = g_inputParam[OMEGA];

  M_star  = (Mratio*CONST_Msun/UNIT_MASS);
  Edd     = (2.6e-5*(Lratio)*(1.0/Mratio));
  Omega_star  = Omega*sqrt((8.0/27.0)*UNIT_G*M_star);

  /* - Rotational frequency (orbit and frame)) - */
  omega_fr = Omega_star;

  /* - Distance from star - */
  r2 = EXPAND(x1*x1, + x2*x2, + x3*x3);
  r = sqrt(r2);

  /* - Gravity outside bodies - */
  gg = -UNIT_G*(M_star - Edd)/r/r;

  /* - Gravity inside bodies - */
  g_in = -(4.0/3.0)*CONST_PI*UNIT_G*v[RHO];     

  /* - Coriolis and centrifugal forces - */
  #if DIMENSIONS == 2
  Fin_x1 = omega_fr*omega_fr*x1;
  Fin_x2 = 0.0;
  #endif
  #if DIMENSIONS == 3
  Fin_x1 = omega_fr*omega_fr*x1 + 2.0*omega_fr*v[VX2];
  Fin_x2 = omega_fr*omega_fr*x2 - 2.0*omega_fr*v[VX1];
  #endif

  if (r >= 1.0){ /* - External gravity + centrifugal + coriolis - */
    g[IDIR] = gg*x1/r + Fin_x1;
    g[JDIR] = gg*x2/r + Fin_x2;
    g[KDIR] = gg*x3/r;
  } else if (r < 1.0) { /* - Star interal gravity - */
    g[IDIR] = g_in*x1 + Fin_x1;
    g[JDIR] = g_in*x2 + Fin_x2;
    g[KDIR] = g_in*x3;
  }
}
#endif
/*================================================================================*/
