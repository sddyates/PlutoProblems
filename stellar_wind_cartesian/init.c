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

double TwoDimensionalInterp(const Data *d, RBox *box, Grid *grid, 
                            int i, int j, int k, double *x1, double *x2, double r, double ddr);
double ThreeDimensionalInterp(const Data *d, RBox *box, Grid *grid, 
                            int i, int j, int k, double *x1, double *x2, double *x3,
                            double r, double ddr);
double VelocityStellarSurface2D(const Data *d, RBox *box, Grid *grid, 
                            int i, int j, int k, double *x1, double *x2, double *x3, double r);
double VelocityStellarSurface3D(const Data *d, RBox *box, Grid *grid, 
                            int i, int j, int k, double *x1, double *x2, double *x3, double r);

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
  v_inf     = 3.0*v_esc,//v_esc * sqrt((a/(1.0-a))),                                      
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

  double Cs_p, Mratio, Lratio, T, mu, a, b, Q, a_eff, beta, Omega, shell, M_star;
  double Edd, L, c, M_dot, ke, Omega2, A, Bcgs, cs, Bq, v_esc, v_inf;

  double       bx = 0.0,                by = 0.0,                bz = 0.0;
  double      bxp = 0.0,               byp = 0.0,               bzp = 0.0;
  double        x = 0.0,                 y = 0.0,                 z = 0.0;
  double       xp = 0.0,                yp = 0.0,                zp = 0.0;
  double      *x1 = grid[IDIR].x,      *x2 = grid[JDIR].x,      *x3 = grid[KDIR].x;
  double     *dx1 = grid[IDIR].dx,    *dx2 = grid[JDIR].dx,    *dx3 = grid[KDIR].dx;
  double   ***vx1 = d->Vc[VX1],     ***vx2 = d->Vc[VX2],     ***vx3 = d->Vc[VX3];
  double       dx = 0.0,                dy = 0.0,                dz = 0.0;

  double P00 = 0.0, P01 = 0.0, P11 = 0.0, P12 = 0.0, P22 = 0.0; 
  double dvdr = 0.0, dvdr2 = 0.0, Beq = 0.0, vel_mag;
  double vr = 0.0, vr2 = 0.0, r = 0.0, r2 = 0.0, theta = 0.0, phi = 0.0;
  double dr2 = 0.0, dr = 0.0, ddr = 0.0, rb = 0.0, rb2 = 0.0, rp = 0.0, rp2 = 0.0, dI2 = 0.0, dI = 0.0, vrI[2];
  double gL = 0.0, nu2_c = 0.0, B = 0.0, sigma = 0.0, f = 0.0, vv = 0.0;

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
  v_inf   = 3.0*v_esc;//v_esc * sqrt((a/(1.0-a)));

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

      dx = x1[i+1] - x1[i];
      dy = x2[j+1] - x2[j];
      dz = x3[k+1] - x3[k];
      dr2  = EXPAND(dx*dx,+dy*dy,+dz*dz);
      dr   = sqrt(dr2);

      #if DIMENSIONS == 2
      if (r < 1.0 && r + 5.0*ddr > 1.0) {
        vel_mag = TwoDimensionalInterp(d, box, grid, i, j, k, x1, x2, r, dr);
      } else {
        vel_mag = 0.0;
      }
      #endif

      #if DIMENSIONS == 3
      if (r < 1.0 && r + 5.0*ddr > 1.0) {
        vel_mag = ThreeDimensionalInterp(d, box, grid, i, j, k, x1, x2, x3, r, dr);
      } else {
        vel_mag = 0.0;
      }
      #endif

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
      //dr2  = EXPAND(dx1[i]*dx1[i],+dx2[j]*dx2[j],+dx3[k]*dx3[k]);
      //dr   = sqrt(dr2);

      dx = x1[i+1] - x1[i];
      dy = x2[j+1] - x2[j];
      dz = x3[k+1] - x3[k];
      dr2  = EXPAND(dx*dx,+dy*dy,+dz*dz);
      dr   = sqrt(dr2);

      for (p = 0; p < 2; p++) {

        if (p == 0) {
          ddr = -0.9*dr;
        } else {
          ddr = 0.9*dr;
        }

        /* - Find nearest nighbours and record coordinates - */
        #if DIMENSIONS == 2
        vrI[p] = TwoDimensionalInterp(d, box, grid, i, j, k, x1, x2, r, ddr);
        #endif

        #if DIMENSIONS == 3
        vrI[p] = ThreeDimensionalInterp(d, box, grid, i, j, k, x1, x2, x3, r, ddr);
        #endif
      }

      vr    = EXPAND(vx1[k][j][i]*x1[i]/r, + vx2[k][j][i]*x2[j]/r, + vx3[k][j][i]*x3[k]/r);
      P00   = vrI[0];
      P11   = vr;
      P22   = vrI[1];
      P01   = ((0.5*fabs(ddr)*P00)+(0.5*fabs(ddr)*P11))/fabs(ddr);
      P12   = ((0.5*fabs(ddr)*P11)+(0.5*fabs(ddr)*P22))/fabs(ddr);
      dvdr  = fabs(((1.0/12.0)*P00)-((2.0/3.0)*P01)+((2.0/3.0)*P12) -((1.0/12.0)*P22))/fabs(0.5*ddr);
      //dvdr = fabs((vrI[1] - vrI[0])/(2.0*ddr));
      nu2_c = (1.0-(1.0/(r*r)));
      B     = ((d->Vc[RHO][k][j][i])*Q*c*ke);
      sigma = (r/fabs(vr))*(dvdr)-1.0;
      f     = ((pow(1.0+sigma,1.0+a)-pow(1.0+sigma*nu2_c,1.0+a))/((1.0+a)*(1.0-nu2_c)*sigma*pow(1.0+sigma,a)));
      gL = fabs(f*A*pow(r,-2)*pow(dvdr/B,a));

      #if AMR_ON == YES
      #if DIMENSIONS == 3
      if (i < 4 || j < 4 || k < 4) {
        gL = 0.0;
      } else if (i > IEND-3 || j > JEND-3 || k > KEND-3){
        gL = 0.0;
      }
      #endif
      #if DIMENSIONS == 2
      if (i < 4 || j < 4) {
        gL = 0.0;
      } else if (i > IEND-3 || j > JEND-3){
        gL = 0.0;
        //printf("gl = %f, f = %f, A = %f, r = %f, dvdr = %f, B = %f, a = %f \n", gL, f, A, r, dvdr, B, a);
      }
      #endif
      #endif

      D_EXPAND(vx1[k][j][i] += (gL*x1[i]/r)*g_dt;,
               vx2[k][j][i] += (gL*x2[j]/r)*g_dt;,
               vx3[k][j][i] += (gL*x3[k]/r)*g_dt;)

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

double TwoDimensionalInterp(const Data *d, RBox *box, Grid *grid, 
                            int i, int j, int k, double *x1, double *x2, double r, double ddr)
{
  /*
   *       _____________________________
   *  j+1 |              |              | 
   *      |              |              |
   *      |              |     *        |
   *      |              |              |
   *      |              |              |
   *    j |______________|______________|
   *      |              |              | 
   *      |              |              |
   *      |              |              |
   *      |              |              |
   *      |              |              |
   *  j-1 |______________|______________|
   *  
   *    i - 1          i            i + 1
   *  
   *  
   *  yb 3               4
   *      |```|`````````|
   *      |   |         |
   *   yI |___|_________|
   *      |   |         |
   *      |   |         |
   *      |___|_________|
   *  ya 1    xI         2
   *      xa            xb
   *
   * The interpolation points are always between xa, xb and ya, yb.
   *       
   * ddr is the distance forwards (backwards) from the point 
   * that the velocity gradient is needed to give the place 
   * to perform the interpolation.
   *
   *  r -> r±ddr
   *
   * |'''''''''''''''''''''''''|
   * |                         |
   * |                         |
   * |      * r+ddr            |
   * |     /|                  |
   * |    / |                  |
   * |   /  | r+ddr*cos(theta) |
   * |  /   |                  |
   * | /    |                  |
   * |/_____|__________________|
   * r     r+ddr*sin(theta)
   *
   * xI and yI are the interpolation componets. 
   * 
   */

  int u, s;
  double Ntot; // total volume of interpolation "space".
  double vrI; // final interpolated radial-velocity. 
  double vI[2]; // interpolated velocity components.
  double N[4]; // Areas used to waight nearest neighbours (NN).
  double V[2][4]; // Array to hold velocity componants at NN.
  double xa, xb; // bracketing x values.
  double ya, yb; // bracketing y values.
  double theta = atan2(x1[i],x2[j]); // convert to polar from Cartesian.
  double xI = (r+ddr)*sin(theta);
  double yI = (r+ddr)*cos(theta);

  /*
   * Eath of the following if statments checks which quadrent 
   * the interpolation point is in and gets the components 
   * of the velocity and the bracketing x and y values.
   */
  if (xI > x1[i] && yI > x2[j]){
    xa = x1[i]; xb = x1[i+1];
    ya = x2[j]; yb = x2[j+1];
    for (u = VX1; u < VX2+1; u++) {
      V[u-1][0] = d->Vc[u][k][j][i];
      V[u-1][1] = d->Vc[u][k][j][i+1];
      V[u-1][2] = d->Vc[u][k][j+1][i];
      V[u-1][3] = d->Vc[u][k][j+1][i+1];      
    }
  } else if (xI < x1[i] && yI > x2[j]){
    xa = x1[i-1]; xb = x1[i];
    ya = x2[j]; yb = x2[j+1];
    for (u = VX1; u < VX2+1; u++) {
      V[u-1][0] = d->Vc[u][k][j][i-1];
      V[u-1][1] = d->Vc[u][k][j][i];
      V[u-1][2] = d->Vc[u][k][j+1][i-1];
      V[u-1][3] = d->Vc[u][k][j+1][i];
    }
  } else if (xI < x1[i] && yI < x2[j]){
    xa = x1[i-1]; xb = x1[i];
    ya = x2[j-1]; yb = x2[j];
    for (u = VX1; u < VX2+1; u++) {
      V[u-1][0] = d->Vc[u][k][j-1][i-1];
      V[u-1][1] = d->Vc[u][k][j-1][i];
      V[u-1][2] = d->Vc[u][k][j][i-1];
      V[u-1][3] = d->Vc[u][k][j][i];
    }
  } else if (xI > x1[i] && yI < x2[j]){
    xa = x1[i]; xb = x1[i+1];
    ya = x2[j-1]; yb = x2[j];
    for (u = VX1; u < VX2+1; u++) {
      V[u-1][0] = d->Vc[u][k][j-1][i];
      V[u-1][1] = d->Vc[u][k][j-1][i+1];
      V[u-1][2] = d->Vc[u][k][j][i];
      V[u-1][3] = d->Vc[u][k][j][i+1];
    }
  }

  // Find total volume.
  N[0] = (xb - xI)*(yb - yI);
  N[1] = (xI - xa)*(yb - yI);
  N[2] = (xb - xI)*(yI - ya);
  N[3] = (xI - xa)*(yI - ya);
  Ntot = N[0] + N[1] + N[2] + N[3];
  // Normalise volumes by total.
  N[0] /= Ntot; 
  N[1] /= Ntot; 
  N[2] /= Ntot; 
  N[3] /= Ntot;

  // ==========================================
  // vI contains the interpolated velovities 
  // 
  // vI[0] = x velocity 
  // vI[1] = y "
  //
  // V contains the velocity componants at 
  // each of the sample points.
  // 
  // V[0, :] = x velocity at each sample point.
  // V[1, :] = y velocity at each sample point.
  //
  // N contains the waighting volumes for the 
  // the corner velocities.
  // ==========================================

  vI[0] = 0.0; vI[1] = 0.0;
  for (s = 0; s < 4; s++) {
    vI[0] += V[0][s]*N[s];
    vI[1] += V[1][s]*N[s];
  }
  vrI = (vI[0]*x1[i] + vI[1]*x2[j])/r;
  return vrI;

}


double ThreeDimensionalInterp(const Data *d, RBox *box, Grid *grid, 
                            int i, int j, int k, double *x1, double *x2, double *x3, double r, double ddr)
{
  int u, s;
  double Ntot;
  double vrI;
  double vI[3];
  double N[8];
  double V[3][8];
  double xa, xb;
  double ya, yb;
  double za, zb;
  double phi = atan2(x2[j],x1[i]);
  double theta = acos(x3[k]/r);
  double xI = (r+ddr)*sin(theta)*cos(phi);
  double yI = (r+ddr)*sin(theta)*sin(phi);
  double zI = (r+ddr)*cos(theta);
  int tag_if;

  if (xI > x1[i] && yI > x2[j] && zI > x3[k]){
    tag_if = 1;
    xa = x1[i]; xb = x1[i+1];
    ya = x2[j]; yb = x2[j+1];
    za = x3[k]; zb = x3[k+1];
    for (u = VX1; u < VX3+1; u++) {
      V[u-1][0] = d->Vc[u][k][j][i];
      V[u-1][1] = d->Vc[u][k][j][i+1];
      V[u-1][2] = d->Vc[u][k][j+1][i];
      V[u-1][3] = d->Vc[u][k][j+1][i+1];
      V[u-1][4] = d->Vc[u][k+1][j][i];
      V[u-1][5] = d->Vc[u][k+1][j][i+1];
      V[u-1][6] = d->Vc[u][k+1][j+1][i];
      V[u-1][7] = d->Vc[u][k+1][j+1][i+1];
    }
  }
  if (xI < x1[i] && yI > x2[j] && zI > x3[k]){
    tag_if = 2;
    xa = x1[i-1]; xb = x1[i];
    ya = x2[j]; yb = x2[j+1];
    za = x3[k]; zb = x3[k+1];
    for (u = VX1; u < VX3+1; u++) {
      V[u-1][0] = d->Vc[u][k][j][i-1];
      V[u-1][1] = d->Vc[u][k][j][i];
      V[u-1][2] = d->Vc[u][k][j+1][i-1];
      V[u-1][3] = d->Vc[u][k][j+1][i];
      V[u-1][4] = d->Vc[u][k+1][j][i-1];
      V[u-1][5] = d->Vc[u][k+1][j][i];
      V[u-1][6] = d->Vc[u][k+1][j+1][i-1];
      V[u-1][7] = d->Vc[u][k+1][j+1][i];
    }
  }
  if (xI > x1[i] && yI < x2[j] && zI > x3[k]){
    tag_if = 3;
    xa = x1[i]; xb = x1[i+1];
    ya = x2[j-1]; yb = x2[j];
    za = x3[k]; zb = x3[k+1];
    for (u = VX1; u < VX3+1; u++) {
      V[u-1][0] = d->Vc[u][k][j-1][i];
      V[u-1][1] = d->Vc[u][k][j-1][i+1];
      V[u-1][2] = d->Vc[u][k][j][i];
      V[u-1][3] = d->Vc[u][k][j][i+1];
      V[u-1][4] = d->Vc[u][k+1][j-1][i];
      V[u-1][5] = d->Vc[u][k+1][j-1][i+1];
      V[u-1][6] = d->Vc[u][k+1][j][i];
      V[u-1][7] = d->Vc[u][k+1][j][i+1];
    }
  }
  if (xI < x1[i] && yI < x2[j] && zI > x3[k]){
    tag_if = 4;
    xa = x1[i-1]; xb = x1[i];
    ya = x2[j-1]; yb = x2[j];
    za = x3[k]; zb = x3[k+1];
    for (u = VX1; u < VX3+1; u++) {
      V[u-1][0] = d->Vc[u][k][j-1][i-1];
      V[u-1][1] = d->Vc[u][k][j-1][i];
      V[u-1][2] = d->Vc[u][k][j][i-1];
      V[u-1][3] = d->Vc[u][k][j][i];
      V[u-1][4] = d->Vc[u][k+1][j-1][i-1];
      V[u-1][5] = d->Vc[u][k+1][j-1][i];
      V[u-1][6] = d->Vc[u][k+1][j][i-1];
      V[u-1][7] = d->Vc[u][k+1][j][i];
    }
  }

  // zI < zP 
  if (xI > x1[i] && yI > x2[j] && zI < x3[k]){
    tag_if = 5;
    xa = x1[i]; xb = x1[i+1];
    ya = x2[j]; yb = x2[j+1];
    za = x3[k-1]; zb = x3[k];
    for (u = VX1; u < VX3+1; u++) {
      V[u-1][0] = d->Vc[u][k-1][j][i];
      V[u-1][1] = d->Vc[u][k-1][j][i+1];
      V[u-1][2] = d->Vc[u][k-1][j+1][i];
      V[u-1][3] = d->Vc[u][k-1][j+1][i+1];
      V[u-1][4] = d->Vc[u][k][j][i];
      V[u-1][5] = d->Vc[u][k][j][i+1];
      V[u-1][6] = d->Vc[u][k][j+1][i];
      V[u-1][7] = d->Vc[u][k][j+1][i+1];
    }
  }
  if (xI < x1[i] && yI > x2[j] && zI < x3[k]){
    tag_if = 6;
    xa = x1[i-1]; xb = x1[i];
    ya = x2[j]; yb = x2[j+1];
    za = x3[k-1]; zb = x3[k];
    for (u = VX1; u < VX3+1; u++) {
      V[u-1][0] = d->Vc[u][k-1][j][i-1];
      V[u-1][1] = d->Vc[u][k-1][j][i];
      V[u-1][2] = d->Vc[u][k-1][j+1][i-1];
      V[u-1][3] = d->Vc[u][k-1][j+1][i];
      V[u-1][4] = d->Vc[u][k][j][i-1];
      V[u-1][5] = d->Vc[u][k][j][i];
      V[u-1][6] = d->Vc[u][k][j+1][i-1];
      V[u-1][7] = d->Vc[u][k][j+1][i];
    }
  }
  if (xI > x2[i] && yI < x2[j] && zI < x3[k]){
    tag_if = 7;
    xa = x1[i]; xb = x1[i+1];
    ya = x2[j-1]; yb = x2[j];
    za = x3[k-1]; zb = x3[k];
    for (u = VX1; u < VX3+1; u++) {
      V[u-1][0] = d->Vc[u][k-1][j-1][i];
      V[u-1][1] = d->Vc[u][k-1][j-1][i+1];
      V[u-1][2] = d->Vc[u][k-1][j][i];
      V[u-1][3] = d->Vc[u][k-1][j][i+1];
      V[u-1][4] = d->Vc[u][k][j-1][i];
      V[u-1][5] = d->Vc[u][k][j-1][i+1];
      V[u-1][6] = d->Vc[u][k][j][i];
      V[u-1][7] = d->Vc[u][k][j][i+1];
    }
  }
  if (xI < x1[i] && yI < x2[j] && zI < x3[k]){
    tag_if = 8;
    xa = x1[i-1]; xb = x1[i];
    ya = x2[j-1]; yb = x2[j];
    za = x3[k-1]; zb = x3[k];
    for (u = VX1; u < VX3+1; u++) {
      V[u-1][0] = d->Vc[u][k-1][j-1][i-1];
      V[u-1][1] = d->Vc[u][k-1][j-1][i];
      V[u-1][2] = d->Vc[u][k-1][j][i-1];
      V[u-1][3] = d->Vc[u][k-1][j][i];
      V[u-1][4] = d->Vc[u][k][j-1][i-1];
      V[u-1][5] = d->Vc[u][k][j-1][i];
      V[u-1][6] = d->Vc[u][k][j][i-1];
      V[u-1][7] = d->Vc[u][k][j][i];
    }
  }

  // Find total volume.
  N[0] = (xb - xI)*(yb - yI)*(zb - zI);
  N[1] = (xI - xa)*(yb - yI)*(zb - zI);
  N[2] = (xb - xI)*(yI - ya)*(zb - zI);
  N[3] = (xI - xa)*(yI - ya)*(zb - zI);
  N[4] = (xb - xI)*(yb - yI)*(zI - za);
  N[5] = (xI - xa)*(yb - yI)*(zI - za);
  N[6] = (xb - xI)*(yI - ya)*(zI - za);
  N[7] = (xI - xa)*(yI - ya)*(zI - za);
  Ntot = N[0] + N[1] + N[2] + N[3] + N[4] + N[5] + N[6] + N[7];
  // Normalise volumes by total.
  N[0] /= Ntot; 
  N[1] /= Ntot; 
  N[2] /= Ntot; 
  N[3] /= Ntot;
  N[4] /= Ntot; 
  N[5] /= Ntot; 
  N[6] /= Ntot; 
  N[7] /= Ntot;
         
  // ==========================================
  // vI contains the interpolated velovities 
  // 
  // vI[0] = x velocity 
  // vI[1] = y "
  // vI[2] = z "
  //
  // V contains the velocity componants at 
  // each of the sample points.
  // 
  // V[0][:] = x velocity at each sample point.
  // V[1][:] = y velocity at each sample point.
  // V[2][:] = z velocity at each sample point.
  //
  // N contains the waighting volumes for the 
  // the corner velocities.
  // ==========================================

  vI[0] = 0.0; vI[1] = 0.0; vI[2] = 0.0;
  for (s = 0; s < 8; s++) {
    vI[0] += V[0][s]*N[s];
    vI[1] += V[1][s]*N[s];
    vI[2] += V[2][s]*N[s];
  }
  vrI = (vI[0]*x1[i] + vI[1]*x2[j] + vI[2]*x3[k])/r;

  if(isnan(vrI)){
    printf("tag_if=%i \n",tag_if);
    printf("vrI=%f, vI[0]=%f, vI[1]=%f, vI[2]=%f Ntot=%f \n",vrI,vI[0],vI[1],vI[2],Ntot);
    printf("N0=%f, N1=%f, N2=%f, N3=%f, N4=%f, N5=%f, N6=%f, N7=%f \n",N[0], N[1], N[2], N[3], N[4], N[5], N[6], N[7]);
    printf("xa=%f, xI=%f, xb=%f \n", xa, xI, xb);
    printf("ya=%f, yI=%f, yb=%f \n", ya, yI, yb);
    printf("za=%f, zI=%f, zb=%f \n", za, zI, zb);
    printf("\n");
  }
  return vrI;
}


double VelocityStellarSurface2D(const Data *d, RBox *box, Grid *grid, 
                            int i, int j, int k, double *x1, double *x2, double *x3, double r){

  int l;
  double vel_x1, vel_x2, vel_x3;
  double r_testx, r_testy, r_testz, vel_mag;
  double dx, dy, dz, dr2, dr, ddr;

  dx = x1[i+1] - x1[i];
  dy = x2[j+1] - x2[j];
  dr2  = EXPAND(dx*dx,+dy*dy,+dz*dz);
  dr   = sqrt(dr2);
  ddr = 0.9*dr;

  if (r < 1.0 && r + ddr > 1.0) {
    vel_mag = TwoDimensionalInterp(d, box, grid, i, j, k, x1, x2, r, ddr);
  }


  for (l=1; l<2; l++){

    // Upper right.
    if (x1[i] > 0.0 && x2[j] > 0.0) {
      r_testx = D_EXPAND(x1[i+l]*x1[i+l], + x2[j]*x2[j], + x3[k]*x3[k]);
      r_testx = sqrt(r_testx);
      r_testy = D_EXPAND(x1[i]*x1[i], + x2[j+l]*x2[j+l], + x3[k]*x3[k]);
      r_testy = sqrt(r_testy);
      if (r < 1.0 && (r_testx > 1.0 || r_testy > 1.0)) {
        vel_mag = TwoDimensionalInterp(d, box, grid, i, j, k, x1, x2, r, ddr);
      }
/*
      if (r < 1.0 && r_testx > 1.0 && r_testy > 1.0) {
        // these are the cells at the boundary.
        D_EXPAND(vel_x1 = d->Vc[VX1][k][j][i+l];,
                 vel_x2 = d->Vc[VX2][k][j+l][i];,
                 vel_x3 = d->Vc[VX3][k][j][i];)   
        vel_mag = sqrt(vel_x1*vel_x1 + vel_x2*vel_x2);
      } else if (r < 1.0 && r_testx < 1.0 && r_testy > 1.0) {
        D_EXPAND(vel_x1 = d->Vc[VX1][k][j][i];,
                 vel_x2 = d->Vc[VX2][k][j+l][i];,
                 vel_x3 = d->Vc[VX3][k][j][i];)          
        vel_mag = sqrt(vel_x1*vel_x1 + vel_x2*vel_x2);
      } else if (r < 1.0 && r_testx > 1.0 && r_testy < 1.0) {
        D_EXPAND(vel_x1 = d->Vc[VX1][k][j][i+l];,
                 vel_x2 = d->Vc[VX2][k][j][i];,
                 vel_x3 = d->Vc[VX3][k][j][i];)          
        vel_mag = sqrt(vel_x1*vel_x1 + vel_x2*vel_x2);
      }
*/
    }
    // Lower right.
    else if (x1[i] > 0.0 && x2[j] < 0.0) {
      r_testx = D_EXPAND(x1[i+l]*x1[i+l], + x2[j]*x2[j], + x3[k]*x3[k]);
      r_testx = sqrt(r_testx);
      r_testy = D_EXPAND(x1[i]*x1[i], + x2[j-l]*x2[j-l], + x3[k]*x3[k]);
      r_testy = sqrt(r_testy);
      if (r < 1.0 && (r_testx > 1.0 || r_testy > 1.0)) {
        vel_mag = TwoDimensionalInterp(d, box, grid, i, j, k, x1, x2, r, ddr);
      }
/*
      if (r < 1.0 && r_testx > 1.0 && r_testy > 1.0) {
        // these are the cells at the boundary.
        D_EXPAND(vel_x1 = d->Vc[VX1][k][j][i+l];,
                 vel_x2 = d->Vc[VX2][k][j-l][i];,
                 vel_x3 = d->Vc[VX3][k][j][i];)          
        vel_mag = sqrt(vel_x1*vel_x1 + vel_x2*vel_x2);
      } else if (r < 1.0 && r_testx < 1.0 && r_testy > 1.0) {
        D_EXPAND(vel_x1 = d->Vc[VX1][k][j][i];,
                 vel_x2 = d->Vc[VX2][k][j-l][i];,
                 vel_x3 = d->Vc[VX3][k][j][i];)          
        vel_mag = sqrt(vel_x1*vel_x1 + vel_x2*vel_x2);
      } else if (r < 1.0 && r_testx > 1.0 && r_testy < 1.0) {
        D_EXPAND(vel_x1 = d->Vc[VX1][k][j][i+l];,
                 vel_x2 = d->Vc[VX2][k][j][i];,
                 vel_x3 = d->Vc[VX3][k][j][i];)          
        vel_mag = sqrt(vel_x1*vel_x1 + vel_x2*vel_x2);
      }
*/
    }
    // lower left.
    else if (x1[i] < 0.0 && x2[j] < 0.0) {
      r_testx = D_EXPAND(x1[i-l]*x1[i-l], + x2[j]*x2[j], + x3[k]*x3[k]);
      r_testx = sqrt(r_testx);
      r_testy = D_EXPAND(x1[i]*x1[i], + x2[j-l]*x2[j-l], + x3[k]*x3[k]);
      r_testy = sqrt(r_testy);
      if (r < 1.0 && (r_testx > 1.0 || r_testy > 1.0)) {
        vel_mag = TwoDimensionalInterp(d, box, grid, i, j, k, x1, x2, r, ddr);
      }
/*
      if (r < 1.0 && r_testx > 1.0 && r_testy > 1.0) {
        // these are the cells at the boundary.
        D_EXPAND(vel_x1 = d->Vc[VX1][k][j][i-l];,
                 vel_x2 = d->Vc[VX2][k][j-l][i];,
                 vel_x3 = d->Vc[VX3][k][j][i];)          
        vel_mag = sqrt(vel_x1*vel_x1 + vel_x2*vel_x2);
      } else if (r < 1.0 && r_testx < 1.0 && r_testy > 1.0) {
        D_EXPAND(vel_x1 = d->Vc[VX1][k][j][i];,
                 vel_x2 = d->Vc[VX2][k][j-l][i];,
                 vel_x3 = d->Vc[VX3][k][j][i];)          
        vel_mag = sqrt(vel_x1*vel_x1 + vel_x2*vel_x2);
      } else if (r < 1.0 && r_testx > 1.0 && r_testy < 1.0) {
        D_EXPAND(vel_x1 = d->Vc[VX1][k][j][i-l];,
                 vel_x2 = d->Vc[VX2][k][j][i];,
                 vel_x3 = d->Vc[VX3][k][j][i];)          
        vel_mag = sqrt(vel_x1*vel_x1 + vel_x2*vel_x2);
      }
*/
    }
    // Upper left.
    else if (x1[i] < 0.0 && x2[j] > 0.0) {
      r_testx = D_EXPAND(x1[i-l]*x1[i-l], + x2[j]*x2[j], + x3[k]*x3[k]);
      r_testx = sqrt(r_testx);
      r_testy = D_EXPAND(x1[i]*x1[i], + x2[j+l]*x2[j+l], + x3[k]*x3[k]);
      r_testy = sqrt(r_testy);
      if (r < 1.0 && (r_testx > 1.0 || r_testy > 1.0)) {
        vel_mag = TwoDimensionalInterp(d, box, grid, i, j, k, x1, x2, r, ddr);
      }
/*
      if (r < 1.0 && r_testx > 1.0 && r_testy > 1.0) {
        // these are the cells at the boundary.
        D_EXPAND(vel_x1 = d->Vc[VX1][k][j][i-l];,
                 vel_x2 = d->Vc[VX2][k][j+l][i];,
                 vel_x3 = d->Vc[VX3][k][j][i];)          
        vel_mag = sqrt(vel_x1*vel_x1 + vel_x2*vel_x2);
      } else if (r < 1.0 && r_testx < 1.0 && r_testy > 1.0) {
        D_EXPAND(vel_x1 = d->Vc[VX1][k][j][i];,
                 vel_x2 = d->Vc[VX2][k][j+l][i];,
                 vel_x3 = d->Vc[VX3][k][j][i];)          
        vel_mag = sqrt(vel_x1*vel_x1 + vel_x2*vel_x2);
      } else if (r < 1.0 && r_testx > 1.0 && r_testy < 1.0) {
        D_EXPAND(vel_x1 = d->Vc[VX1][k][j][i-l];,
                 vel_x2 = d->Vc[VX2][k][j][i];,
                 vel_x3 = d->Vc[VX3][k][j][i];)          
        vel_mag = sqrt(vel_x1*vel_x1 + vel_x2*vel_x2);
      }
*/
    }
  }

  return fabs(vel_mag);
}


double VelocityStellarSurface3D(const Data *d, RBox *box, Grid *grid, 
                            int i, int j, int k, double *x1, double *x2, double *x3, double r){

  int l;
  double vel_x1, vel_x2, vel_x3;
  double r_testx, r_testy, r_testz, vel_mag;
  double dx, dy, dz, dr2, dr, ddr;

  dx = x1[i+1] - x1[i];
  dy = x2[j+1] - x2[j];
  dz = x3[k+1] - x3[k];
  dr2  = EXPAND(dx*dx,+dy*dy,+dz*dz);
  dr   = sqrt(dr2);
  ddr = 0.9*dr;

  for (l=1; l<2; l++){
    // Top of sphere.
    if (x1[i] > 0.0 && x2[j] > 0.0 && x3[k] > 0.0) {
      r_testx = D_EXPAND(x1[i+l]*x1[i+l], + x2[j]*x2[j], + x3[k]*x3[k]);
      r_testx = sqrt(r_testx);
      r_testy = D_EXPAND(x1[i]*x1[i], + x2[j+l]*x2[j+l], + x3[k]*x3[k]);
      r_testy = sqrt(r_testy);
      r_testz = D_EXPAND(x1[i]*x1[i], + x2[j]*x2[j], + x3[k+l]*x3[k+l]);
      r_testz = sqrt(r_testz);
      if (r < 1.0 && (r_testx > 1.0 || r_testy > 1.0 || r_testz > 1.0)) {
        vel_mag = ThreeDimensionalInterp(d, box, grid, i, j, k, x1, x2, x3, r, ddr);
      }
    }
    else if (x1[i] > 0.0 && x2[j] < 0.0 && x3[k] > 0.0) {
      r_testx = D_EXPAND(x1[i+l]*x1[i+l], + x2[j]*x2[j], + x3[k]*x3[k]);
      r_testx = sqrt(r_testx);
      r_testy = D_EXPAND(x1[i]*x1[i], + x2[j-l]*x2[j-l], + x3[k]*x3[k]);
      r_testy = sqrt(r_testy);
      r_testz = D_EXPAND(x1[i]*x1[i], + x2[j]*x2[j], + x3[k+l]*x3[k+l]);
      r_testz = sqrt(r_testz);
      if (r < 1.0 && (r_testx > 1.0 || r_testy > 1.0 || r_testz > 1.0)) {
        vel_mag = ThreeDimensionalInterp(d, box, grid, i, j, k, x1, x2, x3, r, ddr);
      }
    }
    else if (x1[i] < 0.0 && x2[j] < 0.0 && x3[k] > 0.0) {
      r_testx = D_EXPAND(x1[i-l]*x1[i-l], + x2[j]*x2[j], + x3[k]*x3[k]);
      r_testx = sqrt(r_testx);
      r_testy = D_EXPAND(x1[i]*x1[i], + x2[j-l]*x2[j-l], + x3[k]*x3[k]);
      r_testy = sqrt(r_testy);
      r_testz = D_EXPAND(x1[i]*x1[i], + x2[j]*x2[j], + x3[k+l]*x3[k+l]);
      r_testz = sqrt(r_testz);
      if (r < 1.0 && (r_testx > 1.0 || r_testy > 1.0 || r_testz > 1.0)) {
        vel_mag = ThreeDimensionalInterp(d, box, grid, i, j, k, x1, x2, x3, r, ddr);
      }
    }
    else if (x1[i] < 0.0 && x2[j] > 0.0 && x3[k] > 0.0) {
      r_testx = D_EXPAND(x1[i-l]*x1[i-l], + x2[j]*x2[j], + x3[k]*x3[k]);
      r_testx = sqrt(r_testx);
      r_testy = D_EXPAND(x1[i]*x1[i], + x2[j+l]*x2[j+l], + x3[k]*x3[k]);
      r_testy = sqrt(r_testy);
      r_testz = D_EXPAND(x1[i]*x1[i], + x2[j]*x2[j], + x3[k+l]*x3[k+l]);
      r_testz = sqrt(r_testz);
      if (r < 1.0 && (r_testx > 1.0 || r_testy > 1.0 || r_testz > 1.0)) {
        vel_mag = ThreeDimensionalInterp(d, box, grid, i, j, k, x1, x2, x3, r, ddr);
      }
    }
    // Bottom of sphere.
    if (x1[i] > 0.0 && x2[j] > 0.0 && x3[k] < 0.0) {
      r_testx = D_EXPAND(x1[i+l]*x1[i+l], + x2[j]*x2[j], + x3[k]*x3[k]);
      r_testx = sqrt(r_testx);
      r_testy = D_EXPAND(x1[i]*x1[i], + x2[j+l]*x2[j+l], + x3[k]*x3[k]);
      r_testy = sqrt(r_testy);
      r_testz = D_EXPAND(x1[i]*x1[i], + x2[j]*x2[j], + x3[k-l]*x3[k-l]);
      r_testz = sqrt(r_testz);
      if (r < 1.0 && (r_testx > 1.0 || r_testy > 1.0 || r_testz > 1.0)) {
        vel_mag = ThreeDimensionalInterp(d, box, grid, i, j, k, x1, x2, x3, r, ddr);
      }
    }
    else if (x1[i] > 0.0 && x2[j] < 0.0 && x3[k] < 0.0) {
      r_testx = D_EXPAND(x1[i+l]*x1[i+l], + x2[j]*x2[j], + x3[k]*x3[k]);
      r_testx = sqrt(r_testx);
      r_testy = D_EXPAND(x1[i]*x1[i], + x2[j-l]*x2[j-l], + x3[k]*x3[k]);
      r_testy = sqrt(r_testy);
      r_testz = D_EXPAND(x1[i]*x1[i], + x2[j]*x2[j], + x3[k-l]*x3[k-l]);
      r_testz = sqrt(r_testz);
      if (r < 1.0 && (r_testx > 1.0 || r_testy > 1.0 || r_testz > 1.0)) {
        vel_mag = ThreeDimensionalInterp(d, box, grid, i, j, k, x1, x2, x3, r, ddr);
      }
    }
    else if (x1[i] < 0.0 && x2[j] < 0.0 && x3[k] < 0.0) {
      r_testx = D_EXPAND(x1[i-l]*x1[i-l], + x2[j]*x2[j], + x3[k]*x3[k]);
      r_testx = sqrt(r_testx);
      r_testy = D_EXPAND(x1[i]*x1[i], + x2[j-l]*x2[j-l], + x3[k]*x3[k]);
      r_testy = sqrt(r_testy);
      r_testz = D_EXPAND(x1[i]*x1[i], + x2[j]*x2[j], + x3[k-l]*x3[k-l]);
      r_testz = sqrt(r_testz);
      if (r < 1.0 && (r_testx > 1.0 || r_testy > 1.0 || r_testz > 1.0)) {
        vel_mag = ThreeDimensionalInterp(d, box, grid, i, j, k, x1, x2, x3, r, ddr);
      }
    }
    else if (x1[i] < 0.0 && x2[j] > 0.0 && x3[k] < 0.0) {
      r_testx = D_EXPAND(x1[i-l]*x1[i-l], + x2[j]*x2[j], + x3[k]*x3[k]);
      r_testx = sqrt(r_testx);
      r_testy = D_EXPAND(x1[i]*x1[i], + x2[j+l]*x2[j+l], + x3[k]*x3[k]);
      r_testy = sqrt(r_testy);
      r_testz = D_EXPAND(x1[i]*x1[i], + x2[j]*x2[j], + x3[k-l]*x3[k-l]);
      r_testz = sqrt(r_testz);
      if (r < 1.0 && (r_testx > 1.0 || r_testy > 1.0 || r_testz > 1.0)) {
        vel_mag = ThreeDimensionalInterp(d, box, grid, i, j, k, x1, x2, x3, r, ddr);
      }
    }
  }
  return vel_mag;
}

