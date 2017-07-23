/*================================================================================*/
/*
	This init file initialises an isothermal parker wind around an sun like 
	star with a jutier like exoplanet with an orbital radius of 0.047 AU. 
	Due to the close proximity of the planet to its hoast star the planets 
	surface is being iradiated by intence UV radiation which induces a 
	planetary wind. The planet, like the star, has a dipolar magnetosphere 
	and hence the planetary wind is initilised in the same way as the stellar 
	wind.

	The model is described in detail in the paper:
	Matsakos T., Uribe1 A., Königl A., 2015, A&A, 578, A6
*/
/*================================================================================*/
#include "pluto.h"

void ParkerVelocity(double *parker, double cs, double v_esc, double r, 
                    double rc, double R, double RH, double P);

void Init (double *v, double x1, double x2, double x3)
{

  int i;
  double cs, rho, pre, omega, omega_fr;
  double Rs, v_esc, ga;
  double r, r2, theta, phi, rp, rp2;
  double M, B0, M_dot, T;
  double P, RH, rc;
  double kb, mp, sphere;
  double parker[3];
  double omega_orb;
  double xp, yp, zp;
  double bx, by, bz;
  double bxp, byp, bzp;
  double beta;

  // quantity | value | units

  // adiabatic index
  g_smallPressure = 1.0e-5; // Small value for pressure fix.
  g_gamma = 1.05;

  // general quatities
  kb        = 1.38064852e-16; // [erg K-1]
  mp        = 1.67262171e-24; // [g]

  // Planet & Star properties
  T = g_inputParam[T_star];                                           
  M = g_inputParam[M_star]*CONST_Msun/M0;                                    
  Rs = 1.0;//6.955e+10*(1.0/UNIT_LENGTH);                         
  B0 = g_inputParam[B_star]/Bc;                                             
  RH = g_inputParam[RHO_star]/UNIT_DENSITY;                              
  v_esc = sqrt(2.0*UNIT_G*M/Rs); // Star escape speed
  cs = sqrt((2.0*kb*T)/mp)/UNIT_VELOCITY;
  omega = g_inputParam[Omega]*2.67e-6*t0;
  P = cs*cs*RH/g_gamma;

  // find radius from star and planet
  r2       = EXPAND(x1*x1,+x2*x2,+x3*x3);
  r        = sqrt(r2);

  // coordinate transforms
  theta    = acos(x3/r);
  phi      = atan2(x2,x1);

  // The critical point
  rc = UNIT_G*M/(2.0*cs*cs);

  beta = 0.0;
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


  // Set the density, pressure, velocity of the stellar wind 
  // and stellar interiar.
  if (r <= 0.5*Rs) {

    EXPAND(v[VX1] = 0.0;,
           v[VX2] = 0.0;,
           v[VX3] = 0.0;)
    v[PRS] = P + (2.0/3.0)*CONST_PI*UNIT_G*RH*RH*(Rs*Rs-r*r);
    v[RHO] = RH;

#if PHYSICS == MHD 
#if BACKGROUND_FIELD == NO
#if DIMENSIONS == 2
    bx = 0.0;
    by = B0*16.0;
    bxp = bx*cos(beta) + by*sin(beta);;
    byp = -bx*sin(beta) + by*cos(beta);
    EXPAND(v[BX1] = bxp;,
           v[BX2] = byp;,
           v[BX3] = 0.0;)
#endif
#if DIMENSIONS == 3
    bx = 0.0;
    by = 0.0;
    bz = 16.0*B0;
    bxp = bx*cos(beta) + bz*sin(beta);
    byp = by;
    bzp = -bx*sin(beta) + bz*cos(beta);
    EXPAND(v[BX1] = bxp;,
           v[BX2] = byp;,
           v[BX3] = bzp;)
#endif 
#endif
#endif

  } else if(r <= Rs && r > 0.5*Rs){

    EXPAND(v[VX1] = 0.0;,
           v[VX2] = 0.0;,
           v[VX3] = 0.0;)
    v[PRS] = P + (2.0/3.0)*CONST_PI*UNIT_G*RH*RH*(Rs*Rs-r*r);
    v[RHO] = RH;

#if PHYSICS == MHD
#if BACKGROUND_FIELD == NO
#if DIMENSIONS == 2
    bx = 3.0*xp*yp*B0*pow(rp,-5);
    by = (3.0*pow(yp,2)-pow(rp,2))*B0*pow(rp,-5);
    bxp = bx*cos(beta) + by*sin(beta);;
    byp = -bx*sin(beta) + by*cos(beta);
    EXPAND(v[BX1] = bxp;,
           v[BX2] = byp;,
           v[BX3] = 0.0;)
#endif
#if DIMENSIONS == 3
    bx = 3.0*xp*zp*B0*pow(rp,-5);
    by = 3.0*yp*zp*B0*pow(rp,-5);
    bz = (3.0*pow(zp,2)-pow(rp,2))*B0*pow(rp,-5);
    bxp = bx*cos(beta) + bz*sin(beta);
    byp = by;
    bzp = -bx*sin(beta) + bz*cos(beta);
    EXPAND(v[BX1] = bxp;,
           v[BX2] = byp;,
           v[BX3] = bzp;)
#endif
#endif
#endif

  } else if (r > Rs){

    ParkerVelocity(parker, cs, v_esc, r, rc, Rs, RH, P);
    EXPAND(v[VX1] = sin(theta)*(parker[0]*cos(phi)+sin(phi)*r*(omega+omega));,
           v[VX2] = sin(theta)*(parker[0]*sin(phi)-cos(phi)*r*(omega+omega));,
           v[VX3] = parker[0]*cos(theta);)
    v[PRS] = parker[1];
    v[RHO] = parker[2];

#if PHYSICS == MHD
#if BACKGROUND_FIELD == NO
#if DIMENSIONS == 2
    bx = 3.0*xp*yp*B0*pow(rp,-5);
    by = (3.0*pow(yp,2)-pow(rp,2))*B0*pow(rp,-5);
    bxp = bx*cos(beta) + by*sin(beta);;
    byp = -bx*sin(beta) + by*cos(beta);
    EXPAND(v[BX1] = bxp;,
           v[BX2] = byp;,
           v[BX3] = 0.0;)
#endif
#if DIMENSIONS == 3
    bx = 3.0*xp*zp*B0*pow(rp,-5);
    by = 3.0*yp*zp*B0*pow(rp,-5);
    bz = (3.0*pow(zp,2)-pow(rp,2))*B0*pow(rp,-5);
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

  int nv;
  for (nv = 0; nv < NVAR; nv++) {
    if (isnan(v[nv])) {
      printf("x1=%e, x2=%e, x3=%e \n", x1, x2, x3);
      printf("v[RHO]=%e \n", v[RHO]);
      printf("v[PRS]=%e \n", v[PRS]);
      printf("v[VX1]=%e \n", v[VX1]);
      printf("v[VX2]=%e \n", v[VX2]);
      printf("v[VX3]=%e \n", v[VX3]);
      printf("v[BX1]=%e \n", v[BX1]);
      printf("v[BX2]=%e \n", v[BX2]);
      printf("v[BX3]=%e \n", v[BX3]);
    }
  }

}  

/*================================================================================*/

/*================================================================================*/
void Analysis (const Data *d, Grid *grid)
{
}
/*================================================================================*/

/*================================================================================*/
#if BACKGROUND_FIELD == YES
void BackgroundField (double x1, double x2, double x3, double *B0)                                                    
{                                                                                                                     

  double Bs, r2, r, R, rp, rp2;
  double xp, yp, zp;
  double bx, by, bz;
  double bxp, byp, bzp;
  double beta;
                                                                                                          
  R  = 1.0;//g_inputParam[R_star]*CONST_Rsun/UNIT_LENGTH;
  Bs = g_inputParam[B_star]/Bc;
  r2 = EXPAND(x1*x1,+x2*x2,+x3*x3);
  r  = sqrt(r2);

  beta = 0.0;
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

  if (r <= 0.5*R){
#if DIMENSIONS == 2
    bx = 0.0;
    by = B0*16.0;
    bxp = bx*cos(beta) + by*sin(beta);;
    byp = -bx*sin(beta) + by*cos(beta);
    bzp = 0.0;
#endif
#if DIMENSIONS == 3
    bx = 0.0;
    by = 0.0;
    bz = 16.0*B0;
    bxp = bx*cos(beta) + bz*sin(beta);
    byp = by;
    bzp = -bx*sin(beta) + bz*cos(beta);
#endif 
  } else if (r > 0.5*R) {
#if DIMENSIONS == 2
    bx = 3.0*xp*yp*B0*pow(rp,-5);
    by = (3.0*pow(yp,2)-pow(rp,2))*B0*pow(rp,-5);
    bxp = bx*cos(beta) + by*sin(beta);;
    byp = -bx*sin(beta) + by*cos(beta);
    bzp
#endif
#if DIMENSIONS == 3
    bx = 3.0*xp*zp*B0*pow(rp,-5);
    by = 3.0*yp*zp*B0*pow(rp,-5);
    bz = (3.0*pow(zp,2)-pow(rp,2))*B0*pow(rp,-5);
    bxp = bx*cos(beta) + bz*sin(beta);
    byp = by;
    bzp = -bx*sin(beta) + bz*cos(beta);
#endif

  }

  B0[0] = bxp;  
  B0[1] = byp;
  B0[2] = bzp;       

}
#endif
/*================================================================================*/


/*================================================================================*/
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid) 
{
  int i, j, k, o;
  double *x1, *x2, *x3;
  double r, cs, omega;
  double R, v_esc, ga;
  double r2, theta, phi, rp, rp2;
  double M, B0, M_dot, T;
  double P, RH;
  double rc; 
  double kb, mp;
  double parker[3];
  double omega_orb;
  double xp, yp, zp;
  double bx, by, bz;
  double bxp, byp, bzp;
  double beta;
 
  // - quantity | value | units.

  // - adiabatic index.
  g_gamma   = 1.05;

  // general quatities.
  kb        = 1.38064852e-16; // [erg K-1]
  mp        = 1.67262171e-24; // [g]

  // Planet & Star properties.
  T = g_inputParam[T_star];                                           
  M = g_inputParam[M_star]*CONST_Msun/M0;                                    
  R = 1.0;//g_inputParam[R_star]*CONST_Rsun/UNIT_LENGTH;                         
  B0 = g_inputParam[B_star]/Bc;                                             
  RH = g_inputParam[RHO_star]/UNIT_DENSITY;                              
  v_esc    = sqrt(2.0*UNIT_G*M/R);
  cs       = sqrt((2.0*kb*T)/mp)/UNIT_VELOCITY;
  P        = cs*cs*RH/g_gamma;
  rc = UNIT_G*M/(2.0*cs*cs);
  omega = g_inputParam[Omega]*2.67e-6*t0;
  beta = 0.0;
  beta *= 0.0174532925;

  x1 = grid[IDIR].xgc;
  x2 = grid[JDIR].xgc;
  x3 = grid[KDIR].xgc;


  if (side == 0) {
    TOT_LOOP(k,j,i){	

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

      /* - Radial quantities - */
      r2 = EXPAND(x1[i]*x1[i],+x2[j]*x2[j],+x3[k]*x3[k]);
      r  = sqrt(r2);

      /* - coordinate transforms - */
      theta = acos(x3[k]/r);
      phi   = atan2(x2[j],x1[i]);


      // Set the density, pressure, velocity of the stellar wind 
      // and stellar interiar.
      if (r <= 0.5*R){

        EXPAND(d->Vc[VX1][k][j][i] = 0.0;,
               d->Vc[VX2][k][j][i] = 0.0;,
               d->Vc[VX3][k][j][i] = 0.0;)
        d->Vc[PRS][k][j][i] = P + (2.0/3.0)*CONST_PI*UNIT_G*RH*RH*(R*R-r*r);
        d->Vc[RHO][k][j][i] = RH;

#if PHYSICS == MHD
#if BACKGROUND_FIELD == NO
#if DIMENSIONS == 2
        bx = 0.0;
        by = 16.0*B0;
        bxp = bx*cos(beta) + by*sin(beta);
        byp = -bx*sin(beta) + by*cos(beta);
        EXPAND(d->Vc[BX1][k][j][i] = bxp;,
               d->Vc[BX2][k][j][i] = byp;,
               d->Vc[BX3][k][j][i] = 0.0;)
#endif
#if DIMENSIONS == 3
        bx = 0.0;
        by = 0.0;
        bz = 16.0*B0;
        bxp = bx*cos(beta) + bz*sin(beta);
        byp = by;
        bzp = -bx*sin(beta) + bz*cos(beta);
        EXPAND(d->Vc[BX1][k][j][i] = bxp;,
               d->Vc[BX2][k][j][i] = byp;,
               d->Vc[BX3][k][j][i] = bzp;)
#endif 
#endif
#endif

        d->flag[k][j][i]   |= FLAG_INTERNAL_BOUNDARY;

      } else if (r > 0.5*R && r <= R) {

        EXPAND(d->Vc[VX1][k][j][i] = 0.0;,
               d->Vc[VX2][k][j][i] = 0.0;,
               d->Vc[VX3][k][j][i] = 0.0;)
        d->Vc[PRS][k][j][i] = P + (2.0/3.0)*CONST_PI*UNIT_G*RH*RH*(R*R-r*r);
        d->Vc[RHO][k][j][i] = RH;

#if PHYSICS == MHD
#if BACKGROUND_FIELD == NO
#if DIMENSIONS == 2	
        bx = 3.0*xp*yp*B0*pow(rp,-5);
        by = (3.0*pow(yp,2)-pow(rp,2))*B0*pow(rp,-5);
        bxp = bx*cos(beta) + by*sin(beta);
        byp = -bx*sin(beta) + by*cos(beta);
        EXPAND(d->Vc[BX1][k][j][i] = bxp;,
               d->Vc[BX2][k][j][i] = byp;,
               d->Vc[BX3][k][j][i] = 0.0;)
#endif
#if DIMENSIONS == 3
        bx = 3.0*xp*zp*B0*pow(rp,-5);
        by = 3.0*yp*zp*B0*pow(rp,-5);
        bz = (3.0*pow(zp,2)-pow(rp,2))*B0*pow(rp,-5);
        bxp = bx*cos(beta) + bz*sin(beta);
        byp = by;
        bzp = -bx*sin(beta) + bz*cos(beta);
        EXPAND(d->Vc[BX1][k][j][i] = bxp;,
               d->Vc[BX2][k][j][i] = byp;,
               d->Vc[BX3][k][j][i] = bzp;)
#endif
#endif
#endif
        d->flag[k][j][i]   |= FLAG_INTERNAL_BOUNDARY;

      } else if (r > R && r <= 1.5*R){

        ParkerVelocity(parker, cs, v_esc, r, rc, R, RH, P);
        EXPAND(d->Vc[VX1][k][j][i] = sin(theta)*(parker[0]*cos(phi)+sin(phi)*
                                       r*(omega+omega));,
               d->Vc[VX2][k][j][i] = sin(theta)*(parker[0]*sin(phi)-cos(phi)*
                                       r*(omega+omega));,
               d->Vc[VX3][k][j][i] = parker[0]*cos(theta);)
        d->Vc[PRS][k][j][i] = parker[1];
        d->Vc[RHO][k][j][i] = parker[2];
        d->flag[k][j][i]   |= FLAG_INTERNAL_BOUNDARY;

      }

      if (d->Vc[PRS][k][j][i] < g_smallPressure){
        d->Vc[PRS][k][j][i] = g_smallPressure;
      }

      int nv;
      for (nv = 0; nv < NVAR; nv++) {
        if (isnan(d->Vc[nv][k][j][i])) {
          printf("x1=%e, x2=%e, x3=%e \n", x1[i], x2[j], x3[k]);
          printf("xp=%e, yp=%e, zp=%e \n", xp, yp, zp);
          printf("bx=%e, by=%e, bz=%e \n", bx, by, bz);
          printf("bxp=%e, byp=%e, bzp=%e \n", bxp, byp, bzp);
          printf("B0=%e, beta=%e \n", B0, beta);
          printf("d->Vc[RHO][k][j][i]=%e \n", d->Vc[RHO][k][j][i]);
          printf("d->Vc[PRS][k][j][i]=%e \n", d->Vc[PRS][k][j][i]);
          printf("d->Vc[VX1][k][j][i]=%e \n", d->Vc[VX1][k][j][i]);
          printf("d->Vc[VX2][k][j][i]=%e \n", d->Vc[VX2][k][j][i]);
          printf("d->Vc[VX3][k][j][i]=%e \n", d->Vc[VX3][k][j][i]);
          printf("d->Vc[BX1][k][j][i]=%e \n", d->Vc[BX1][k][j][i]);
          printf("d->Vc[BX2][k][j][i]=%e \n", d->Vc[BX2][k][j][i]);
          printf("d->Vc[BX3][k][j][i]=%e \n", d->Vc[BX3][k][j][i]);
        }
      }

    }
  }
}
/*================================================================================*/

/*================================================================================*/
#if BODY_FORCE != NO
void BodyForceVector(double *v, double *g, double x1, double x2, double x3)
{

  double M, r;
  double r2, RH, R, omega;
  double Fcentr_s_x1, Fcentr_s_x2, Fcor_s_x1, Fcor_s_x2, Fin_s_x1, Fin_s_x2;
  double Fcentr_p_x1, Fcentr_p_x2, Fcor_p_x1, Fcor_p_x2, Fin_p_x1, Fin_p_x2;
  double Fin_x1, Fin_x2, gs, gp, gs_in, gp_in;
  double omega_orb;

  M = g_inputParam[M_star]*CONST_Msun/M0;
  RH = g_inputParam[RHO_star]/UNIT_DENSITY;
  R = 1.0;//g_inputParam[R_star]*CONST_Rsun/UNIT_LENGTH;

  // Rotational frequency (orbit and frame).
  omega = g_inputParam[Omega]*2.67e-6*t0;
  // Distance from star.
  r2 = EXPAND(x1*x1, + x2*x2, + x3*x3);
  r = sqrt(r2);
  // Gravity outside bodies.
  gs = -UNIT_G*M/r/r;
  // Gravity inside bodies.
  gs_in = -(4.0/3.0)*CONST_PI*UNIT_G*RH;     

  // Coriolis and centrifugal forces.
#if DIMENSIONS == 3
  Fin_x1 = omega*omega*x1 + 2.0*omega*v[VX2];
  Fin_x2 = omega*omega*x2 - 2.0*omega*v[VX1];
#endif
#if DIMENSIONS == 2 
  Fin_x1 = 0.0;
  Fin_x2 = 0.0;
#endif

  if (r > R){ // External gravity + centrifugal + coriolis.
    g[IDIR] = gs*x1/r + Fin_x1;
    g[JDIR] = gs*x2/r + Fin_x2;
    g[KDIR] = gs*x3/r;
  } else if (r < R) { // Star interal gravity.
    g[IDIR] = gs_in*x1;
    g[JDIR] = gs_in*x2;
    g[KDIR] = gs_in*x3;
  }
}
#endif
/*================================================================================*/

/*================================================================================*/
/* 
  In order to obtain the velocity as a function of radial 
  distance from either the star or planet it is nessacery 
  to solve parkers equation and find its roots. There are 
  two roots one for the subsonic and one for the supersonic 
  regions. the transition between the two rigines is given 
  by the critical point rc (see above expressions). The 
  roots are found using the Newton raphson tschnique.
*/
/*================================================================================*/
void ParkerVelocity(double *parker, double cs, double v_esc, double r, 
                    double rc, double R, double RH, double P)
{

  int o;
  double lambda, eta, b;
  double psi, step, fn, dfn;

  lambda = 0.5*pow(v_esc/cs,2);
  eta = r/R;
  b = -3.0 - 4.0*log(lambda/2.0) + 4.0*log(eta) + 2.0*(lambda/eta);
  o = 0;

  if (r <= rc) {
    psi = 2.e-8;
  } else {
    psi = 2.5;
  }
  do {
    fn = -psi + log(psi) + b;
    dfn = -1.0 + (1.0/psi);
    step = fn/dfn;
    psi = psi - step;
    o++;
  } while (fabs(step) > 1.0e-8 && o < 100000);
  if (isnan(psi)){
    psi = 0.0;
  }

  parker[0] = cs*sqrt(psi);
  parker[1] = P*exp(lambda*(R/r-1.0)-0.5*pow(parker[0]/cs,2));
  parker[2] = (RH/P)*parker[1];
}

/*================================================================================*/

