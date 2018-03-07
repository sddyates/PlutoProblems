/*============================================================================*/
/*
	This init file initialises an isothermal parker wind around an sun like 
	star with a jutier like exoplanet with an orbital radius of 0.047 AU. 
	Due to the close proximity of the planet to its hoast star the planets 
	surface is being iradiated by intence UV radiation which induces a 
	planetary wind. The planet, like the star, has a dipolar magnetosphere 
	and hence the planetary wind is initilised in the same way as the 
        stellar wind.

	The model is described in detail in the paper:
	Matsakos T., Uribe1 A., Königl A., 2015, A&A, 578, A6
*/
/*============================================================================*/
#include "pluto.h"

void ParkerVelocity(double *, double, double, double, 
                    double, double, double, double);

void planet(double *, double *, double, double, double);

void VelCarToSph(double *, double, double, double, double, double, double, double);

void Init(double *v, double x1, double x2, double x3)
{


  int i;
  double css, csp, rhos, rhop, Mp, omegas, omegap, omega_fr;
  double Rs, Rp, a, B0p, v_escs, v_escp, ga;
  double rp2, rp, thetap, phip;
  double rs, rs2, thetas, phis;
  double Ms, B0s, vs, M_dots, M_dotp, Ts, Tp, omega_orb;
  double Ps, Pp, RHs, RHp, rcs, rcp;
  double kb, mp, sphere;
  double parker[3];
  double br, btheta, bphi, r, theta, phi;

  double x, y, z;
  double xp, yp, zp;
  double vx1, vx2, vx3;
  double vr, vtheta, vphi;
  double vx, vy, vz;

  double bpx, bpy, bpz;
  double bsx, bsy, bsz;
  double bsx1, bsx2, bsx3;
  double bpx1, bpx2, bpx3;

  // quantity | value | units

  // adiabatic index

  g_smallPressure = 1.0e-5; // Small value for pressure fix.
  g_gamma = 1.05;

  // general quatities

  a = g_inputParam[seperation]*214.9394693836;
  kb        = 1.38064852e-16; // [erg K-1]
  mp        = 1.67262171e-24; // [g]

  // Planet & Star properties

  Tp        = g_inputParam[planet_temperature];                                         
  Ts        = g_inputParam[star_temperature];                                           

  Mp        = g_inputParam[planet_mass]*0.0009543*CONST_Msun/UNIT_MASS;                                  
  Ms        = g_inputParam[star_mass]*CONST_Msun/UNIT_MASS;                                    

  Rp        = g_inputParam[planet_radius]*0.10045*CONST_Rsun/UNIT_LENGTH;                       
  Rs        = g_inputParam[star_radius]*CONST_Rsun/UNIT_LENGTH;                         

  B0p       = g_inputParam[planet_Bfield]/UNIT_B;                                           
  B0s       = g_inputParam[star_Bfield]/UNIT_B;                                             

  RHp       = g_inputParam[planet_surface_rho]/UNIT_DENSITY;                            
  RHs       = g_inputParam[star_surface_rho]/UNIT_DENSITY;                              

  v_escp    = sqrt(2.0*UNIT_G*Mp/Rp); // Planet escape speed
  v_escs    = sqrt(2.0*UNIT_G*Ms/Rs); // Star escape speed
  csp       = sqrt((2.0*kb*Tp)/mp)/UNIT_VELOCITY;
  css       = sqrt((2.0*kb*Ts)/mp)/UNIT_VELOCITY;

  omega_orb = sqrt(UNIT_G*Ms/pow(a,3));
  omegas    = omega_orb;      
  omegap    = omega_orb;
  omega_fr  = omega_orb;

  Pp        = csp*csp*RHp/g_gamma;                                                           
  Ps        = css*css*RHs/g_gamma;

  // The critical point
  rcs = UNIT_G*Ms/(2.0*css*css);
  rcp = UNIT_G*Mp/(2.0*csp*csp);

#if ROTATING_FRAME == YES
  g_OmegaZ = omega_orb;
#endif

  /* - The switch to Cartesian coordinates. - */

  // Convert to Cartesian centered on the star.
  x = x1*sin(x2)*cos(x3);
  y = x1*sin(x2)*sin(x3);
  z = x1*cos(x2);

  // Convert to Cartesian centered on the planet.
  xp = x - a;
  yp = y;
  zp = z;

  // Convert to spherical centered on the planet. 
  rp2 = EXPAND(xp*xp, + yp*yp, + zp*zp);
  rp = sqrt(rp2);
  thetap = acos(zp/rp);
  phip = atan2(yp, xp);

  // Set the density, pressure, velocity of the stellar wind 
  // and planetary interiar.
  if (rp <= Rp) {

    vx = 0.0;
    vy = 0.0;
    vz = 0.0;
    v[PRS] = Pp + (2.0/3.0)*CONST_PI*UNIT_G*RHp*RHp*(Rp*Rp-rp*rp);
    v[RHO] = RHp;

  } else if (rp > Rp && rp <= 10.0*Rp) {

    ParkerVelocity(parker, csp, v_escp, rp, rcp, Rp, RHp, Pp);
    v[PRS] = parker[1];
    v[RHO] = parker[2];

    // Cartesian velocity components.
    vx = sin(thetap)*(parker[0]*cos(phip) + 
         sin(x3)*x1*omega_fr - 
         sin(phip)*rp*omegap);

    vy = sin(thetap)*(parker[0]*sin(phip) - 
         cos(x3)*x1*omega_fr + 
         a*omega_orb + 
         cos(phip)*rp*omegap);

    vz = parker[0]*cos(thetap);

    EXPAND(v[VX1] = vx*sin(x2)*cos(x3) + vy*sin(x2)*sin(x3) + vz*cos(x2);,
           v[VX2] = vx*cos(x2)*cos(x3) + vy*cos(x2)*sin(x3) - vz*sin(x2);,
           v[VX3] = -vx*sin(x3)        + vy*cos(x3);)
  }

  if (rp > Rp && rp <= 1.5*Rp){
    v[TRC] = 1.0;
  }

  // Set the density, pressure, velocity of the stellar wind 
  // and stellar interiar.
  if(x1 <= Rs){

    EXPAND(v[VX1] = 0.0;,
           v[VX2] = 0.0;,
           v[VX3] = 0.0;)
    v[PRS] = Ps + (2.0/3.0)*CONST_PI*UNIT_G*RHs*RHs*(Rs*Rs-x1*x1);
    v[RHO] = RHs;

  } else if (x1 > Rs && rp > 10.0*Rp){

    ParkerVelocity(parker, css, v_escs, x1, rcs, Rs, RHs, Ps);

    vx = sin(x2)*(parker[0]*cos(x3) + sin(x3)*x1*(omega_fr + omegas));
    vy = sin(x2)*(parker[0]*sin(x3) - cos(x3)*x1*(omega_fr + omegas));
    vz = parker[0]*cos(x2);

    EXPAND(v[VX1] = vx*sin(x2)*cos(x3) + vy*sin(x2)*sin(x3) + vz*cos(x2);,
           v[VX2] = vx*cos(x2)*cos(x3) + vy*cos(x2)*sin(x3) - vz*sin(x2);,
           v[VX3] = -vx*sin(x3)        + vy*cos(x3);)
    v[PRS] = parker[1];
    v[RHO] = parker[2];
  }
  
#if PHYSICS == MHD                                   
#if BACKGROUND_FIELD == NO

  // Define Planet magnetic field in Cartesian.
  if (rp <= 0.5*Rp){

    bpx = 0.0;
    bpy = 0.0;
    bpz = 16.0*B0p;

    // Convert To spherical basis form Cartesian.
    bpx1 = bpx*sin(x2)*cos(x3) + bpy*sin(x2)*sin(x3) + bpz*cos(x2);
    bpx2 = bpx*cos(x2)*cos(x3) + bpy*cos(x2)*sin(x3) - bpz*sin(x2);
    bpx3 = -bpx*sin(x3)        + bpy*cos(x3);

    EXPAND(v[BX1] = bpx1;,            
           v[BX2] = bpx2;,                 
           v[BX3] = bpx3;) 

  } else if (rp > 0.5*Rp && rp <= Rp) {

    bpx = 3.0*xp*zp*B0p*pow(Rp, 3)*pow(rp, -5);
    bpy = 3.0*zp*yp*B0p*pow(Rp, 3)*pow(rp, -5);
    bpz = (3.0*zp*zp - rp*rp)*B0p*pow(Rp, 3)*pow(rp, -5);

    // Convert To spherical basis form Cartesian.
    bpx1 = bpx*sin(x2)*cos(x3) + bpy*sin(x2)*sin(x3) + bpz*cos(x2);
    bpx2 = bpx*cos(x2)*cos(x3) + bpy*cos(x2)*sin(x3) - bpz*sin(x2);
    bpx3 = -bpx*sin(x3)        + bpy*cos(x3);

    EXPAND(v[BX1] = bpx1;,            
           v[BX2] = bpx2;,                 
           v[BX3] = bpx3;) 

  } else if (rp > Rp) {

    bpx = 3.0*xp*zp*B0p*pow(Rp, 3)*pow(rp, -5);
    bpy = 3.0*zp*yp*B0p*pow(Rp, 3)*pow(rp, -5);
    bpz = (3.0*zp*zp - rp*rp)*B0p*pow(Rp, 3)*pow(rp, -5);

    // Convert To spherical basis form Cartesian.
    bpx1 = bpx*sin(x2)*cos(x3) + bpy*sin(x2)*sin(x3) + bpz*cos(x2);
    bpx2 = bpx*cos(x2)*cos(x3) + bpy*cos(x2)*sin(x3) - bpz*sin(x2);
    bpx3 = -bpx*sin(x3)        + bpy*cos(x3);

    bsx = 3.0*x*z*B0s*pow(Rs, 3)*pow(x1, -5);
    bsy = 3.0*z*y*B0s*pow(Rs, 3)*pow(x1, -5);
    bsz = (3.0*z*z - x1*x1)*B0s*pow(Rs, 3)*pow(x1, -5);

    bsx1 = bsx*sin(x2)*cos(x3) + bsy*sin(x2)*sin(x3) + bsz*cos(x2);
    bsx2 = bsx*cos(x2)*cos(x3) + bsy*cos(x2)*sin(x3) - bsz*sin(x2);
    bsx3 = -bsx*sin(x3)        + bsy*cos(x3);

    EXPAND(v[BX1] = bpx1 + bsx1;,            
           v[BX2] = bpx2 + bsx2;,                 
           v[BX3] = bpx3 + bsx3;) 
  }
    
#endif
#if BACKGROUND_FIELD == YES
  EXPAND(v[BX1] = 0.0;,            
         v[BX2] = 0.0;,                 
         v[BX3] = 0.0;) 
#endif
#endif

}  
/*============================================================================*/

/*============================================================================*/
void Analysis (const Data *d, Grid *grid)
{

}
/*============================================================================*/

/*================================================================================*/
#if BACKGROUND_FIELD == YES
void BackgroundField (double x1, double x2, double x3, double *B0)                                                    
{                                                                                                                     

  double a;
  double x, y, z;
  double xp, yp, zp;
  double rp, rp2;

  double bpx, bpy, bpz;
  double bsx, bsy, bsz;
  double bpx1, bpx2, bpx3;
  double bsx1, bsx2, bsx3;

  double Rp, B0p;
  double Rs, B0s;

  a = g_inputParam[seperation]*214.9394693836;
  Rp = g_inputParam[planet_radius]*0.10045*CONST_Rsun/UNIT_LENGTH;                       
  Rs = g_inputParam[star_radius]*CONST_Rsun/UNIT_LENGTH;                         

  B0p = g_inputParam[planet_Bfield]/UNIT_B;                                           
  B0s = g_inputParam[star_Bfield]/UNIT_B;                                             


  // Convert to Cartesian centered on the star.
  x = x1*sin(x2)*cos(x3);
  y = x1*sin(x2)*sin(x3);
  z = x1*cos(x2);

  // Convert to Cartesian centered on the planet.
  xp = x - a;
  yp = y;
  zp = z;
  rp2 = EXPAND(xp*xp, + yp*yp, + zp*zp);
  rp = sqrt(rp2);
                                                                                                         
  // Define Planet magnetic field in Cartesian.
  if (rp <= 0.5*Rp){

    bpx = 0.0;
    bpy = 0.0;
    bpz = 16.0*B0p;

    // Convert To spherical basis form Cartesian.
    bpx1 = bpx*sin(x2)*cos(x3) + bpy*sin(x2)*sin(x3) + bpz*cos(x2);
    bpx2 = bpx*cos(x2)*cos(x3) + bpy*cos(x2)*sin(x3) - bpz*sin(x2);
    bpx3 = -bpx*sin(x3)        + bpy*cos(x3);

    B0[0] = bpx1;            
    B0[1] = bpx2;                 
    B0[2] = bpx3;

  } else if (rp > 0.5*Rp && rp <= Rp) {

    bpx = 3.0*xp*zp*B0p*pow(Rp, 3)*pow(rp, -5);
    bpy = 3.0*zp*yp*B0p*pow(Rp, 3)*pow(rp, -5);
    bpz = (3.0*zp*zp - rp*rp)*B0p*pow(Rp, 3)*pow(rp, -5);

    // Convert To spherical basis form Cartesian.
    bpx1 = bpx*sin(x2)*cos(x3) + bpy*sin(x2)*sin(x3) + bpz*cos(x2);
    bpx2 = bpx*cos(x2)*cos(x3) + bpy*cos(x2)*sin(x3) - bpz*sin(x2);
    bpx3 = -bpx*sin(x3)        + bpy*cos(x3);

    B0[0] = bpx1;            
    B0[1] = bpx2;                 
    B0[2] = bpx3;

  } else if (rp > Rp) {

    bpx = 3.0*xp*zp*B0p*pow(Rp, 3)*pow(rp, -5);
    bpy = 3.0*zp*yp*B0p*pow(Rp, 3)*pow(rp, -5);
    bpz = (3.0*zp*zp - rp*rp)*B0p*pow(Rp, 3)*pow(rp, -5);

    // Convert To spherical basis form Cartesian.
    bpx1 = bpx*sin(x2)*cos(x3) + bpy*sin(x2)*sin(x3) + bpz*cos(x2);
    bpx2 = bpx*cos(x2)*cos(x3) + bpy*cos(x2)*sin(x3) - bpz*sin(x2);
    bpx3 = -bpx*sin(x3)        + bpy*cos(x3);

    bsx = 3.0*x*z*B0s*pow(Rs, 3)*pow(x1, -5);
    bsy = 3.0*z*y*B0s*pow(Rs, 3)*pow(x1, -5);
    bsz = (3.0*z*z - x1*x1)*B0s*pow(Rs, 3)*pow(x1, -5);

    bsx1 = bsx*sin(x2)*cos(x3) + bsy*sin(x2)*sin(x3) + bsz*cos(x2);
    bsx2 = bsx*cos(x2)*cos(x3) + bsy*cos(x2)*sin(x3) - bsz*sin(x2);
    bsx3 = -bsx*sin(x3)        + bsy*cos(x3);

    B0[0] = bpx1 + bsx1;
    B0[1] = bpx2 + bsx2;
    B0[2] = bpx3 + bsx3; 

  } else {
    print("Error at x1=%f, x2=%f, x3=%f \n", x1, x2, x3);
  }

}
#endif
/*================================================================================*/


/*================================================================================*/
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid) 
{

  int i, j, k, o, ghost;

  double mu, a, omegas, omegap, omega_fr, omega_orb;

  double Ms, Ts, Ps, RHs, Rs, rcs, B0s, css, v_escs; 
  double Mp, Tp, Pp, RHp, Rp, rcp, B0p, csp, v_escp;

  double kb, mp;
  double parker[3];

  double rp, rp2, thetap, phip;
  double x, y, z;
  double xp, yp, zp;
  double vx, vy, vz;
  double bpx1, bpx2, bpx3;
  double bpx, bpy, bpz;

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

  // quantity | value | units

  // adiabatic index

  g_smallPressure = 1.0e-5; // Small value for pressure fix.
  g_gamma = 1.05;

  // general quatities

  a = g_inputParam[seperation]*214.9394693836;
  kb        = 1.38064852e-16; // [erg K-1]
  mp        = 1.67262171e-24; // [g]

  // Planet & Star properties

  Tp        = g_inputParam[planet_temperature];                                         
  Ts        = g_inputParam[star_temperature];                                           

  Mp        = g_inputParam[planet_mass]*0.0009543*CONST_Msun/UNIT_MASS;                                  
  Ms        = g_inputParam[star_mass]*CONST_Msun/UNIT_MASS;                                    

  Rp        = g_inputParam[planet_radius]*0.10045*CONST_Rsun/UNIT_LENGTH;                       
  Rs        = g_inputParam[star_radius]*CONST_Rsun/UNIT_LENGTH;                         

  B0p       = g_inputParam[planet_Bfield]/UNIT_B;                                           
  B0s       = g_inputParam[star_Bfield]/UNIT_B;                                             

  RHp       = g_inputParam[planet_surface_rho]/UNIT_DENSITY;                            
  RHs       = g_inputParam[star_surface_rho]/UNIT_DENSITY;                              

  v_escp    = sqrt(2.0*UNIT_G*Mp/Rp); // Planet escape speed
  v_escs    = sqrt(2.0*UNIT_G*Ms/Rs); // Star escape speed
  csp       = sqrt((2.0*kb*Tp)/mp)/UNIT_VELOCITY;
  css       = sqrt((2.0*kb*Ts)/mp)/UNIT_VELOCITY;

  omega_orb = sqrt(UNIT_G*Ms/pow(a,3));
  omegas    = omega_orb;      
  omegap    = omega_orb;
  omega_fr  = omega_orb;

  Pp        = csp*csp*RHp/g_gamma;                                                           
  Ps        = css*css*RHs/g_gamma;

  rcs = UNIT_G*Ms/(2.0*css*css);
  rcp = UNIT_G*Mp/(2.0*csp*csp);

  x1 = grid[IDIR].xgc;
  x2 = grid[JDIR].xgc;
  x3 = grid[KDIR].xgc;

  ghost = (NX1_TOT - NX1)/2;

  if(side == X1_BEG){
    if(box->vpos == CENTER){
      BOX_LOOP(box,k,j,i){ 
  
        ParkerVelocity(parker, css, v_escs, Rs, rcs, Rs, RHs, Ps);
        
        EXPAND(vx1[k][j][i] = parker[0];,
               vx2[k][j][i] = 0.0;,
               vx3[k][j][i] = 0.0;)


        prs[k][j][i] = parker[1];          
        rho[k][j][i] = parker[2]; 

#if PHYSICS == MHD   
#if BACKGROUND_FIELD == NO
       EXPAND(bx1[k][j][i] = B0s*pow(x1[i], -3)*cos(x2[j]);,            
               bx2[k][j][i] = (B0s/2.0)*pow(x1[i], -3)*sin(x2[j]);,                 
               bx3[k][j][i] = 0.0;)                                                   
#endif
#if BACKGROUND_FIELD == YES
        EXPAND(bx1[k][j][i] = 0.0;,
               bx2[k][j][i] = 0.0;,
               bx3[k][j][i] = 0.0;)
#endif
#endif                                                            


      }
    }
  }

  if (side == 0) {
    TOT_LOOP(k,j,i){

      // Convert to Cartesian centered on the star.
      x = x1[i]*sin(x2[j])*cos(x3[k]);
      y = x1[i]*sin(x2[j])*sin(x3[k]);
      z = x1[i]*cos(x2[j]);

      // Convert to Cartesian centered on the planet.
      xp = x - a;
      yp = y;
      zp = z;

      // Convert to spherical centered on the planet. 
      rp2 = EXPAND(xp*xp, + yp*yp, + zp*zp);
      rp = sqrt(rp2);
      thetap = acos(zp/rp);
      phip = atan2(yp, xp);

      /* - Perform calculations in Cartesian coordinates. - */

      // Set the density, pressure, velocity of the stellar wind 
      // and planetary interiar.


      if (rp <= 0.5*Rp) {

        EXPAND(vx1[k][j][i] = 0.0;,
               vx2[k][j][i] = 0.0;,
               vx3[k][j][i] = 0.0;)
        prs[k][j][i] = Pp + (2.0/3.0)*CONST_PI*UNIT_G*RHp*RHp*(Rp*Rp-rp*rp);
        rho[k][j][i] = RHp;

#if PHYSICS == MHD                                   
#if BACKGROUND_FIELD == NO
        bpx = 0.0;
        bpy = 0.0;
        bpz = 16.0*B0p;

        // Convert To spherical basis form Cartesian.
        bpx1 = bpx*sin(x2[j])*cos(x3[k]) + bpy*sin(x2[j])*sin(x3[k]) + bpz*cos(x2[j]);
        bpx2 = bpx*cos(x2[j])*cos(x3[k]) + bpy*cos(x2[j])*sin(x3[k]) - bpz*sin(x2[j]);
        bpx3 = -bpx*sin(x3[k])           + bpy*cos(x3[k]);

        EXPAND(bx1[k][j][i] = bpx1;,            
               bx2[k][j][i] = bpx2;,                 
               bx3[k][j][i] = bpx3;) 
#endif
#if BACKGROUND_FIELD == YES
        EXPAND(bx1[k][j][i] = 0.0;,            
               bx2[k][j][i] = 0.0;,                 
               bx3[k][j][i] = 0.0;) 
#endif
#endif

        d->flag[k][j][i]   |= FLAG_INTERNAL_BOUNDARY;

      } else if (rp > 0.5*Rp && rp <= Rp) {

        EXPAND(vx1[k][j][i] = 0.0;,
               vx2[k][j][i] = 0.0;,
               vx3[k][j][i] = 0.0;)
        prs[k][j][i] = Pp + (2.0/3.0)*CONST_PI*UNIT_G*RHp*RHp*(Rp*Rp-rp*rp);
        rho[k][j][i] = RHp;

#if PHYSICS == MHD                                   
#if BACKGROUND_FIELD == NO
        bpx = 3.0*xp*zp*B0p*pow(Rp, 3)*pow(rp, -5);
        bpy = 3.0*zp*yp*B0p*pow(Rp, 3)*pow(rp, -5);
        bpz = (3.0*zp*zp - rp*rp)*B0p*pow(Rp, 3)*pow(rp, -5);

        // Convert To spherical basis form Cartesian.
        bpx1 = bpx*sin(x2[j])*cos(x3[k]) + bpy*sin(x2[j])*sin(x3[k]) + bpz*cos(x2[j]);
        bpx2 = bpx*cos(x2[j])*cos(x3[k]) + bpy*cos(x2[j])*sin(x3[k]) - bpz*sin(x2[j]);
        bpx3 = -bpx*sin(x3[k])           + bpy*cos(x3[k]);

        EXPAND(bx1[k][j][i] = bpx1;,            
               bx2[k][j][i] = bpx2;,                 
               bx3[k][j][i] = bpx3;) 
#endif
#if BACKGROUND_FIELD == YES
  EXPAND(bx1[k][j][i] = 0.0;,            
         bx2[k][j][i] = 0.0;,                 
         bx3[k][j][i] = 0.0;) 
#endif
#endif
    
        d->flag[k][j][i]   |= FLAG_INTERNAL_BOUNDARY;

      } else if (rp > Rp && rp <= 1.5*Rp) {

        ParkerVelocity(parker, csp, v_escp, rp, rcp, Rp, RHp, Pp);
        prs[k][j][i] = parker[1];
        rho[k][j][i] = parker[2];

        // Cartesian velocity components.
        vx = sin(thetap)*(parker[0]*cos(phip) + 
             sin(x3[k])*x1[i]*omega_fr - 
             sin(phip)*rp*omegap);

        vy = sin(thetap)*(parker[0]*sin(phip) - 
             cos(x3[k])*x1[i]*omega_fr + 
             a*omega_orb + 
             cos(phip)*rp*omegap);

        vz = parker[0]*cos(thetap);

        EXPAND(vx1[k][j][i] = vx*sin(x2[j])*cos(x3[k]) + vy*sin(x2[j])*sin(x3[k]) + vz*cos(x2[j]);,
               vx2[k][j][i] = vx*cos(x2[j])*cos(x3[k]) + vy*cos(x2[j])*sin(x3[k]) - vz*sin(x2[j]);,
               vx3[k][j][i] =-vx*sin(x3[k])            + vy*cos(x3[k]);)

      }

      if (rp > Rp && rp <= 1.5*Rp){
        d->Vc[TRC][k][j][i] = 1.0;
      }

    }
  }

}
/*================================================================================*/

/*================================================================================*/
#if BODY_FORCE != NO
void BodyForceVector(double *v, double *g, double x1, double x2, double x3)
{

  double x, y, z;
  double xp, yp, zp;
  double gx1, gx2, gx3;
  double gpx, gpy, gpz;
  double gpx1, gpx2, gpx3; 

  double Ms, gs, gs_in;
  double Mp, gp, gp_in;
  double RHs, RHp, Rp, Rs;
  double a, omega_orb, omega_fr, rp, rp2;

  Ms = g_inputParam[star_mass]*CONST_Msun/UNIT_MASS;
  RHs = g_inputParam[star_surface_rho]/UNIT_DENSITY;
  Rs = g_inputParam[star_radius]*CONST_Rsun/UNIT_LENGTH;
  Mp = g_inputParam[planet_mass]*0.0009543*CONST_Msun/UNIT_MASS;
  RHp = g_inputParam[planet_surface_rho]/UNIT_DENSITY;
  Rp = g_inputParam[planet_radius]*0.10045*CONST_Rsun/UNIT_LENGTH;
  a = g_inputParam[seperation]*214.9394693836;

  // Rotational frequency (orbit and frame).
  omega_orb = sqrt(UNIT_G*Ms/pow(a, 3));
  omega_fr = omega_orb;

  // Convert to Cartesian.
  x = x1*sin(x2)*cos(x3);
  y = x1*sin(x2)*sin(x3);
  z = x1*cos(x2);
  xp = x - a;
  yp = y;
  zp = z;
  rp2 = EXPAND(xp*xp, + yp*yp, + zp*zp);
  rp = sqrt(rp2);

  // Gravity outside bodies.
  gs = -UNIT_G*Ms/x1/x1;
  gp = -UNIT_G*Mp/rp/rp;

  // Gravity inside bodies.
  gs_in = -(4.0/3.0)*CONST_PI*UNIT_G*RHs;     
  gp_in = -(4.0/3.0)*CONST_PI*UNIT_G*RHp;

  if (x1 > Rs && rp > Rp){ // External gravity + centrifugal + coriolis.
    gx1 = gs*x/x1 + gp*xp/rp;
    gx2 = gs*y/x1 + gp*yp/rp;
    gx3 = gs*z/x1 + gp*zp/rp;
  } else if (x1 < Rs) { // Star interal gravity.
    gx1 = gs_in*x;
    gx2 = gs_in*y;
    gx3 = gs_in*z;
  } else if (rp < Rp) { // Planet interal gravity.
    gx1 = gp_in*xp;
    gx2 = gp_in*yp;
    gx3 = gp_in*zp;
  }

  g[IDIR] = gx1*sin(x2)*cos(x3) + gx2*sin(x2)*sin(x3) + gx3*cos(x2);
  g[JDIR] = gx1*cos(x2)*cos(x3) + gx2*cos(x2)*sin(x3) - gx3*sin(x2);
  g[KDIR] = -gx1*sin(x3)        + gx2*cos(x3);

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

  //printf("cs=%e, v_esc=%e, r=%e, rc=%e, R=%e, RH=%e, P=%e \n", cs, v_esc, r, rc, R, RH, P);

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

  //printf("vel=%e, prs=%e, rho=%e \n", parker[0], parker[1], parker[2]);

}

/*================================================================================*/

