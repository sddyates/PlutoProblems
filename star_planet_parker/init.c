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
  double css, csp, rhos, rhop, pres, prep, Mp, omegas, omegap, omega_fr;
  double Rs, Rp, a, B0p, v_escs, v_escp, ga;
  double rp2, rp, rs, rs2, thetas, phis, thetap, phip;
  double Ms, B0s, vs, M_dots, M_dotp, Ts, Tp, omega_orb;
  double Ps, Pp, RHs, RHp, rcs, rcp;
  double kb, mp, sphere;
  double parker[3];
  double br, btheta, bphi, r, theta, phi;
  int n_pole = 2, bx1, bx2, bx3, bp0, bp1, bp2;

  double a11, a12, a13;
  double a21, a22, a23;
  double a31, a32, a33;

  // quantity | value | units

  // adiabatic index

  g_smallPressure = 1.0e-5; // Small value for pressure fix.
  g_gamma = 1.05;

  // general quatities

  a = g_inputParam[seperation]*1.49597892e+11/UNIT_LENGTH;
  kb        = 1.38064852e-16; // [erg K-1]
  mp        = 1.67262171e-24; // [g]

  // Planet & Star properties

  Tp        = g_inputParam[planet_temperature];                                         
  Ts        = g_inputParam[star_temperature];                                           

  Mp        = g_inputParam[planet_mass]*0.0009543*CONST_Msun/M0;                                  
  Ms        = g_inputParam[star_mass]*CONST_Msun/M0;                                    

  Rp        = g_inputParam[planet_radius]*0.10045*CONST_Rsun/UNIT_LENGTH;                       
  Rs        = g_inputParam[star_radius]*CONST_Rsun/UNIT_LENGTH;                         

  B0p       = g_inputParam[planet_Bfield]/Bc;                                           
  B0s       = g_inputParam[star_Bfield]/Bc;                                             

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

  // find radius from star and planet
  rs2       = EXPAND(x1*x1,+x2*x2,+x3*x3);
  rs        = sqrt(rs2);
  rp2       = EXPAND((x1-a)*(x1-a),+x2*x2,+x3*x3);
  rp        = sqrt(rp2);

  // coordinate transforms
  thetap    = acos(x3/rp);
  phip      = atan2(x2,(x1-a));
  thetas    = acos(x3/rs);
  phis      = atan2(x2,x1);

  // The critical point
  rcs = UNIT_G*Ms/(2.0*css*css);
  rcp = UNIT_G*Mp/(2.0*csp*csp);

  sphere = 10.0;

  // Set the density, pressure, velocity of the stellar wind 
  // and planetary interiar.
  if (rp <= Rp) {

    EXPAND(v[VX1] = 0.0;,
           v[VX2] = 0.0;,
           v[VX3] = 0.0;)
    v[PRS] = Pp + (2.0/3.0)*CONST_PI*UNIT_G*RHp*RHp*(Rp*Rp-rp*rp);
    v[RHO] = RHp;

  } else if (rp > Rp && rp <= sphere*Rp) {

    ParkerVelocity(parker, csp, v_escp, rp, rcp, Rp, RHp, Pp);
    EXPAND(v[VX1] = sin(thetap)*(parker[0]*cos(phip)+sin(phis)*rs*omega_fr-sin(phip)*rp*omegap);,
           v[VX2] = sin(thetap)*(parker[0]*sin(phip)-cos(phis)*rs*omega_fr+a*omega_orb + cos(phip)*rp*omegap);,
           v[VX3] = parker[0]*cos(thetap);)
    v[PRS] = parker[1];
    v[RHO] = parker[2];
  }

  if (rp > Rp && rp <= 1.5*Rp){
    v[TRC+1] = 1.0;
  }

  // Set the density, pressure, velocity of the stellar wind 
  // and stellar interiar.
  if(rs <= Rs){

    EXPAND(v[VX1] = 0.0;,
           v[VX2] = 0.0;,
           v[VX3] = 0.0;)
    v[PRS] = Ps + (2.0/3.0)*CONST_PI*UNIT_G*RHs*RHs*(Rs*Rs-rs*rs);
    v[RHO] = RHs;

  } else if (rs > Rs && rp > sphere*Rp){

    ParkerVelocity(parker, css, v_escs, rs, rcs, Rs, RHs, Ps);
    EXPAND(v[VX1] = sin(thetas)*(parker[0]*cos(phis)+sin(phis)*rs*(omega_fr+omegas));,
           v[VX2] = sin(thetas)*(parker[0]*sin(phis)-cos(phis)*rs*(omega_fr+omegas));,
           v[VX3] = parker[0]*cos(thetas);)
    v[PRS] = parker[1];
    v[RHO] = parker[2];
  }
  
  if (rs > Rs && rs <= 1.5*Rs){
    v[TRC] = 1.0;
  }

  // Set magnetic fields of the star and planet as the sum of two 
  // dipole fields.
#if PHYSICS == MHD
#if BACKGROUND_FIELD == NO
  if (n_pole == 2){

    if (rs <= 0.5*Rs){

      EXPAND(v[BX1] = 0.0;,
             v[BX2] = 0.0;,
             v[BX3] = 16.0*B0s;)

    } else if (rs > 0.5*Rs && rs <= Rs) {

      EXPAND(v[BX1] = 3.0*x1*x3*B0s*pow(Rs, 3)*pow(rs, -5);,
             v[BX2] = 3.0*x3*x2*B0s*pow(Rs, 3)*pow(rs, -5);,
             v[BX3] = (3.0*x3*x3 - rs*rs)*B0s*pow(Rs, 3)*pow(rs, -5);)

    } else if (rp <= 0.5*Rp){

      EXPAND(v[BX1] = 0.0;,
             v[BX2] = 0.0;,
             v[BX3] = 16.0*B0p;) 

    } else if (rp > 0.5*Rp && rp <= Rp) {

      EXPAND(v[BX1] = 3.0*(x1 - a)*x3*B0p*pow(Rp, 3)*pow(rp, -5);,
             v[BX2] = 3.0*x3*x2*B0p*pow(Rp, 3)*pow(rp, -5);,
             v[BX3] = (3.0*x3*x3 - rp*rp)*B0p*pow(Rp, 3)*pow(rp, -5);)

    } else if (rs > Rs && rp > Rp) {

      EXPAND(v[BX1] = 3.0*x1*x3*B0s*pow(Rs, 3)*pow(rs, -5) + 
               3.0*(x1 - a)*x3*B0p*pow(Rp, 3)*pow(rp, -5);,
             v[BX2] = 3.0*x3*x2*B0s*pow(Rs, 3)*pow(rs, -5) + 
               3.0*x3*x2*B0p*pow(Rp, 3)*pow(rp, -5);,
             v[BX3] = (3.0*x3*x3 - rs*rs)*B0s*pow(Rs, 3)*pow(rs, -5) + 
               (3.0*x3*x3 - rp*rp)*B0p*pow(Rp, 3)*pow(rp, -5);)
  
    }
  }

  if(n_pole == 4){

    if (rs <= 0.5*Rs){

      EXPAND(v[BX1] = 0.0;,  
             v[BX2] = 0.0;,
             v[BX3] = 16.0*B0s;)

    } else if (rs > 0.5*Rs && rs <= Rs) {

      r = sqrt(x1*x1 + x2*x2 + x3*x3);
      theta = acos(x3/r);
      phi = atan2(x2, x1);

      a11 = sin(theta)*cos(phi); a12 = cos(theta)*cos(phi); a13 = -sin(phi);
      a21 = sin(theta)*sin(phi); a22 = cos(theta)*sin(phi); a23 = cos(phi);
      a31 = cos(theta);          a32 = -sin(theta);         a33 = 0.0;

      br = B0s*pow(Rs/r, 4)*(3.0*pow(cos(theta), 2) - 1.0);
      btheta = 2.0*B0s*pow(Rs/r, 4)*cos(theta)*sin(theta);
      bphi = 0.0;

      bx1 = br*a11 + btheta*a12;
      bx2 = br*a21 + btheta*a22;
      bx3 = br*a31 + btheta*a32;

      EXPAND(v[BX1] = bx1;,
             v[BX2] = bx2;,
             v[BX3] = bx3;)

    } else if (rp <= 0.5*Rp){

      EXPAND(v[BX1] = 0.0;,
             v[BX2] = 0.0;,
             v[BX3] = 16.0*B0p;)       

    } else if (rp > 0.5*Rp && rp <= Rp) {

      EXPAND(v[BX1] = 3.0*(x1 - a)*x3*B0p*pow(Rp, 3)*pow(rp, -5);,
             v[BX2] = 3.0*x3*x2*B0p*pow(Rp, 3)*pow(rp, -5);,
             v[BX3] = (3.0*x3*x3 - rp*rp)*B0p*pow(Rp, 3)*pow(rp, -5);)

    } else if (rs > Rs && rp > Rp) {

      r = sqrt(x1*x1 + x2*x2 + x3*x3);
      theta = acos(x3/r);
      phi = atan2(x2, x1);

      a11 = sin(theta)*cos(phi); a12 = cos(theta)*cos(phi); a13 = -sin(phi);
      a21 = sin(theta)*sin(phi); a22 = cos(theta)*sin(phi); a23 = cos(phi);
      a31 = cos(theta);          a32 = -sin(theta);         a33 = 0.0;

      br = B0s*pow(Rs/r, 4)*(3.0*pow(cos(theta), 2) - 1.0);
      btheta = 2.0*B0s*pow(Rs/r, 4)*cos(theta)*sin(theta);
      bphi = 0.0;

      bx1 = br*a11 + btheta*a12;
      bx2 = br*a21 + btheta*a22;
      bx3 = br*a31 + btheta*a32;

      EXPAND(v[BX1] = bx1 + 3.0*(x1 - a)*x3*B0p*pow(Rp, 3)*pow(rp, -5);,
             v[BX2] = bx2 + 3.0*x3*x2*B0p*pow(Rp, 3)*pow(rp, -5);,
             v[BX3] = bx3 + (3.0*x3*x3 - rp*rp)*B0p*pow(Rp, 3)*pow(rp, -5);)
    }
      
  }

#endif
#endif
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

  double B0p, B0s, rs2, rs, rp2, rp, a, Rs, Rp;
  double br, btheta, bphi, r, theta, phi;
  double bx1, bx2, bx3;
  int n_pole = 2;

  double a11, a12, a13;
  double a21, a22, a23;
  double a31, a32, a33;

  Rp  = g_inputParam[planet_radius]*0.10045*CONST_Rsun/UNIT_LENGTH;
  Rs  = g_inputParam[star_radius]*CONST_Rsun/UNIT_LENGTH;
  a   = g_inputParam[seperation]*1.49597892e+11/UNIT_LENGTH;
  B0p = g_inputParam[planet_Bfield]/Bc;
  B0s = g_inputParam[star_Bfield]/Bc;
  rs2 = EXPAND(x1*x1,+x2*x2,+x3*x3);
  rs  = sqrt(rs2);
  rp2 = EXPAND((x1-a)*(x1-a),+x2*x2,+x3*x3);
  rp  = sqrt(rp2);


  if(n_pole == 2){
 
    if (rs <= 0.5*Rs){

      B0[0] = 0.0;  
      B0[1] = 0.0;
      B0[2] = 16.0*B0s;       

    } else if (rs > 0.5*Rs && rs <= Rs) {

      B0[0] = 3.0*x1*x3*B0s*pow(Rs, 3)*pow(rs, -5);
      B0[1] = 3.0*x3*x2*B0s*pow(Rs, 3)*pow(rs, -5);
      B0[2] = (3.0*x3*x3 - rs*rs)*B0s*pow(Rs, 3)*pow(rs, -5);

    } else if (rp <= 0.5*Rp){

      B0[0] = 0.0;
      B0[1] = 0.0;
      B0[2] = 16.0*B0p;       

    } else if (rp > 0.5*Rp && rp <= Rp) {

      B0[0] = 3.0*(x1 - a)*x3*B0p*pow(Rp, 3)*pow(rp, -5);
      B0[1] = 3.0*x3*x2*B0p*pow(Rp, 3)*pow(rp, -5);
      B0[2] = (3.0*x3*x3 - rp*rp)*B0p*pow(Rp, 3)*pow(rp, -5);       

    } else if (rs > Rs && rp > Rp) {

      B0[0] = 3.0*x1*x3*B0s*pow(Rs, 3)*pow(rs, -5) + 
          3.0*(x1 - a)*x3*B0p*pow(Rp, 3)*pow(rp, -5);
      B0[1] = 3.0*x3*x2*B0s*pow(Rs, 3)*pow(rs, -5) + 
          3.0*x3*x2*B0p*pow(Rp, 3)*pow(rp, -5);
      B0[2] = (3.0*x3*x3 - rs*rs)*B0s*pow(Rs, 3)*pow(rs, -5) + 
          (3.0*x3*x3 - rp*rp)*B0p*pow(Rp, 3)*pow(rp, -5);

    }
  }

  if(n_pole == 4){
    if (rs <= 0.5*Rs){

      B0[0] = 0.0;  
      B0[1] = 0.0;
      B0[2] = 16.0*B0s;       

    } else if (rs > 0.5*Rs && rs <= Rs) {

      r = sqrt(x1*x1 + x2*x2 + x3*x3);
      theta = acos(x3/r);
      phi = atan2(x2, x1);

      a11 = sin(theta)*cos(phi); a12 = cos(theta)*cos(phi); a13 = -sin(phi);
      a21 = sin(theta)*sin(phi); a22 = cos(theta)*sin(phi); a23 = cos(phi);
      a31 = cos(theta);          a32 = -sin(theta);         a33 = 0.0;

      br = B0s*pow(Rs/r, 4)*(3.0*pow(cos(theta), 2) - 1.0);
      btheta = 2.0*B0s*pow(Rs/r, 4)*cos(theta)*sin(theta);
      bphi = 0.0;

      bx1 = br*a11 + btheta*a12;
      bx2 = br*a21 + btheta*a22;
      bx3 = br*a31 + btheta*a32;

      B0[0] = bx1;
      B0[1] = bx2;
      B0[2] = bx3;

    } else if (rp <= 0.5*Rp){

      B0[0] = 0.0;
      B0[1] = 0.0;
      B0[2] = 16.0*B0p;       

    } else if (rp > 0.5*Rp && rp <= Rp) {

      B0[0] = 3.0*(x1 - a)*x3*B0p*pow(Rp, 3)*pow(rp, -5);
      B0[1] = 3.0*x3*x2*B0p*pow(Rp, 3)*pow(rp, -5);
      B0[2] = (3.0*x3*x3 - rp*rp)*B0p*pow(Rp, 3)*pow(rp, -5);       

    } else if (rs > Rs && rp > Rp) {

      r = sqrt(x1*x1 + x2*x2 + x3*x3);
      theta = acos(x3/r);
      phi = atan2(x2, x1);

      a11 = sin(theta)*cos(phi); a12 = cos(theta)*cos(phi); a13 = -sin(phi);
      a21 = sin(theta)*sin(phi); a22 = cos(theta)*sin(phi); a23 = cos(phi);
      a31 = cos(theta);          a32 = -sin(theta);         a33 = 0.0;

      br = B0s*pow(Rs/r, 4)*(3.0*pow(cos(theta), 2) - 1.0);
      btheta = 2.0*B0s*pow(Rs/r, 4)*cos(theta)*sin(theta);
      bphi = 0.0;

      bx1 = br*a11 + btheta*a12;
      bx2 = br*a21 + btheta*a22;
      bx3 = br*a31 + btheta*a32;

      B0[0] = bx1 + 3.0*(x1 - a)*x3*B0p*pow(Rp, 3)*pow(rp, -5);
      B0[1] = bx2 + 3.0*x3*x2*B0p*pow(Rp, 3)*pow(rp, -5);
      B0[2] = bx3 + (3.0*x3*x3 - rp*rp)*B0p*pow(Rp, 3)*pow(rp, -5);       

    }
      
  }

}
#endif
/*================================================================================*/


/*================================================================================*/
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid) 
{
  int i, j, k, o;
  double *x1, *x2, *x3;
  double r, css, csp, Mp, omegas, omegap, omega_fr;
  double Rs, Rp, a, B0p, v_escs, v_escp, ga;
  double rp2, rp, rs, rs2, thetas, phis, thetap, phip;
  double Ms, B0s, vs, M_dots, M_dotp, Ts, Tp;
  double omega_orb;
  double Ps, Pp, RHs, RHp;
  double rcs, rcp; 
  double kb, mp;
  double parker[3];

  double theta, phi;
  double br, btheta, bphi;
  double bx1, bx2, bx3;

  double pole = 2;
  double a11, a12, a13;
  double a21, a22, a23;
  double a31, a32, a33;

  // - quantity | value | units.

  // - adiabatic index.

  g_gamma   = 1.05;

  // general quatities.

  a         = g_inputParam[seperation]*1.49597892e+11/UNIT_LENGTH;
  kb        = 1.38064852e-16; // [erg K-1]
  mp        = 1.67262171e-24; // [g]

  // Planet & Star properties.

  Tp        = g_inputParam[planet_temperature];
  Ts        = g_inputParam[star_temperature];

  Mp        = g_inputParam[planet_mass]*0.0009543*CONST_Msun/M0;
  Ms        = g_inputParam[star_mass]*CONST_Msun/M0;

  Rp        = g_inputParam[planet_radius]*0.10045*CONST_Rsun/UNIT_LENGTH;
  Rs        = g_inputParam[star_radius]*CONST_Rsun/UNIT_LENGTH;

  B0p       = g_inputParam[planet_Bfield]/Bc;
  B0s       = g_inputParam[star_Bfield]/Bc;

  RHp       = g_inputParam[planet_surface_rho]/UNIT_DENSITY;
  RHs       = g_inputParam[star_surface_rho]/UNIT_DENSITY;

  v_escp    = sqrt(2.0*UNIT_G*Mp/Rp);
  v_escs    = sqrt(2.0*UNIT_G*Ms/Rs);

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

  if (side == 0) {
    TOT_LOOP(k,j,i){	

      /* - Radial quantities - */
      rs2 = EXPAND(x1[i]*x1[i],+x2[j]*x2[j],+x3[k]*x3[k]);
      rs  = sqrt(rs2);
      rp2 = EXPAND((x1[i]-a)*(x1[i]-a),+x2[j]*x2[j],+x3[k]*x3[k]);
      rp  = sqrt(rp2);

      /* - coordinate transforms - */
      thetap = acos(x3[k]/rp);
      phip   = atan2(x2[j],(x1[i]-a));
      thetas = acos(x3[k]/rs);
      phis   = atan2(x2[j],x1[i]);

      // Set the density, pressure, velocity of the planetary wind 
      // and planetary interiar.
      if (rp <= 0.5*Rp){

        EXPAND(d->Vc[VX1][k][j][i] = 0.0;,
               d->Vc[VX2][k][j][i] = 0.0;,
               d->Vc[VX3][k][j][i] = 0.0;)
        d->Vc[PRS][k][j][i] = Pp + (2.0/3.0)*CONST_PI*UNIT_G*RHp*RHp*(Rp*Rp-rp*rp);
        d->Vc[RHO][k][j][i] = RHp;
#if PHYSICS == MHD
#if BACKGROUND_FIELD == NO
        EXPAND(d->Vc[BX1][k][j][i] = 0.0;,
               d->Vc[BX2][k][j][i] = 0.0;,
               d->Vc[BX3][k][j][i] = 16.0*B0p;) 
#endif
#if BACKGROUND_FIELD == YES
        EXPAND(d->Vc[BX1][k][j][i] = 0.0;,
               d->Vc[BX2][k][j][i] = 0.0;,
               d->Vc[BX3][k][j][i] = 0.0;) 
#endif
#endif
        d->flag[k][j][i]   |= FLAG_INTERNAL_BOUNDARY;

      } else if (rp > 0.5*Rp && rp <= Rp) {

        EXPAND(d->Vc[VX1][k][j][i] = 0.0;,
               d->Vc[VX2][k][j][i] = 0.0;,
               d->Vc[VX3][k][j][i] = 0.0;)
        d->Vc[PRS][k][j][i] = Pp + (2.0/3.0)*CONST_PI*UNIT_G*RHp*RHp*(Rp*Rp-rp*rp);
        d->Vc[RHO][k][j][i] = RHp;
#if PHYSICS == MHD
#if BACKGROUND_FIELD == NO
        EXPAND(d->Vc[BX1][k][j][i] = 3.0*(x1[i] - a)*x3[k]*B0p*pow(Rp, 3)*pow(rp, -5);,
               d->Vc[BX2][k][j][i] = 3.0*x3[k]*x2[j]*B0p*pow(Rp, 3)*pow(rp, -5);,
               d->Vc[BX3][k][j][i] = (3.0*x3[k]*x3[k] - rp*rp)*B0p*pow(Rp, 3)*pow(rp, -5);)
#endif
#if BACKGROUND_FIELD == YES
        EXPAND(d->Vc[BX1][k][j][i] = 0.0;,
               d->Vc[BX2][k][j][i] = 0.0;,
               d->Vc[BX3][k][j][i] = 0.0;) 
#endif
#endif
        d->flag[k][j][i]   |= FLAG_INTERNAL_BOUNDARY;

      } else if (rp > Rp && rp <= 1.5*Rp){

        ParkerVelocity(parker, csp, v_escp, rp, rcp, Rp, RHp, Pp);
        EXPAND(d->Vc[VX1][k][j][i] = sin(thetap)*(parker[0]*cos(phip)+sin(phis)*
                                       rs*omega_fr-sin(phip)*rp*omegap);,
               d->Vc[VX2][k][j][i] = sin(thetap)*(parker[0]*sin(phip)-cos(phis)*
                                       rs*omega_fr+a*omega_orb + cos(phip)*rp*omegap);,
               d->Vc[VX3][k][j][i] = parker[0]*cos(thetap);) 
        d->Vc[PRS][k][j][i] = parker[1];
        d->Vc[RHO][k][j][i] = parker[2];
        //d->flag[k][j][i]   |= FLAG_INTERNAL_BOUNDARY;
      }

      if (rp > Rp && rp <= 1.5*Rp){
        d->Vc[TRC+1][k][j][i] = 1.0;
      }

      // Set the density, pressure, velocity of the stellar wind 
      // and stellar interiar.
      if (rs <= 0.5*Rs){

        EXPAND(d->Vc[VX1][k][j][i] = 0.0;,
               d->Vc[VX2][k][j][i] = 0.0;,
               d->Vc[VX3][k][j][i] = 0.0;)
        d->Vc[PRS][k][j][i] = Ps + (2.0/3.0)*CONST_PI*UNIT_G*RHs*RHs*(Rs*Rs-rs*rs);
        d->Vc[RHO][k][j][i] = RHs;
#if PHYSICS == MHD
#if BACKGROUND_FIELD == NO
        EXPAND(d->Vc[BX1][k][j][i] = 0.0;,
               d->Vc[BX2][k][j][i] = 0.0;,
               d->Vc[BX3][k][j][i] = 16.0*B0s;)
#endif
#if BACKGROUND_FIELD == YES
        EXPAND(d->Vc[BX1][k][j][i] = 0.0;,
               d->Vc[BX2][k][j][i] = 0.0;,
               d->Vc[BX3][k][j][i] = 0.0;) 
#endif
#endif
        d->flag[k][j][i]   |= FLAG_INTERNAL_BOUNDARY;

      } else if (rs > 0.5*Rs && rs <= Rs) {

        EXPAND(d->Vc[VX1][k][j][i] = 0.0;,
               d->Vc[VX2][k][j][i] = 0.0;,
               d->Vc[VX3][k][j][i] = 0.0;)
        d->Vc[PRS][k][j][i] = Ps + (2.0/3.0)*CONST_PI*UNIT_G*RHs*RHs*(Rs*Rs-rs*rs);
        d->Vc[RHO][k][j][i] = RHs;
#if PHYSICS == MHD
#if BACKGROUND_FIELD == NO

        if (pole == 4){

          r = sqrt(x1[i]*x1[i] + x2[j]*x2[j] + x3[k]*x3[k]);
          theta = acos(x3[k]/r);
          phi = atan2(x2[j], x1[i]);

          a11 = sin(theta)*cos(phi); a12 = cos(theta)*cos(phi); a13 = -sin(phi);
          a21 = sin(theta)*sin(phi); a22 = cos(theta)*sin(phi); a23 = cos(phi);
          a31 = cos(theta);          a32 = -sin(theta);         a33 = 0.0;

          br = B0s*pow(Rs/r, 4)*(3.0*pow(cos(theta), 2) - 1.0);
          btheta = 2.0*B0s*pow(Rs/r, 4)*cos(theta)*sin(theta);
          bphi = 0.0;

          bx1 = br*a11 + btheta*a12;
          bx2 = br*a21 + btheta*a22;
          bx3 = br*a31 + btheta*a32;

          EXPAND(d->Vc[BX1][k][j][i] = bx1;,
                 d->Vc[BX2][k][j][i] = bx2;,
                 d->Vc[BX3][k][j][i] = bx3;)

        }else if (pole == 2){

          EXPAND(d->Vc[BX1][k][j][i] = 3.0*x1[i]*x3[k]*B0s*pow(Rs, 3)*pow(rs, -5);,
                 d->Vc[BX2][k][j][i] = 3.0*x3[k]*x2[j]*B0s*pow(Rs, 3)*pow(rs, -5);,
                 d->Vc[BX3][k][j][i] = (3.0*x3[k]*x3[k] - rs*rs)*B0s*pow(Rs, 3)*pow(rs, -5);)

        }
#endif
#if BACKGROUND_FIELD == YES
        EXPAND(d->Vc[BX1][k][j][i] = 0.0;,
               d->Vc[BX2][k][j][i] = 0.0;,
               d->Vc[BX3][k][j][i] = 0.0;) 
#endif
#endif
        d->flag[k][j][i]   |= FLAG_INTERNAL_BOUNDARY;

      } else if (rs > Rs && rs <= 1.5*Rs){

        ParkerVelocity(parker, css, v_escs, rs, rcs, Rs, RHs, Ps);
        EXPAND(d->Vc[VX1][k][j][i] = sin(thetas)*(parker[0]*cos(phis)+sin(phis)*
                                       rs*(omega_fr+omegas));,
               d->Vc[VX2][k][j][i] = sin(thetas)*(parker[0]*sin(phis)-cos(phis)*
                                       rs*(omega_fr+omegas));,
               d->Vc[VX3][k][j][i] = parker[0]*cos(thetas);)
        d->Vc[PRS][k][j][i] = parker[1];
        d->Vc[RHO][k][j][i] = parker[2];
        d->Vc[TRC][k][j][i] = 1.0;
        //d->flag[k][j][i]   |= FLAG_INTERNAL_BOUNDARY;

      }

      if (rs > Rs && rs <= 1.5*Rs){
        d->Vc[TRC][k][j][i] = 1.0;
      }

      if (d->Vc[PRS][k][j][i] < g_smallPressure){
        d->Vc[PRS][k][j][i] = g_smallPressure;
      }

    }
  }
}
/*================================================================================*/

/*================================================================================*/
#if BODY_FORCE != NO
void BodyForceVector(double *v, double *g, double x1, double x2, double x3)
{

  double Ms, Mp, r, a;
  double rs2, rs, rp2, rp, RHp, RHs, Rs, Rp, omega_fr, omega_orb;
  double Fcentr_s_x1, Fcentr_s_x2, Fcor_s_x1, Fcor_s_x2, Fin_s_x1, Fin_s_x2;
  double Fcentr_p_x1, Fcentr_p_x2, Fcor_p_x1, Fcor_p_x2, Fin_p_x1, Fin_p_x2;
  double Fin_x1, Fin_x2, gs, gp, gs_in, gp_in;

  Ms = g_inputParam[star_mass]*CONST_Msun/M0;
  RHs = g_inputParam[star_surface_rho]/UNIT_DENSITY;
  Rs = g_inputParam[star_radius]*CONST_Rsun/UNIT_LENGTH;
  Mp = g_inputParam[planet_mass]*0.0009543*CONST_Msun/M0;
  RHp = g_inputParam[planet_surface_rho]/UNIT_DENSITY;
  Rp = g_inputParam[planet_radius]*0.10045*CONST_Rsun/UNIT_LENGTH;
  a = g_inputParam[seperation]*1.49597892e+11/UNIT_LENGTH;

  // Rotational frequency (orbit and frame).
  omega_orb = sqrt(UNIT_G*Ms/pow(a, 3));
  omega_fr = omega_orb;

  // Distance from star.
  rs2 = EXPAND(x1*x1, + x2*x2, + x3*x3);
  rs = sqrt(rs2);

  // Distance from planet.
  rp2 = EXPAND((x1-a)*(x1-a), + x2*x2, + x3*x3);
  rp = sqrt(rp2);

  // Gravity outside bodies.
  gs = -UNIT_G*Ms/rs/rs;
  gp = -UNIT_G*Mp/rp/rp;

  // Gravity inside bodies.
  gs_in = -(4.0/3.0)*CONST_PI*UNIT_G*RHs;     
  gp_in = -(4.0/3.0)*CONST_PI*UNIT_G*RHp;

  // Coriolis and centrifugal forces.
  Fin_x1 = omega_fr*omega_fr*x1 + 2.0*omega_fr*v[VX2];
  Fin_x2 = omega_fr*omega_fr*x2 - 2.0*omega_fr*v[VX1];

  if (rs > Rs && rp > Rp){ // External gravity + centrifugal + coriolis.
    g[IDIR] = gs*x1/rs + gp*(x1 - a)/rp + Fin_x1;
    g[JDIR] = gs*x2/rs + gp*x2/rp + Fin_x2;
    g[KDIR] = gs*x3/rs + gp*x3/rp;
  } else if (rs < Rs) { // Star interal gravity.
    g[IDIR] = gs_in*x1;
    g[JDIR] = gs_in*x2;
    g[KDIR] = gs_in*x3;
  } else if (rp < Rp) { // Planet interal gravity.
    g[IDIR] = gp_in*(x1 - a);
    g[JDIR] = gp_in*x2;
    g[KDIR] = gp_in*x3;
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

