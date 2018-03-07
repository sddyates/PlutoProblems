/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

void Init (double *v, double x1, double x2, double x3)
/*
*/
{
    // set the maximum cooling rate and the adiabatic 
    // index (g_gamma = 1 is isothermal).
    g_maxCoolingRate = 0.1;
    g_gamma = 1.05;

    double r1, r2, c = 3.0e+5;
    double displace1, displace2, a1, b1, a2, b2, P1, P2;
    double r_max1, r_min1, r_max2, r_min2;
    double sep, e;

    double Rratio1, Mratio1, Lratio1, T1, mu1, alpha1, Q1, a_eff1; 
    double Rs_1, M_star1, Edd1, L1, M_dot1, cs1, v_esc1, v01;
    double Rratio2, Mratio2, Lratio2, T2, mu2, alpha2, Q2, a_eff2; 
    double Rs_2, M_star2, Edd2, L2, M_dot2, cs2, v_esc2, v02;

    // The mass-loss rates for the stars are set accourding to the equation 
    // from Petit et al. 2013. This method means that the mass-loss rates 
    // are calculated self consitantly from the the base stellar paroperties 
    // and CAK theory.

    // -------------------------------------------------------------------------
    // Get the binary separation and eccentricity.
    sep = g_inputParam[SEP];
    e = g_inputParam[Ep];
    // -------------------------------------------------------------------------

    // -------------------------------------------------------------------------
    // Parameters for star 1.
    Rratio1 = g_inputParam[R1_RATIO], // Radius of star 1.      [Rsun]
    Mratio1 = g_inputParam[M1_RATIO], // Mass of star 1.        [Msun]
    Lratio1 = g_inputParam[L1_RATIO], // Luminosity of star 1.  [Lsun]
    T1      = g_inputParam[TT1s],     // Temperature of star 1. [K]
    mu1     = g_inputParam[MU1s],     // Mean molecular waight of star 1.
    alpha1  = g_inputParam[AA1s],     // CAK alpha parameter star 1.
    Q1      = g_inputParam[QQ1s],     // CAK Q parameter star 1.
    a_eff1  = g_inputParam[aa_eff1s], // CAK effective alpha parameter star 1.
    // Convert to code units.
    Rs_1    = (Rratio1*CONST_Rsun/UNIT_LENGTH),
    M_star1 = (Mratio1*CONST_Msun/UNIT_MASS),                                  
    Edd1    = (2.6e-5*(Lratio1)*(1.0/Mratio1)),                                 
    L1      = (Lratio1*L_sun/UNIT_L),                                          
    // Derive mass-loss rate and sound, escape, and terminal velocities 
    // in code units.
    M_dot1  = pow(1.0+a_eff1,-(1.0/a_eff1)) * a_eff1 * pow(1.0-a_eff1,-1)*
             (L1/(c*c))*pow(((Q1*Edd1)*pow(1.0-Edd1,-1)),pow(a_eff1,-1)-1.0),   
    cs1     = sqrt(UNIT_kB*T1/(mu1*(CONST_AH/UNIT_MASS)*CONST_amu)),            
    v_esc1  = sqrt(2.0*UNIT_G*M_star1*(1.0-Edd1)),                              
    v01     = v_esc1 * sqrt((alpha1/(1.0-alpha1)));   
    // -------------------------------------------------------------------------

    // -------------------------------------------------------------------------
    // Parameters for star 2.
    Rratio2 = g_inputParam[R2_RATIO], // Radius of star 2.      [Rsun]
    Mratio2 = g_inputParam[M2_RATIO], // Mass of star 2.        [Msun]
    Lratio2 = g_inputParam[L2_RATIO], // Luminosity of star 2.  [Lsun]
    T2      = g_inputParam[TT2s],     // Temperature of star 2. [K]
    mu2     = g_inputParam[MU2s],     // Mean molecular waight of star 2.
    alpha2  = g_inputParam[AA2s],     // CAK alpha parameter star 2.
    Q2      = g_inputParam[QQ2s],     // CAK Q parameter star 2.
    a_eff2  = g_inputParam[aa_eff2s], // CAK effective alpha parameter star 2.
    // Convert to code units.
    Rs_2    = (Rratio2*CONST_Rsun/UNIT_LENGTH),
    M_star2 = (Mratio2*CONST_Msun/UNIT_MASS),                                  
    Edd2    = (2.6e-5*(Lratio2)*(1.0/Mratio2)),                                 
    L2      = (Lratio2*L_sun/UNIT_L),                                          
    // Derive mass-loss rate and sound, escape, and terminal velocities
    // in code units.
    M_dot2  = pow(1.0+a_eff2,-(1.0/a_eff2)) * a_eff2 * pow(1.0-a_eff2,-1)*
             (L2/(c*c))*pow(((Q2*Edd2)*pow(1.0-Edd2,-1)),pow(a_eff2,-1)-1.0),   
    cs2     = sqrt(UNIT_kB*T2/(mu2*(CONST_AH/UNIT_MASS)*CONST_amu)),            
    v_esc2  = sqrt(2.0*UNIT_G*M_star2*(1.0-Edd2)),                              
    v02     = v_esc2 * sqrt((alpha2/(1.0-alpha2)));   
    // -------------------------------------------------------------------------

    // Set the minimum temperature for the simulation.
    g_minCoolingTemp = (T1 + T2)/2.0;

    // Set the initial positions of the two stars.
  
    // -------------------------------------------------------------------------
    // Star 1 (right)
   
    // Determin the max and min distance of the star during its orbit.
    r_max1 = sep/(1.0+(M_star1/M_star2));
    r_min1 = r_max1*((1-e)/(1+e));

    displace1 = r_min1; // the star starts off at closses point to other star.

    // Set distance r for star 1.
    r1 = EXPAND((x1-displace1)*(x1-displace1), + x2*x2, + x3*x3);
    r1 = sqrt(r1);

    // Set internal values.
    if (r1 <= Rs_1){
        v[RHO] = M_dot1/(4.0*CONST_PI*Rs_1*Rs_1*v01);
        EXPAND(v[VX1] = 0.0;,
               v[VX2] = 0.0;,
               v[VX3] = 0.0;)
        #if HAVE_ENERGY
        v[PRS] = (v[RHO]*T1/(KELVIN*mu1));
        #endif
    // Set external wind values for the star.
    } else if (r1 > Rs_1 && x1 > 0.0){
        v[RHO] = M_dot1/(4.0*CONST_PI*r1*r1*v01);
        EXPAND(v[VX1] = v01*(x1-displace1)/r1;,
               v[VX2] = v01*x2/r1;,
               v[VX3] = v01*x3/r1;)
        #if HAVE_ENERGY
        v[PRS] = (v[RHO]*T1/(KELVIN*mu1));
        #endif
    }   
    // -------------------------------------------------------------------------

    // -------------------------------------------------------------------------
    // Star 2 (left) 

    // Determin the max and min distance of the star during its orbit.
    r_max2 = sep/(1.0+(M_star2/M_star1));
    r_min2 = r_max2*((1-e)/(1+e));

    displace2 = -r_min2; // the star starts off at closses point to other star.

    // Set distance r for star.
    r2 = EXPAND((x1-displace2)*(x1-displace2), + x2*x2, + x3*x3);
    r2 = sqrt(r2);

    // Set internal values.
    if (r2 <= Rs_2){
        v[RHO] = M_dot2/(4.0*CONST_PI*Rs_2*Rs_2*v02);
        EXPAND(v[VX1] = 0.0;,
               v[VX2] = 0.0;,
               v[VX3] = 0.0;)
        #if HAVE_ENERGY
        v[PRS] = (v[RHO]*T2/(KELVIN*mu2));
        #endif
    // Set external wind values for the star.
    } else if (r2 > Rs_2 && x1 < 0.0){
        v[RHO] = M_dot2/(4.0*CONST_PI*r2*r2*v02);
        EXPAND(v[VX1] = v02*(x1-displace2)/r2;,
               v[VX2] = v02*x2/r2;,
               v[VX3] = v02*x3/r2;)
        #if HAVE_ENERGY
        v[PRS] = (v[RHO]*T2/(KELVIN*mu2));
        #endif
    }
    // -------------------------------------------------------------------------

    #if PHYSICS == MHD || PHYSICS == RMHD
    v[BX1] = 0.0;
    v[BX2] = 0.0;
    v[BX3] = 0.0;
    v[AX1] = 0.0;
    v[AX2] = 0.0;
    v[AX3] = 0.0;
    #endif
    
}

/*================================================================================*/
void Analysis (const Data *d, Grid *grid)
{
}
/*================================================================================*/



#if PHYSICS == MHD
void BackgroundField (double x1, double x2, double x3, double *B0)
/*
*/
{
    B0[0] = 0.0;
    B0[1] = 0.0;
    B0[2] = 0.0;
}
#endif

void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid) 
/* 
*/
{
    int i, j, k, nv;
    double *x1, *x2, *x3;
    double r1, r2, c= 3.0e+5;
    double function, dfunction, step, E, max_iter, tolerance, o;
    double xs1, ys1, xs2, ys2, Mm, t, rm, xsp, ysp;
    double a1, b1, a2, b2, P;
    double r_max1, r_min1, r_max2, r_min2;
    double sep, e;

    double Rratio1, Mratio1, Lratio1, T1, mu1, alpha1, Q1, a_eff1;
    double Rs_1, M_star1, Edd1, L1, M_dot1, cs1, v_esc1, v01;
    double Rratio2, Mratio2, Lratio2, T2, mu2, alpha2, Q2, a_eff2; 
    double Rs_2, M_star2, Edd2, L2, M_dot2, cs2, v_esc2, v02;

    // The mass-loss rate for the star is set accourding to the equation 
    // from Petit et al. 2013. This method means that the mass-loss rate 
    // is 

    // -------------------------------------------------------------------------
    // Get the binary separation and eccentricity.
    sep = g_inputParam[SEP];
    e = g_inputParam[Ep];
    // -------------------------------------------------------------------------

    // -------------------------------------------------------------------------
    // Parameters for star 1.
    Rratio1 = g_inputParam[R1_RATIO], // Radius of star 1.      [Rsun]
    Mratio1 = g_inputParam[M1_RATIO], // Mass of star 1.        [Msun]
    Lratio1 = g_inputParam[L1_RATIO], // Luminosity of star 1.  [Lsun]
    T1      = g_inputParam[TT1s],     // Temperature of star 1. [K]
    mu1     = g_inputParam[MU1s],     // Mean molecular waight of star 1.
    alpha1  = g_inputParam[AA1s],     // CAK alpha parameter star 1.
    Q1      = g_inputParam[QQ1s],     // CAK Q parameter star 1.
    a_eff1  = g_inputParam[aa_eff1s], // CAK effective alpha parameter star 1.
    // Convert to code units.
    Rs_1    = (Rratio1*CONST_Rsun/UNIT_LENGTH),
    M_star1 = (Mratio1*CONST_Msun/UNIT_MASS),                                  
    Edd1    = (2.6e-5*(Lratio1)*(1.0/Mratio1)),                                 
    L1      = (Lratio1*L_sun/UNIT_L),                                          
    // Derive mass-loss rate and sound, escape, and terminal velocities 
    // in code units.
    M_dot1  = pow(1.0+a_eff1,-(1.0/a_eff1)) * a_eff1 * pow(1.0-a_eff1,-1)*
             (L1/(c*c))*pow(((Q1*Edd1)*pow(1.0-Edd1,-1)),pow(a_eff1,-1)-1.0),   
    cs1     = sqrt(UNIT_kB*T1/(mu1*(CONST_AH/UNIT_MASS)*CONST_amu)),            
    v_esc1  = sqrt(2.0*UNIT_G*M_star1*(1.0-Edd1)),                              
    v01     = v_esc1 * sqrt((alpha1/(1.0-alpha1)));   
    // -------------------------------------------------------------------------

    // -------------------------------------------------------------------------
    // Parameters for star 2.
    Rratio2 = g_inputParam[R2_RATIO], // Radius of star 2.      [Rsun]
    Mratio2 = g_inputParam[M2_RATIO], // Mass of star 2.        [Msun]
    Lratio2 = g_inputParam[L2_RATIO], // Luminosity of star 2.  [Lsun]
    T2      = g_inputParam[TT2s],     // Temperature of star 2. [K]
    mu2     = g_inputParam[MU2s],     // Mean molecular waight of star 2.
    alpha2  = g_inputParam[AA2s],     // CAK alpha parameter star 2.
    Q2      = g_inputParam[QQ2s],     // CAK Q parameter star 2.
    a_eff2  = g_inputParam[aa_eff2s], // CAK effective alpha parameter star 2.
    // Convert to code units.
    Rs_2    = (Rratio2*CONST_Rsun/UNIT_LENGTH),
    M_star2 = (Mratio2*CONST_Msun/UNIT_MASS),                                  
    Edd2    = (2.6e-5*(Lratio2)*(1.0/Mratio2)),                                 
    L2      = (Lratio2*L_sun/UNIT_L),                                          
    // Derive mass-loss rate and sound, escape, and terminal velocities
    // in code units.
    M_dot2  = pow(1.0+a_eff2,-(1.0/a_eff2)) * a_eff2 * pow(1.0-a_eff2,-1)*
             (L2/(c*c))*pow(((Q2*Edd2)*pow(1.0-Edd2,-1)),pow(a_eff2,-1)-1.0),   
    cs2     = sqrt(UNIT_kB*T2/(mu2*(CONST_AH/UNIT_MASS)*CONST_amu)),            
    v_esc2  = sqrt(2.0*UNIT_G*M_star2*(1.0-Edd2)),                              
    v02     = v_esc2 * sqrt((alpha2/(1.0-alpha2)));   
    // -------------------------------------------------------------------------

    x1 = grid[IDIR].x;
    x2 = grid[JDIR].x;
    x3 = grid[KDIR].x;

    t = g_time;

    // orbital parameters for star 1.
    r_max1 = sep/(1.0 + (M_star1/M_star2));
    r_min1 = r_max1*((1.0 - e)/(1.0 + e));
    b1 = sqrt(r_max1*r_min1);
    a1 = (r_max1 + r_min1)/2.0;

    // orbital parameters for star 2.
    r_max2 = sep/(1.0+(M_star2/M_star1));
    r_min2 = r_max2*((1.0 - e)/(1.0 + e));
    b2 = sqrt(r_max2*r_min2);
    a2 = (r_max2 + r_min2)/2.0;

    P = 2.0*CONST_PI*sqrt(pow(a1+a2, 3)/(UNIT_G*(M_star2 + M_star1)));
    Mm = 2.0*CONST_PI*t/P;

    // Loop over the internal grid.
    if (side == 0) {
        DOM_LOOP(k,j,i){

            // -------------------------------------------------------------------------
            // Star 1 (right) 

            // Newton Raphson method for finding E and solvinf for the new x and y 
            // coodinates of star 1.
            o = 0;
            tolerance = 1.0e-8;
            E = 1.0e-3;
            do {
                function  = E - Mm - e*sin(E);
                dfunction = 1.0 - e*cos(E);
                step = function / dfunction;
                E = E - step;
                o++;
            } while (fabs(step) > tolerance && o < max_iter);
            xs1 = a1*(cos(E)-e);
            ys1 = b1*sin(E);
            r1 = EXPAND((x1[i] - xs1)*(x1[i] - xs1), + (x2[j] - ys1)*(x2[j] - ys1), + x3[k]*x3[k]);
            r1 = sqrt(r1);

            // Move star to new x and y position.
            if (r1 <= Rs_1){
                d->Vc[RHO][k][j][i] = M_dot1/(4.0*CONST_PI*Rs_1*Rs_1*v01);
                EXPAND(d->Vc[VX1][k][j][i] = 0.0;,
                       d->Vc[VX2][k][j][i] = 0.0;,
                       d->Vc[VX3][k][j][i] = 0.0;)
                #if HAVE_ENERGY
                d->Vc[PRS][k][j][i] = (d->Vc[RHO][k][j][i]*T1/(KELVIN*mu1));
                #endif
            }

            // Hold constant the surface of star 1.
            if (r1 > Rs_1 && r1 < 2.0*Rs_1){
                d->Vc[RHO][k][j][i] = M_dot1/(4.0*CONST_PI*r1*r1*v01);
                d->Vc[PRS][k][j][i] = (d->Vc[RHO][k][j][i]*T1/(KELVIN*mu1)); 
                EXPAND(d->Vc[VX1][k][j][i] = v01*(x1[i]-xs1)/r1;,
                       d->Vc[VX2][k][j][i] = v01*(x2[j]-ys1)/r1;,
                       d->Vc[VX3][k][j][i] = v01*x3[k]/r1;)
            }

            // -------------------------------------------------------------------------
            // Star 2 (left) 

            // Newton Raphson method for finding E and solvinf for the new x and y 
            // coodinates of star 2.
            o = 0;
            tolerance = 1.0e-8;
            E = 1.0e-3;
            do {
                function  = E - Mm - e*sin(E);
                dfunction = 1.0 - e*cos(E);
                step=function/dfunction;
                E = E - step;
                o++;
            } while (fabs(step) > tolerance && o < max_iter);
            xs2 = a2*(cos(E)-e);
            ys2 = b2*sin(E);
            xsp = xs2*cos(CONST_PI) - ys2*sin(CONST_PI);
            ysp = xs2*sin(CONST_PI) + ys2*cos(CONST_PI);
            r2 = EXPAND((x1[i]-xsp)*(x1[i]-xsp), + (x2[j]-ysp)*(x2[j]-ysp), + x3[k]*x3[k]);
            r2 = sqrt(r2);

            // Move star to new x and y position.
            if (r2 <= Rs_2){
                d->Vc[RHO][k][j][i] = M_dot2/(4.0*CONST_PI*Rs_2*Rs_2*v02);
                d->Vc[VX1][k][j][i] = 0.0;
                d->Vc[VX2][k][j][i] = 0.0;
                //d->Vc[VX3][k][j][i] = 0.0;
                #if HAVE_ENERGY
                d->Vc[PRS][k][j][i] = (d->Vc[RHO][k][j][i]*T2/(KELVIN*mu2));
                #endif
            }

            // Hold constant the surface of star 1.
            if (r2 > Rs_2 && r2 < 2.0*Rs_2){
                d->Vc[RHO][k][j][i] = M_dot2/(4.0*CONST_PI*r2*r2*v02); 
                d->Vc[PRS][k][j][i] = (d->Vc[RHO][k][j][i]*T2/(KELVIN*mu1)); 
                EXPAND(d->Vc[VX1][k][j][i] = v02*(x1[i]-xsp)/r2;,
                       d->Vc[VX2][k][j][i] = v02*(x2[j]-ysp)/r2;,
                       d->Vc[VX3][k][j][i] = v02*x3[k]/r2;)
            }
            // -------------------------------------------------------------------------

        }

    }

}

#if CAK == NO
#if BODY_FORCE != NO
/* ********************************************************************* */
void BodyForceVector(double *v, double *g, double x1, double x2, double x3)
/*

*/
{

    // This function sets the graviational vector for both stars and moves 
    // them in the same way as the stars are moved in the function above.

    // All the variables and methods are the same as in the above 
    // function.

    double M_star1, gg1, Edd1, Mratio1, Lratio1;
    double M_star2, gg2, Edd2, Mratio2, Lratio2;
    double r1, r2, xs1, ys1, xs2, ys2;
    double function, dfunction, step, E, max_iter, tolerance, i;
    double Mm, t, e, sep, P, xsp, ysp, xs, ys;
    double a1, b1, r_max1, r_min1;
    double a2, b2, r_max2, r_min2;

    sep = g_inputParam[SEP];
    e = g_inputParam[Ep];

    Mratio1  = g_inputParam[M1_RATIO];
    M_star1  = (Mratio1*CONST_Msun/UNIT_MASS);
    Lratio1  = g_inputParam[L1_RATIO];
    Edd1     = (2.6e-5*(Lratio1)*(1.0/Mratio1));

    Mratio2  = g_inputParam[M2_RATIO];
    M_star2  = (Mratio2*CONST_Msun/UNIT_MASS);
    Lratio2  = g_inputParam[L2_RATIO];
    Edd2     = (2.6e-5*(Lratio2)*(1.0/Mratio2));

    t = g_time;

    // orbital parameters for star 1.
    r_max1 = sep/(1.0 + (M_star1/M_star2));
    r_min1 = r_max1*((1.0 - e)/(1.0 + e));
    b1 = sqrt(r_max1*r_min1);
    a1 = (r_max1 + r_min1)/2.0;
    // orbital parameters for star 2.
    r_max2 = sep/(1.0+(M_star2/M_star1));
    r_min2 = r_max2*((1.0 - e)/(1.0 + e));
    b2 = sqrt(r_max2*r_min2);
    a2 = (r_max2 + r_min2)/2.0;

    P = 2.0*CONST_PI*sqrt(pow(a1+a2, 3)/(UNIT_G*(M_star2 + M_star1)));
    Mm = 2.0*CONST_PI*t/P;

    // Star 1 (right).
    i = 0;
    tolerance = 1.0e-8;
    E = 1.0e-3;
    do {
        function  = E - Mm - e*sin(E);
        dfunction = 1.0 - e*cos(E);
        step=function/dfunction;
        E = E - step;
        i++;
    } while (fabs(step) > tolerance && i < max_iter);
    xs1 = a1*(cos(E)-e);
    ys1 = b1*sin(E);
    r1 = EXPAND((x1 - xs1)*(x1 - xs1), + (x2 - ys1)*(x2 - ys1), + x3*x3);
    r1 = sqrt(r1);

    // Star 2 (right).
    i = 0;
    tolerance = 1.0e-8;
    E = 1.0e-3;
    do {
        function  = E - Mm - e*sin(E);
        dfunction = 1.0 - e*cos(E);
        step = function/dfunction;
        E = E - step;
        i++;
    } while (fabs(step) > tolerance && i < max_iter);
    xs2 = a2*(cos(E)-e);
    ys2 = b2*sin(E);
    xsp = xs2*cos(CONST_PI) - ys2*sin(CONST_PI);
    ysp = xs2*sin(CONST_PI) + ys2*cos(CONST_PI);
    r2 = EXPAND((x1 - xsp)*(x1 - xsp), + (x2 - ysp)*(x2 - ysp), + x3*x3);
    r2 = sqrt(r2);

    gg1 = -UNIT_G*M_star1*(1.0 - Edd1)/r1/r1;
    gg2 = -UNIT_G*M_star2*(1.0 - Edd2)/r2/r2;

    g[IDIR] = gg1*(x1-xs1)/r1 + gg2*(x1-xsp)/r2;
    g[JDIR] = gg1*(x2-ys1)/r1 + gg2*(x2-ysp)/r2;
    g[KDIR] = gg1*x3/r1 + gg2*x3/r2;

}
double BodyForcePotential(double x1, double x2, double x3)
{
  return 0.0;
}
#endif
#endif

