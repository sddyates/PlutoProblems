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

void Init (double *v, double x1, double x2, double x3)
{
    int i, max_iter;
    double r, css, csp, rhos, rhop, pres, prep, Mp, omegas, omegap, omega_fr;
    double Rs, Rp, rotational_frequency, a, B0p, v_escs, v_escp, ga;
    double rp2, rp, rs, rs2, thetas, phis, thetap, phip, Ms, B0s, vs, M_dots, M_dotp, Ts, Tp;
    double lambas, etas, lambap, etap, b, psi, tolerance, v_ws, v_wp, psi_old;
    double v_ws_x1, v_ws_x2, v_ws_x3, v_wp_x1, v_wp_x2, v_wp_x3, omega_orb;
    double Periods, Periodp, PRSs, PRSp, Ps, Pp, RHOs, RHOp, RHs, RHp;
    double Fcentr_s_x1, Fcentr_s_x2, Fcor_s_x1, Fcor_s_x2, Fin_s_x1, Fin_s_x2;
    double Fcentr_p_x1, Fcentr_p_x2, Fcor_p_x1, Fcor_p_x2, Fin_p_x1, Fin_p_x2;
    double function, dfunction, step, rcs, rcp;

    g_smallPressure = 1.0e-5; /**< Small value for pressure fix. */

    /* - quantity | value | units - */
    /* - adiabatic index - */
    g_gamma   = 1.05;
    /* - general quatities - */
    a         = g_inputParam[seperation]*1.49597892e+11/UNIT_LENGTH;
    max_iter  = 100000;
    /* - Planet & Star properties - */                                                    
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
    v_escp    = sqrt(2.0*UNIT_G*Mp/Rp); /* - Planet escape speed - */
    v_escs    = sqrt(2.0*UNIT_G*Ms/Rs); /* - Star escape speed - */
    csp       = sqrt((2.0*Tp)/KELVIN); /* - Stellar sound speed at base - */
    css       = sqrt((2.0*Ts)/KELVIN); /* - Planet sound speed at base - */
    omega_orb = sqrt(UNIT_G*Ms/pow(a,3));
    omegas    = 0.0;//omega_orb;      
    omegap    = omega_orb;
    omega_fr  = 0.0;//omega_orb;
    #if ROTATING_FRAME == YES
        g_OmegaZ = omega_orb;
    #endif
    Pp        = csp*csp*RHp/g_gamma;                                                           
    Ps        = css*css*RHs/g_gamma;
    /* - find radius from star and planet - */
    rs2       = EXPAND(x1*x1,+x2*x2,+x3*x3);
    rs        = sqrt(rs2);
    rp2       = EXPAND((x1-a)*(x1-a),+x2*x2,+x3*x3);
    rp        = sqrt(rp2);
    /* - coordinate transforms - */
    thetap    = acos(x3/rp);
    phip      = atan2(x2,(x1-a));
    thetas    = acos(x3/rs);
    phis      = atan2(x2,x1);
    /* - The critical point - */
    rcs = UNIT_G*Ms/(2.0*css*css);
    rcp = UNIT_G*Mp/(2.0*csp*csp);

    /*****************************************************************/
    /* - In order to obtain the velocity as a function of radial 
         distance from either the star or planet it is nessacery 
         to solve parkers equation and find its roots. There are 
         two roots one for the subsonic and one for the supersonic 
         regions. the transition between the two rigines is given 
         by the critical point rc (see above expressions). The 
         roots are found using the Newton raphson tschnique.        - */
    /******************************************************************/

    /* - Set the density, pressure, velocity of the stellar wind 
         and planetary interiar. - */

    if (rp <= Rp) {
        EXPAND(v[VX1] = 0.0;,
               v[VX2] = 0.0;,
               v[VX3] = 0.0;)
        v[PRS] = Pp + (2.0/3.0)*CONST_PI*UNIT_G*RHp*RHp*(Rp*Rp-rp*rp);
        v[RHO] = RHp;
    } else if (rp <= 10.0*Rp && rp > Rp) {
        // Newton raphson method planet.
        lambap = 0.5*pow(v_escp/csp,2);
        etap = rp/Rp;
        b      = -3.0 - 4.0*log(lambap/2.0) + 4.0*log(etap)+2.0*(lambap/etap);
        i      = 0;
        tolerance = 1.0e-8;

        if (rp <= rcp) {
            psi = 2.e-8;
        } else {
            psi = 2.5;
        }
        do {
            function  = -psi + log(psi) + b;
            dfunction = -1.0 + (1.0/psi);
            step=function/dfunction;
            psi=psi-step;
            i++;
        } while (fabs(step) > tolerance && i < max_iter);

        v_wp = csp*sqrt(psi);

        if(isnan(v_wp)){
          v_wp = 0.0;
        }

        EXPAND(v_wp_x1 = sin(thetap)*(v_wp*cos(phip)+sin(phis)*rs*omega_fr-sin(phip)*rp*omegap);,
               v_wp_x2 = sin(thetap)*(v_wp*sin(phip)-cos(phis)*rs*omega_fr+a*omega_orb + cos(phip)*rp*omegap);,
               v_wp_x3 = v_wp*cos(thetap);)
        PRSp = Pp*exp(lambap*(Rp/rp-1.0)-0.5*pow(v_wp/csp,2));
        RHOp = (RHp/Pp)*PRSp;
        EXPAND(v[VX1] = v_wp_x1;,
               v[VX2] = v_wp_x2;,
               v[VX3] = v_wp_x3;)
        v[PRS] = PRSp;
        v[RHO] = RHOp;
    }

    /* - Set the density, pressure, velocity of the stellar wind 
         and stellar interiar. - */
    if(rs <= Rs){
        EXPAND(v[VX1] = 0.0;,
               v[VX2] = 0.0;,
               v[VX3] = 0.0;)
        v[PRS] = Ps + (2.0/3.0)*CONST_PI*UNIT_G*RHs*RHs*(Rs*Rs-rs*rs);
        v[RHO] = RHs;
    } else if (rs > Rs && rp > 10.0*Rp){
        // Newton raphson method Star.
        lambas = 0.5*pow(v_escs/css,2);
        etas   = rs/Rs;
        b      = -3.0 - 4.0*log(lambas/2.0) + 4.0*log(etas) + 2.0*(lambas/etas);
        i      = 0;
        tolerance = 1.0e-8;

        if (rs <= rcs) {
            psi = 2.e-8;
        } else {
            psi = 2.5;
        }
        do {
            function  = -psi + log(psi) + b;
            dfunction = -1.0 + (1.0/psi);
            step = function/dfunction;
            psi = psi - step;
            i++;
        } while (fabs(step) > tolerance && i < max_iter);

        v_ws = css*sqrt(psi);

        if(isnan(v_ws)){
          v_ws = 0.0;
        }

        EXPAND(v_ws_x1 = sin(thetas)*(v_ws*cos(phis)+sin(phis)*rs*(omega_fr+omegas));,
               v_ws_x2 = sin(thetas)*(v_ws*sin(phis)-cos(phis)*rs*(omega_fr+omegas));,
               v_ws_x3 = v_ws*cos(thetas);)
        PRSs = Ps*exp(lambas*(Rs/rs-1.0)-0.5*pow(v_ws/css,2));
        RHOs = (RHs/Ps)*PRSs;
        EXPAND(v[VX1] = v_ws_x1;,
               v[VX2] = v_ws_x2;,
               v[VX3] = v_ws_x3;)
        v[PRS] = PRSs;
        v[RHO] = RHOs;
    }
  
    /* - Set magnetic fields of the star and planet as the sum of two 
         dipole fields. - */
    #if PHYSICS == MHD
    #if BACKGROUND_FIELD == NO
    if (rs <= 0.5*Rs){
        EXPAND(v[BX1] = 0.0;,
               v[BX2] = 0.0;,
               v[BX3] = 16.0*B0s;)
    } else if (rs > 0.5*Rs && rs <= 1.0*Rs) {
        EXPAND(v[BX1] = 3.0*x1*x3*B0s*pow(Rs, 3)*pow(rs, -5);,
               v[BX2] = 3.0*x3*x2*B0s*pow(Rs, 3)*pow(rs, -5);,
               v[BX3] = (3.0*x3*x3 - rs*rs)*B0s*pow(Rs, 3)*pow(rs, -5);)
    } else if (rp <= 0.5*Rp){
        EXPAND(v[BX1] = 0.0;,
               v[BX2] = 0.0;,
               v[BX3] = 16.0*B0p;) 
    } else if (rp > 0.5*Rp && rp <= 10.0*Rp) {
        EXPAND(v[BX1] = 3.0*(x1 - a)*x3*B0p*pow(Rp, 3)*pow(rp, -5);,
               v[BX2] = 3.0*x3*x2*B0p*pow(Rp, 3)*pow(rp, -5);,
               v[BX3] = (3.0*x3*x3 - rp*rp)*B0p*pow(Rp, 3)*pow(rp, -5);)
    } else if (rs > 1.0*Rs && rp > 10.0*Rp) {
        EXPAND(v[BX1] = 3.0*x1*x3*B0s*pow(Rs, 3)*pow(rs, -5) + 3.0*(x1 - a)*x3*B0p*pow(Rp, 3)*pow(rp, -5);,
               v[BX2] = 3.0*x3*x2*B0s*pow(Rs, 3)*pow(rs, -5) + 3.0*x3*x2*B0p*pow(Rp, 3)*pow(rp, -5);,
               v[BX3] = (3.0*x3*x3 - rs*rs)*B0s*pow(Rs, 3)*pow(rs, -5) + (3.0*x3*x3 - rp*rp)*B0p*pow(Rp, 3)*pow(rp, -5);)
    }
    #endif
    #endif
}
/*================================================================================*/

/*================================================================================*/
void Analysis (const Data *d, Grid *grid){}
/*================================================================================*/

/*================================================================================*/
#if BACKGROUND_FIELD == YES
void BackgroundField (double x1, double x2, double x3, double *B0)                                                    
{                                                                                                                     
    double B0p, B0s, rs2, rs, rp2, rp, a, Rs, Rp;
                                                                                                          
    Rp  = g_inputParam[planet_radius]*0.10045*CONST_Rsun/UNIT_LENGTH;
    Rs  = g_inputParam[star_radius]*CONST_Rsun/UNIT_LENGTH;
    a   = g_inputParam[seperation]*1.49597892e+11/UNIT_LENGTH;
    B0p = g_inputParam[planet_Bfield]/Bc;
    B0s = g_inputParam[star_Bfield]/Bc;
    rs2 = EXPAND(x1*x1,+x2*x2,+x3*x3);
    rs  = sqrt(rs2);
    rp2 = EXPAND((x1-a)*(x1-a),+x2*x2,+x3*x3);
    rp  = sqrt(rp2);

    if (rs <= 0.5*Rs){
        B0[0] = 0.0;  
        B0[1] = 0.0;
        B0[2] = 16.0*B0s;       
    } else if (rs > 0.5*Rs && rs <= 1.0*Rs) {
        B0[0] = 3.0*x1*x3*B0s*pow(Rs, 3)*pow(rs, -5);                                                                                                        
        B0[1] = 3.0*x3*x2*B0s*pow(Rs, 3)*pow(rs, -5);
        B0[2] = (3.0*x3*x3 - rs*rs)*B0s*pow(Rs, 3)*pow(rs, -5);       
    } else if (rp <= 0.5*Rp){
        B0[0] = 0.0;                                                                                                        
        B0[1] = 0.0;
        B0[2] = 16.0*B0p;       
    } else if (rp > 0.5*Rp && rp <= 1.0*Rp) {
        B0[0] = 3.0*(x1 - a)*x3*B0p*pow(Rp, 3)*pow(rp, -5);                                                                                                        
        B0[1] = 3.0*x3*x2*B0p*pow(Rp, 3)*pow(rp, -5);
        B0[2] = (3.0*x3*x3 - rp*rp)*B0p*pow(Rp, 3)*pow(rp, -5);       
    } else if (rs > 1.0*Rs && rp > 1.0*Rp) {
        B0[0] = 3.0*x1*x3*B0s*pow(Rs, 3)*pow(rs, -5) + 3.0*(x1 - a)*x3*B0p*pow(Rp, 3)*pow(rp, -5);
        B0[1] = 3.0*x3*x2*B0s*pow(Rs, 3)*pow(rs, -5) + 3.0*x3*x2*B0p*pow(Rp, 3)*pow(rp, -5);
        B0[2] = ((3.0*x3*x3 - rs*rs)*B0s*pow(Rs, 3)*pow(rs, -5) + (3.0*x3*x3 - rp*rp)*B0p*pow(Rp, 3)*pow(rp, -5));       
    }
}
#endif
/*================================================================================*/


/*================================================================================*/
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid) 
{
    int l, i, j, k, nv, o, max_iter;
    double *x1, *x2, *x3;
    double r, css, csp, rhos, rhop, pres, prep, Mp, omegas, omegap, omega_fr;
    double Rs, Rp, rotational_frequency, a, B0p, v_escs, v_escp, ga;
    double rp2, rp, rs, rs2, thetas, phis, thetap, phip, Ms, B0s, vs, M_dots, M_dotp, Ts, Tp;
    double lambas, etas, lambap, etap, b, psi, tolerance, v_ws, v_wp, psi_old;
    double omega_orb, PRSs, RHOs;
    double Periods, Periodp, PRSp, Ps, Pp, RHOp, RHs, RHp;
    double Fcentr_s_x1, Fcentr_s_x2, Fcor_s_x1, Fcor_s_x2, Fin_s_x1, Fin_s_x2;
    double Fcentr_p_x1, Fcentr_p_x2, Fcor_p_x1, Fcor_p_x2, Fin_p_x1, Fin_p_x2;
    double function, dfunction, step, rcs, rcp; 
    double v_ws_x1, v_ws_x2, v_ws_x3, v_wp_x1, v_wp_x2, v_wp_x3;
 
    /* - quantity | value | units - */
    /* - adiabatic index - */
    g_gamma   = 1.05;
    /* - general quatities - */
    a         = g_inputParam[seperation]*1.49597892e+11/UNIT_LENGTH;
    /* - Planet & Star properties - */
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
    v_escp    = sqrt(2.0*UNIT_G*Mp/Rp); /* - Planet escape speed - */
    v_escs    = sqrt(2.0*UNIT_G*Ms/Rs); /* - Star escape speed - */
    csp       = sqrt((2.0*Tp)/KELVIN); /* - Stellar sound speed at base - */
    css       = sqrt((2.0*Ts)/KELVIN); /* - Planet sound speed at base - */
    omega_orb = sqrt(UNIT_G*Ms/pow(a,3));
    omegas    = 0.0;//omega_orb;
    omegap    = omega_orb;
    omega_fr  = 0.0;//omega_orb;
    Pp        = csp*csp*RHp/g_gamma;
    Ps        = css*css*RHs/g_gamma;
    /* - The critical point - */
    rcs = UNIT_G*Ms/(2.0*css*css);
    rcp = UNIT_G*Mp/(2.0*csp*csp);
    /* - Grid points - */
    x1 = grid[IDIR].x;
    x2 = grid[JDIR].x;
    x3 = grid[KDIR].x;

    if (side == 0) {    /* - Hold stellar and planetary wind bases constant - */

        TOT_LOOP(k,j,i){	

            if (d->Vc[PRS][k][j][i] < g_smallPressure){
                d->Vc[PRS][k][j][i] = g_smallPressure;
            }

            /* - Radial quantities - */
            rs2 = EXPAND(x1[i]*x1[i],+x2[j]*x2[j],+x3[k]*x3[k]);
            rs  = sqrt(rs2);
            thetas = acos(x3[k]/rs);
            phis   = atan2(x2[j],x1[i]);


            // Newton Raphson method for finding E and solvinf for the new x and y 
            // coodinates of star 1.
            double E, xp, yp, e;
            e = 0.0;
            o = 0;
            tolerance = 1.0e-8;
            E = 1.0e-3;
            do {
                function  = E - g_time*omega_orb - e*sin(E);
                dfunction = 1.0 - e*cos(E);
                step = function / dfunction;
                E = E - step;
                o++;
            } while (fabs(step) > tolerance && o < max_iter);
            xp = a*(cos(E)-e);
            yp = a*sin(E);

            rp2 = EXPAND((x1[i]-xp)*(x1[i]-xp), + (x2[j]-yp)*(x2[j]-yp), + x3[k]*x3[k]);
            rp  = sqrt(rp2);
            // - coordinate transforms
            thetap = acos(x3[k]/rp);
            phip   = atan2(x2[j]-yp, x1[i]-xp);

            if (isnan(thetap)){
              /*printf("xp=%f, yp=%f, rp=%f, thetap=%f, phip=%f \n", xp, yp, rp, thetap, phip);*/
              thetap = CONST_PI/2.0;
            }

/*
            // - coordinate transforms
            rp2 = EXPAND((x1[i]-a)*(x1[i]-a), + x2[j]*x2[j], + x3[k]*x3[k]);
            rp  = sqrt(rp2);
            thetap = acos(x3[k]/rp);
            phip   = atan2(x2[j], x1[i]-a);
            thetas = acos(x3[k]/rs);
            phis   = atan2(x2[j],x1[i]);
*/
            // - Set the density, pressure, velocity of the stellar wind and planetary interiar.
            if (rp < 0.5*Rp){

                EXPAND(d->Vc[VX1][k][j][i] = -a*omega_orb*sin(omega_orb*g_time);,
                       d->Vc[VX2][k][j][i] = a*omega_orb*cos(omega_orb*g_time);,
                       d->Vc[VX3][k][j][i] = 0.0;) 
                d->Vc[PRS][k][j][i] = Pp + (2.0/3.0)*CONST_PI*UNIT_G*RHp*RHp*(Rp*Rp-rp*rp);
                d->Vc[RHO][k][j][i] = RHp;
                #if PYHSICS == MHD
                EXPAND(d->Vc[BX1][k][j][i] = 0.0;,
                       d->Vc[BX2][k][j][i] = 0.0;,
                       d->Vc[BX3][k][j][i] = 16.0*B0p;) 
                #endif

                //d->flag[k][j][i] |= FLAG_INTERNAL_BOUNDARY;

            } else if (rp >= 0.5*Rp && rp < Rp) {

                EXPAND(d->Vc[VX1][k][j][i] = -a*omega_orb*sin(omega_orb*g_time);,
                       d->Vc[VX2][k][j][i] = a*omega_orb*cos(omega_orb*g_time);,
                       d->Vc[VX3][k][j][i] = 0.0;) 
                d->Vc[PRS][k][j][i] = Pp + (2.0/3.0)*CONST_PI*UNIT_G*RHp*RHp*(Rp*Rp-rp*rp);
                d->Vc[RHO][k][j][i] = RHp;
                #if PYHSICS == MHD
                EXPAND(d->Vc[BX1][k][j][i] = 3.0*(x1[i]-xp)*x3*B0p*pow(Rp, 3)*pow(rp, -5);,
                       d->Vc[BX2][k][j][i] = 3.0*x3[k]*(x2[j]-yp)*B0p*pow(Rp, 3)*pow(rp, -5);,
                       d->Vc[BX3][k][j][i] = (3.0*x3*x3 - rp*rp)*B0p*pow(Rp, 3)*pow(rp, -5);) 
                #endif

                //d->flag[k][j][i] |= FLAG_INTERNAL_BOUNDARY;

            } else if (rp <= 1.5*Rp && rp > Rp){

                // Newton raphson method planet 
                lambap = 0.5*pow(v_escp/csp,2);
                etap   = rp/Rp;
                b      = -3.0 - 4.0*log(lambap/2.0) + 4.0*log(etap) + 2.0*(lambap/etap);
                o      = 0;
                tolerance = 1.0e-8;

                if (rp <= rcp) {
                    psi = 2.e-8;
                } else {
                    psi = 2.5;
                }
                do {
                    function  = -psi + log(psi) + b;
                    dfunction = -1.0 + (1.0/psi);
                    step=function/dfunction;
                    psi=psi-step;
                    o++;
                } while (fabs(step) > tolerance && o < max_iter);

                v_wp = csp*sqrt(psi);

                if(isnan(v_wp)){
                  v_wp = 0.0;
                }

                EXPAND(v_wp_x1 = sin(thetap)*(v_wp*cos(phip)+sin(phis)*rs*omega_fr - a*omega_orb*sin(omega_orb*g_time) - sin(phip)*rp*omegap);,
                       v_wp_x2 = sin(thetap)*(v_wp*sin(phip)-cos(phis)*rs*omega_fr + a*omega_orb*cos(omega_orb*g_time) + cos(phip)*rp*omegap);,
                       v_wp_x3 = v_wp*cos(thetap);)
                PRSp = Pp*exp(lambap*(Rp/rp-1.0)-0.5*pow(v_wp/csp,2));
                RHOp = (RHp/Pp)*PRSp;
                EXPAND(d->Vc[VX1][k][j][i] = v_wp_x1;,
                       d->Vc[VX2][k][j][i] = v_wp_x2;,
                       d->Vc[VX3][k][j][i] = v_wp_x3;) 
                d->Vc[PRS][k][j][i] = PRSp;
                d->Vc[RHO][k][j][i] = RHOp;
            }

            /* - Set the density, pressure, velocity of the stellar wind 
                 and stellar interiar. - */
            if(rs <= Rs){
                d->flag[k][j][i] |= FLAG_INTERNAL_BOUNDARY;
            } else if (rs<1.5*Rs){

                /* - Newton raphson method Star - */
                lambas = 0.5*pow(v_escs/css,2);
                etas   = rs/Rs;
                b      = -3.0 - 4.0*log(lambas/2.0) + 4.0*log(etas) + 2.0*(lambas/etas);
                o      = 0;
                tolerance = 1.0e-8;

                if (rs <= rcs) {
                    psi = 2.0e-8;
                } else {
                    psi = 2.5;
                }
                do {
                    function  = -psi + log(psi) + b;
                    dfunction = -1.0 + (1.0/psi);
                    step=function/dfunction;
                    psi=psi-step;
                    o++;
                } while (fabs(step) > tolerance && o < max_iter);

                v_ws = css*sqrt(psi);

                if(isnan(v_wp)){
                  v_wp = 0.0;
                }

                EXPAND(v_ws_x1 = sin(thetas)*(v_ws*cos(phis)+sin(phis)*rs*(omega_fr+omegas));,
                       v_ws_x2 = sin(thetas)*(v_ws*sin(phis)-cos(phis)*rs*(omega_fr+omegas));,
                       v_ws_x3 = v_ws*cos(thetas);)
                PRSs = Ps*exp(lambas*(Rs/rs-1.)-0.5*pow(v_ws/css,2));
                RHOs = (RHs/Ps)*PRSs;
                EXPAND(d->Vc[VX1][k][j][i] = v_ws_x1;,
                       d->Vc[VX2][k][j][i] = v_ws_x2;,
                       d->Vc[VX3][k][j][i] = v_ws_x3;)
                d->Vc[PRS][k][j][i] = PRSs;
                d->Vc[RHO][k][j][i] = RHOs;
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
    /* - Rotational frequency (orbit and frame)) - */
    omega_orb = sqrt(UNIT_G*Ms/pow(a, 3));
    omega_fr = 0.0;//omega_orb;
    /* - Distance from star - */
    rs2 = EXPAND(x1*x1, + x2*x2, + x3*x3);
    rs = sqrt(rs2);
    /* - Distance from planet - */
    int o;
    double tolerance, function, dfunction, E, xp, yp, e, max_iter=1e10, step;
    e = 0.0;
    o = 0;
    tolerance = 1.0e-8;
    E = 1.0e-3;
    do {
      function  = E - g_time*omega_orb - e*sin(E);
      dfunction = 1.0 - e*cos(E);
      step = function / dfunction;
      E = E - step;
      o++;
    } while (fabs(step) > tolerance && o < max_iter);
    xp = a*(cos(E)-e);
    yp = a*sin(E);
    rp2 = EXPAND((x1-xp)*(x1-xp), + (x2-yp)*(x2-yp), + x3*x3);
    rp  = sqrt(rp2);
/*
    rp2 = EXPAND((x1-a)*(x1-a), + x2*x2, + x3*x3);
    rp = sqrt(rp2);
*/
    /* - Gravity outside bodies - */
    gs = -UNIT_G*Ms/rs/rs;
    gp = -UNIT_G*Mp/rp/rp;
    /* - Gravity inside bodies - */
    gs_in = -(4.0/3.0)*CONST_PI*UNIT_G*RHs;     
    gp_in = -(4.0/3.0)*CONST_PI*UNIT_G*RHp;
    /* - Coriolis and centrifugal forces - */
    Fin_x1 = omega_fr*omega_fr*x1 + 2.0*omega_fr*v[VX2];
    Fin_x2 = omega_fr*omega_fr*x2 - 2.0*omega_fr*v[VX1];

    if (rs > Rs && rp > Rp){ /* - External gravity + centrifugal + coriolis - */
        g[IDIR] = gs*x1/rs + gp*(x1-xp)/rp + Fin_x1;
        g[JDIR] = gs*x2/rs + gp*(x2-yp)/rp + Fin_x2;
        g[KDIR] = gs*x3/rs + gp*x3/rp;
    } else if (rs < Rs) { /* - Star internal gravity - */
        g[IDIR] = gs_in*x1 + Fin_x1;
        g[JDIR] = gs_in*x2 + Fin_x1;
        g[KDIR] = gs_in*x3;
    } else if (rp < Rp) { /* - Planet internal gravity - */
        g[IDIR] = gp_in*(x1-xp);
        g[JDIR] = gp_in*(x2-yp);
        g[KDIR] = gp_in*x3;
    }
}
#endif
/*================================================================================*/

