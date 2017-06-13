#include "pluto.h"

/* *************************************************************** */
void ComputeUserVar (const Data *d, Grid *grid)
/*
 *
 *  PURPOSE
 *
 *    Define user-defined output variables
 *
 *
 *
 ***************************************************************** */
{
#if PHYSICS == MHD
#if BACKGROUND_FIELD == YES

  double a, Rp, Rs, B0p, B0s, rs, rs2, rp, rp2;
  double ***Bx, ***By, ***Bz;

  a = g_inputParam[seperation]*1.49597892e+11/UNIT_LENGTH;
  Rp = g_inputParam[planet_radius]*0.10045*CONST_Rsun/UNIT_LENGTH;
  Rs = g_inputParam[star_radius]*CONST_Rsun/UNIT_LENGTH;
  B0p = g_inputParam[planet_Bfield]/Bc;
  B0s = g_inputParam[star_Bfield]/Bc;
  
  DOM_LOOP(k,j,i){

    rs2 = EXPAND(x1[i]*x1[i],+x2[j]*x2[j],+x3[k]*x3[k]);
    rs  = sqrt(rs2);
    rp2 = EXPAND((x1[i]-a)*(x1[i]-a),+x2[j]*x2[j],+x3[k]*x3[k]);
    rp  = sqrt(rp2);

    if (rs <= 0.5*Rs){
      Bx[k][j][i] = 0.0;
      By[k][j][i] = 0.0;
      Bz[k][j][i] = 16.0*B0s;
    } else if (rs > 0.5*Rs && rs <= 1.0*Rs) {

      Bx[k][j][i] = 3.0*x1[i]*x3[k]*B0s*pow(Rs, 3)*pow(rs, -5);
      By[k][j][i] = 3.0*x3[k]*x2[j]*B0s*pow(Rs, 3)*pow(rs, -5);
      Bz[k][j][i] = (3.0*x3[k]*x3[k] - rs*rs)*B0s*pow(Rs, 3)*pow(rs, -5);

    } else if (rp <= 0.5*Rp){

      Bx[k][j][i] = 0.0;
      By[k][j][i] = 0.0;
      Bz[k][j][i] = 16.0*B0p;

    } else if (rp > 0.5*Rp && rp <= sphere*Rp) {

      Bx[k][j][i] = 3.0*(x1[i] - a)*x3*B0p*pow(Rp, 3)*pow(rp, -5);
      By[k][j][i] = 3.0*x3[k]*x2[j]*B0p*pow(Rp, 3)*pow(rp, -5);
      Bz[k][j][i] = (3.0*x3[k]*x3[k] - rp*rp)*B0p*pow(Rp, 3)*pow(rp, -5);

    } else if (rs > 1.0*Rs && rp > sphere*Rp) {

      Bx[k][j][i] = 3.0*x1[i]*x3[k]*B0s*pow(Rs, 3)*pow(rs, -5) + 
          3.0*(x1[i] - a)*x3[k]*B0p*pow(Rp, 3)*pow(rp, -5);
      By[k][j][i] = 3.0*x3[k]*x2[j]*B0s*pow(Rs, 3)*pow(rs, -5) + 
          3.0*x3[k]*x2[j]*B0p*pow(Rp, 3)*pow(rp, -5);
      Bz[k][j][i] = (3.0*x3[k]*x3[k] - rs*rs)*B0s*pow(Rs, 3)*pow(rs, -5) + 
          (3.0*x3[k]*x3[k] - rp*rp)*B0p*pow(Rp, 3)*pow(rp, -5);

    }
  }
#endif
#endif
}
/* ************************************************************* */
void ChangeDumpVar ()
/* 
 *
 * 
 *************************************************************** */
{ 
  Image *image;

}
