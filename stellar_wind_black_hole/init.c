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
/*================================================================================*/
#include "pluto.h"                                                                  
/*================================================================================*/
/*
          Initialise the Grid acouding to C.A.K. steady state Hydro model.          
*/                                                                                  
void Init (double *v, double x1, double x2, double x3){
/*================================================================================*/
 double tyear     = g_inputParam[T_YEAR],                                           
        Mratio    = g_inputParam[M_RATIO],                                          
        L_sun     = g_inputParam[L_SUN],                                            
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
        L         = (Lratio*L_sun*UNIT_L),                                          
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
				r         = 0.0,
				r2        = 0.0,
				rp        = 0.0,
				rb        = 0.0,
				rb2       = 0.0,
				theta     = 0.0;
/*================================================================================*/
	#if EOS == IDEAL
		g_gamma = 1.05;
	#endif
	#if EOS == ISOTHERMAL                                                  
		g_isoSoundSpeed = sqrt(UNIT_kB*T/(mu*(CONST_AH/UNIT_MASS)*CONST_amu));
	#endif
	
	r2 = EXPAND(x1*x1,+x2*x2,+x3*x3);
	r  = sqrt(r2);
	Omega = (0.5*sqrt((8.0*UNIT_G*M_star)/27.0));
	#if EOS == IDEAL                                                              
		v[PRS] = (v[RHO]*T/(KELVIN*mu));
	#endif
	if (r > 1.0) {
		vv = v_inf*pow(1.0 - (1.0/r),b);
		v[RHO] = (M_dot/(4.0*CONST_PI*vv*r2));
		EXPAND(v[VX1] = vv*x1/r;,                                                 
	        v[VX2] = vv*x2/r;,                               
	        v[VX3] = vv*x3/r;)
	} 
	else if (r > 0.9 && r <= 1.0) {
		v[RHO] = (M_dot/(4.0*CONST_PI*(cs/Cs_p)));
		EXPAND(v[VX1] = cs*x1/r;,                                                 
		       v[VX2] = cs*x2/r;,                               
	         v[VX3] = cs*x3/r;)
	}
	else {
		v[RHO] = (M_dot/(4.0*CONST_PI*(cs/Cs_p)));
		EXPAND(v[VX1] = 0.0;,                                                 
	         v[VX2] = 0.0;,                               
		       v[VX3] = 0.0;)
	}        

  rb2 = EXPAND((x1-3.5)*(x1-3.5),+x2*x2,+x3*x3);
  rb  = sqrt(rb2);
  if (rb <= 0.1) {
  	v[RHO] = 1.0e-19/UNIT_DENSITY;
  	EXPAND(v[VX1] = 0.0;,
           v[VX2] = 0.0;,
           v[VX3] = 0.0;)
  }
         

	#if PHYSICS == MHD
  	if (r > 0.5) {
			beta *= 0.0174532925;
			xp = x1*cos(beta) - x2*sin(beta);
			yp = x1*sin(beta) + x2*cos(beta);
			zp = x3;
			rp = sqrt(xp*xp + yp*yp + zp*zp);
			#if DIMENSIONS == 2
				EXPAND(v[BX1] = 3.0*xp*yp*Bq*pow(rp,-5);,            
		  		     v[BX2] = (3.0*pow(yp,2)-pow(rp,2))*Bq*pow(rp,-5);,                 
		    		   v[BX3] = 0.0;)
			#endif
 			#if DIMENSIONS == 3
 				EXPAND(v[BX1] = 3.0*xp*zp*Bq*pow(rp,-5);,            
 				       v[BX2] = 3.0*yp*zp*Bq*pow(rp,-5);,                 
	 			       v[BX3] = (3.0*pow(zp,2)-pow(rp,2))*Bq*pow(rp,-5);)
	 		#endif 
	  	} else {
			#if DIMENSIONS == 2
      	EXPAND(v[BX1] = 0.0;,            
        	     v[BX2] = 16.0*Bq;,                 
          	   v[BX3] = 0.0;)
     	#endif
		 	#if DIMENSIONS == 3
		  	EXPAND(v[BX1] = 0.0;,            
		  			   v[BX2] = 0.0;,                 
		 		  		 v[BX3] = 16.0*Bq;)
		 	#endif
		}                                                 
	#endif

}                                                                          
/*================================================================================*/
void Analysis (const Data *d, Grid *grid){}                               
/*================================================================================*/
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid) {        
/*================================================================================*/
 int i, j, k, p, o, l, m, n, ghost, phase, kprime, jprime;                                                               

 double tyear   = g_inputParam[T_YEAR],
        Cs_p    = g_inputParam[Cs_P],
        Mratio  = g_inputParam[M_RATIO],
        L_sun   = g_inputParam[L_SUN],
        Lratio  = g_inputParam[L_RATIO],
        T       = g_inputParam[TT],
        mu      = g_inputParam[MU],
        a       = g_inputParam[AA],
        b       = g_inputParam[b_law],
        Q       = g_inputParam[QQ],
        a_eff   = g_inputParam[aa_eff],
        beta    = g_inputParam[BB],
        M_star  = (Mratio*CONST_Msun/UNIT_MASS),
        Edd     = (2.6e-5*(Lratio)*(1.0/Mratio)),
        L       = (Lratio*L_sun*UNIT_L),
        c       = 3.0e+5,
        M_dot   = pow(1.0+a_eff,-(1.0/a_eff)) * a_eff * pow(1.0-a_eff,-1)*
                  (L*pow(c,-2))*pow(((Q*Edd)*pow(1.0-Edd,-1)),pow(a_eff,-1)-1.0),
        ke      = ((4.0*CONST_PI*UNIT_G*M_star*c*Edd)/L),
        Omega2  = 10.0*pow(0.5,2)*(8.0/27.0)*UNIT_G*M_star,
        A       = ((1.0/(1.0-a))*((ke*L*Q)/(4.0*CONST_PI*c))),
        Bcgs    = g_inputParam[B_CGS],                                        
        cs      = sqrt(UNIT_kB*T/(mu*(CONST_AH/UNIT_MASS)*CONST_amu)),       
        Bq      = Bcgs/UNIT_B;                                                    

  double        x = 0.0,                 y = 0.0,                 z = 0.0;
	double       xp = 0.0,                yp = 0.0,                zp = 0.0;
	double      *x1 = grid[IDIR].x,      *x2 = grid[JDIR].x,      *x3 = grid[KDIR].x;
  double     *dx1 = grid[IDIR].dx,    *dx2 = grid[JDIR].dx,    *dx3 = grid[KDIR].dx;
	double   ***vx1 = d->Vc[VX1],     ***vx2 = d->Vc[VX2],     ***vx3 = d->Vc[VX3];
	double      Ix1 = 0.0,               Ix2 = 0.0,               Ix3 = 0.0;
	double  vIx1=0.0,                vIx2=0.0,                vIx3=0.0;

	double N0 = 0.0, N1 = 0.0, N2 = 0.0, N3 = 0.0, N4 = 0.0, N5 = 0.0, N6 = 0.0, N7 = 0.0, Ntot = 0.0;
	double P00 = 0.0, P01 = 0.0, P11 = 0.0, P12 = 0.0, P22 = 0.0; 
	double dvdr = 0.0, dvdr2 = 0.0;
	double vr = 0.0, vr2 = 0.0, r = 0.0, r2 = 0.0, rb2 = 0.0, rb = 0.0, theta = 0.0, phi = 0.0;
	double dr2 = 0.0, dr = 0.0, ddr = 0.0, rp = 0.0, dI2 = 0.0, dI = 0.0, vIr[2];
	double gL = 0.0,gLx1 = 0.0, gLx2 = 0.0, gcx1 = 0.0, gcx2 = 0.0, gg = 0.0,gb = 0.0;
  double nu2_c = 0.0, B = 0.0, sigma = 0.0, f = 0.0;
/*================================================================================*/

 beta *= 0.0174532925;
 #if EOS == IDEAL
  g_gamma = 1.05;
 #endif
 #if EOS == ISOTHERMAL                                                  
	g_isoSoundSpeed = sqrt(UNIT_kB*T/(mu*(CONST_AH/UNIT_MASS)*CONST_amu));
 #endif

 if(side == 0){DOM_LOOP(k,j,i){

   r2 = EXPAND(x1[i]*x1[i],+x2[j]*x2[j],+x3[k]*x3[k]);
   r  = sqrt(r2); 
	 xp = x1[i]*cos(beta) - x2[j]*sin(beta);
	 yp = x1[i]*sin(beta) + x2[j]*cos(beta);
	 zp = x3[k];
   rp = sqrt(xp*xp + yp*yp + zp*zp);

	 rb2 = EXPAND((x1[i]-3.5)*(x1[i]-3.5),+x2[j]*x2[j],+x3[k]*x3[k]);
   rb  = sqrt(rb2);
	 if (rb <= 0.1) {
	 	 d->flag[k][j][i] |= FLAG_INTERNAL_BOUNDARY;
	 }

   /* - Set stellar interior to constant - */
   if (r <= 0.9) { 
    d->flag[k][j][i] |= FLAG_INTERNAL_BOUNDARY;
   }

	 /* - Outer layar of star - */
   if (r > 0.9 && r <= 1.1) {
		d->flag[k][j][i] |= FLAG_INTERNAL_BOUNDARY;
/*
    #if EOS != ISOTHERMAL                                                    
    	d->Vc[RHO][k][j][i] = (M_dot/(4.0*CONST_PI*(cs/Cs_p)));
    	d->Vc[PRS][k][j][i] = ((d->Vc[RHO][k][j][i])*T/(KELVIN*mu));
    #endif
    #if EOS == ISOTHERMAL                                                      
    	d->Vc[RHO][k][j][i] = (M_dot/(4.0*CONST_PI*(g_isoSoundSpeed/Cs_p)));
    #endif
*/
/*
		#if PHYSICS == MHD
			#if DIMENSIONS == 2
	 			EXPAND(d->Vc[BX1][k][j][i] = 3.0*xp*yp*Bq*pow(rp,-5);,            
 			  	     d->Vc[BX2][k][j][i] = (3.0*pow(yp,2)-pow(rp,2))*Bq*pow(rp,-5);,                 
			    	   d->Vc[BX3][k][j][i] = 0.0;)
		 	#endif
		 	#if DIMENSIONS == 3
				EXPAND(d->Vc[BX1][k][j][i] = 3.0*xp*yp*Bq*pow(rp,-5);,            
			  	     d->Vc[BX2][k][j][i] = 3.0*yp*zp*Bq*pow(rp,-5);,                 
			    	   d->Vc[BX3][k][j][i] = (3.0*pow(zp,2)-pow(rp,2))*Bq*pow(rp,-5);)
			#endif
		#endif
*/
   }
	 /* - Radiative diriving calculation. - */
	  if (r>1.1 && rb > 0.2) {
			/* - Convert to spherical coordinates - */
			theta = atan2(x1[i],x2[j]);
			/* - Determine infinitesimal ratial distance dr - */
			Ntot = dx1[i]*dx2[j];
			dr2 = EXPAND(dx1[i]*dx1[i],+dx2[j]*dx2[j],+dx3[k]*dx3[k]);
			dr  = sqrt(dr2);
			for (p = 0; p < 2; p++) {
				if (p == 0) {
					ddr = -0.9*dx1[i];
				} else {
					ddr = 0.9*dx1[i];
				}
				/* - Record the coordinates of the interpolation point + - */
					Ix1 = (r+ddr)*sin(theta);
					Ix2 = (r+ddr)*cos(theta);
				/* - Find nearest nighbours and record coordinates - */
				if (Ix1>x1[i] && Ix2>x2[j]){
					N0 = (x1[i+1]-Ix1)*(x2[j+1]-Ix2)/Ntot;
					N1 = (Ix1-x1[i]  )*(x2[j+1]-Ix2)/Ntot;
					N2 = (x1[i+1]-Ix1)*(Ix2-x2[j]  )/Ntot;
					N3 = (Ix1-x1[i]  )*(Ix2-x2[j]  )/Ntot;
					vIx1 = (vx1[k][j][i]*N0+vx1[k][j][i+1]*N1+vx1[k][j+1][i]*N2+vx1[k][j+1][i+1]*N3);
					vIx2 = (vx2[k][j][i]*N0+vx2[k][j][i+1]*N1+vx2[k][j+1][i]*N2+vx2[k][j+1][i+1]*N3);
					vIr[p]  = vIx1*x1[i]/r + vIx2*x2[j]/r;
				} else if (Ix1<x1[i] && Ix2>x2[j]){
					N0 = (x1[i]-Ix1  )*(x2[j+1]-Ix2)/Ntot;
					N1 = (Ix1-x1[i-1])*(x2[j+1]-Ix2)/Ntot;
					N2 = (x1[i]-Ix1  )*(Ix2-x2[j]  )/Ntot;
					N3 = (Ix1-x1[i-1])*(Ix2-x2[j]  )/Ntot;
					vIx1 = (vx1[k][j][i-1]*N0+vx1[k][j][i]*N1+vx1[k][j+1][i-1]*N2+vx1[k][j+1][i]*N3);
					vIx2 = (vx2[k][j][i-1]*N0+vx2[k][j][i]*N1+vx2[k][j+1][i-1]*N2+vx2[k][j+1][i]*N3);
					vIr[p]  = vIx1*x1[i]/r + vIx2*x2[j]/r;
				} else if (Ix1<x1[i] && Ix2<x2[j]){
					N0 = (x1[i]-Ix1  )*(x2[j]-Ix2  )/Ntot;
					N1 = (Ix1-x1[i-1])*(x2[j]-Ix2  )/Ntot;
					N2 = (x1[i]-Ix1  )*(Ix2-x2[j-1])/Ntot;
					N3 = (Ix1-x1[i-1])*(Ix2-x2[j-1])/Ntot;
					vIx1 = (vx1[k][j-1][i-1]*N0+vx1[k][j-1][i]*N1+vx1[k][j][i-1]*N2+vx1[k][j][i]*N3);
					vIx2 = (vx2[k][j-1][i-1]*N0+vx2[k][j-1][i]*N1+vx2[k][j][i-1]*N2+vx2[k][j][i]*N3);
					vIr[p]  = vIx1*x1[i]/r + vIx2*x2[j]/r;
				} else if (Ix1>x1[i] && Ix2<x2[j]){
					N0 = (x1[i+1]-Ix1)*(x2[j]-Ix2  )/Ntot;
					N1 = (Ix1-x1[i]  )*(x2[j]-Ix2  )/Ntot;
					N2 = (x1[i+1]-Ix1)*(Ix2-x2[j-1])/Ntot;
					N3 = (Ix1-x1[i]  )*(Ix2-x2[j-1])/Ntot;
					vIx1 = (vx1[k][j-1][i]*N0+vx1[k][j-1][i+1]*N1+vx1[k][j][i]*N2+vx1[k][j][i+1]*N3);
					vIx2 = (vx2[k][j-1][i]*N0+vx2[k][j-1][i+1]*N1+vx2[k][j][i]*N2+vx2[k][j][i+1]*N3);
					vIr[p]  = vIx1*x1[i]/r + vIx2*x2[j]/r;
				}
			}

			vr = EXPAND(vx1[k][j][i]*x1[i]/r, + vx2[k][j][i]*x2[j]/r, + vx3[k][j][i]*x3[k]/r);
			/* - Use Neville's algorithum and 4th order central finite 
					 diference to calculate the velocity gradient          - */
			P00 = vIr[0];
			P11 = vr;
			P22 = vIr[1];
			P01 = ((0.5*fabs(ddr)*P00)+(0.5*fabs(ddr)*P11))/fabs(ddr);
			P12 = ((0.5*fabs(ddr)*P11)+(0.5*fabs(ddr)*P22))/fabs(ddr);
			dvdr = fabs(((1.0/12.0)*P00)
					   -((2.0/3.0)*P01)
					   +((2.0/3.0)*P12)
					   -((1.0/12.0)*P22))
					   /fabs(0.5*ddr);

			nu2_c = (1.0-(1.0/(r*r)));
			B     = ((d->Vc[RHO][k][j][i])*Q*c*ke);
			sigma = (r/vr)*(dvdr)-1.0;
			f     = ((pow(1.0+sigma,1.0+a)-pow(1.0+sigma*nu2_c,1.0+a))/((1.0+a)*(1.0-nu2_c)*sigma*pow(1.0+sigma,a)));
			gL    = (f*A*pow(r,-2)*pow(dvdr/B,a));
			EXPAND(vx1[k][j][i] += gL*(x1[i]/r)*g_dt;,
						 vx2[k][j][i] += gL*(x2[j]/r)*g_dt;,
						 vx3[k][j][i] += gL*(x3[k]/r)*g_dt;)
		} 
 }}
}                                                                          
/*================================================================================*/

/*================================================================================*/
#if BODY_FORCE != NO
	void BodyForceVector(double *v, double *g, double x1, double x2, double x3)
	{
		double 	Mratio  = g_inputParam[M_RATIO],
						Lratio  = g_inputParam[L_RATIO],
						M_star  = (Mratio*CONST_Msun/UNIT_MASS),
						M_Bhole = (16.0*CONST_Msun/UNIT_MASS),
						Edd     = (2.6e-5*(Lratio)*(1.0/Mratio)),
						sep     = 3.5,
						gg      = 0.0,
						gb      = 0.0,
						r       = 0.0,
						r2      = 0.0,
						rb      = 0.0,
						rb2     = 0.0,
						gcx1    = 0.0,
						gcx2    = 0.0,
						Omega2  = pow(0.5,2)*(8.0/27.0)*UNIT_G*M_star;

		rb2 = EXPAND((x1-sep)*(x1-sep),+x2*x2,+x3*x3);
		rb  = sqrt(rb2);
		gb = -(UNIT_G*M_Bhole)/(rb2);

		r2 = EXPAND(x1*x1, + x2*x2, + x3*x3);
		r  = sqrt(r2);
		gg = -UNIT_G*M_star*(1.0-Edd)/r2;
 
		gcx1 = Omega2*x1 + 2.0*sqrt(Omega2)*v[VX2];
		gcx2 = Omega2*x2 - 2.0*sqrt(Omega2)*v[VX1];
  
		g[IDIR] = gg*x1/r + gb*(x1-sep)/rb + gcx1;
		g[JDIR] = gg*x2/r + gb*x2/rb + gcx2;
		g[KDIR] = gg*x3/r + gb*x3/rb;
	}
#endif
/*================================================================================*/


