#include "main.h"

double mvec(double x, double y)
{
	return(pow((x * x)  + (y * y), 0.5));
}


void check_LSS(double SRRx, double SRRy, double prs, double &tau)
{
// LSS depends on pressure (tauL = lambda * p). 
// Lambda is a sliding friction coefficent, obtained from traction curves
// Habchi, 2018. p. 426
// tau and prs must be in dimension form
	
	double lambda, lss, ntau;

	// we set lambda as constant (0.083), Habchi, 2018. p. 427
	//lambda = 0.083;
	lambda = 0.04;
	lss = prs * lambda;

	if (tau > lss) ntau = lss; else ntau = tau;

	tau = ntau;
	
} // proc


double calc_integrals(dbvector f, double a, double dz, int n, int kval, int q)
	{
		double z_curr, summ, res;
				
		summ = 0.0;

		switch (kval) {

		case 1:
			/* integrate (f) */
			{

				if (q == 0) { /* rectangular rule */
					for (int i = 1; i < n; ++i)
						summ = summ + 0.5 * ((f[i]) + (f[i - 1]));

					res = summ * dz;
				}
				if (q == 1) { /* parabolic trapezoid (3/8) rule */
					summ = (f[0]) + (f[n - 1]);
					for (int i = 1; i < (n - 1); i++) {
						if (i % 3 == 0)
							summ = summ + 2 * (f[i]);
						else
							summ = summ + 3 * (f[i]);
					}

					res = summ * (3.0 / 8.0) * dz;

				}

				return (res);
				break;
			}

		case 2:
			/* integrate (z*f) */
			{
				if (q == 0) { /* rectangular rule */
					for (int i = 1; i < n; ++i) {
						summ = summ + 0.5 * (((a + i * dz) * f[i]) + ((a + (i - 1) * dz) * f[i - 1]));
					}
					res = summ * dz;
				}
				if (q == 1) { /* parabolic trapezoid (3/8) rule */

					/* first z_curr = a; last z_curr = a+(n-1)*dz */
					summ = (a * f[0]) + ((a + (n - 1) * dz) * f[n - 1]);

					for (int i = 1; i < (n - 1); i++) {
						if (i % 3 == 0)
							summ = summ + 2 * ((a + i * dz) * f[i]);
						else
							summ = summ + 3 * ((a + i * dz) * f[i]);
					}

					res = summ * (3.0 / 8.0) * dz;

				}

				return (res);
				break;
			}
		} /* switch */
	} /* proc */


void calc_physical_fields(Domain &D, MM_density MMrho, MM_ls_viscosity MMlsvs, MM_gn_viscosity MMgnvs,
	MM_th_conductivity MMthc, MM_ht_capacity MMhcp, double vel1[2], double vel2[2])
{
	
	int		i, j, l, cnt, cnt_max, Nx, Ny, Nz;
	double	ph, rhoR, muR, etaR, TR, KR, CpR;
	double	a, b, h, dx, dy, dz, dpdx, dpdy;
	double	T_curr, p_curr, rho_curr, z1_curr, z2_curr, d1, d2, z_curr;

	double p1, p2, p3, p4, p5, p6, p7;
	double e1, e2, toler;
	double rtau_zx0, rtau_zy0, taux_curr, tauy_curr, tau_mod, tau0;
	double rtau_mod, rtauL;
	double a11, a12, a22, b1, b2;

	int jend, q;

	/* vectors for storing function values in integrals */
	dbvector intf1, intf2, intf3, intf4, intf5;

	/* vectors for storing tau_zx0, tau_zy0 approx. values in N-R solving procedure */
	dbvector tau_zx0, tau_zy0;

	/* get domain data */
	a	 = D.a;
	b	 = D.b;
	Nx   = D.Nx;
	Ny   = D.Ny;
	dx   = D.dx;
	dy   = D.dy;
	rhoR = D.rhoR;
	muR  = D.muR;
	TR   = D.T0;
	ph	 = D.ph;
	KR	 = D.KR;
	CpR  = D.CpR;
	etaR = D.etaR;
		
	q = 0; /* numerical integration method: 0 - rectangular; 1 - parabolic, (Nz-1) must be a multiple of 3 */
	
	jend = 1; /* line contact is set by default */
	if (D.dkey == 2) jend = Ny;

	for (i = 0; i < Nx; ++i){
		for (j = 0; j < jend; ++j) {
			
			Nz = D.Area[i][j].Nz;
			
			
			d1 = D.Area[i][j].d1;	/* m */
			d2 = D.Area[i][j].d2;	/* m */

			z1_curr = D.Area[i][j].z1 - d1;	/* m */
			z2_curr = D.Area[i][j].z2 + d2;	/* m */

			h = (z2_curr - z1_curr) + D.h0;	/* m */
			dz = h / (Nz - 1);		/* m */
			
			for (l = 0; l < Nz; ++l) {
				
				T_curr = D.Area[i][j].T[l] * TR;
				p_curr = D.Area[i][j].p * ph;
				
				// we store DIMENSION values for all parameters

				/* -- DENSITY -- */
				{
					/* this is a common terms for Tait and Murnaghan eqns */
					p1 = (1 + MMrho.arg4 * (T_curr - MMrho.arg2));
					p2 = 1 / (1 + MMrho.arg5);
					p5 = 1 / MMrho.arg5;

					/* constant */
					if (MMrho.mkey == 0) {

						D.Area[i][j].rho[l] = MMrho.arg1;
					}
					
					/* Dowson & Higginson */
					if (MMrho.mkey == 1) {
						
						D.Area[i][j].rho[l] = MMrho.arg1 * (1 + MMrho.arg2 * p_curr / (1 + MMrho.arg3 * p_curr));
						
					}

					/* Tait */
					if (MMrho.mkey == 2) {

						p3 = 1 + (1 + MMrho.arg5) * p_curr / (MMrho.arg6 * exp(-1.0 * MMrho.arg3 * T_curr));
						if (p_curr == 0) D.Area[i][j].rho[l] = MMrho.arg1;
						else D.Area[i][j].rho[l] = MMrho.arg1 / (p1 * (1 - p2 * log(p3))) ;
					}

					/* Murnaghan */
					if (MMrho.mkey == 3) {

						p4 = 1 + MMrho.arg5 * p_curr / (MMrho.arg6 * exp(-1.0 * MMrho.arg4 * T_curr));
						if (p_curr == 0) D.Area[i][j].rho[l] = MMrho.arg1;
						else D.Area[i][j].rho[l] = (MMrho.arg1 / p1) * pow(p4, p5);
					}

				}

				/* -- LS-VISCOSITY -- */
					{
					/* constant */
					if (MMlsvs.mkey == 0) {
						D.Area[i][j].mu[l] = MMlsvs.arg1;
					}
						
					/* Barus */
					if (MMlsvs.mkey == 1) {

						D.Area[i][j].mu[l] = MMlsvs.arg1 * exp(MMlsvs.arg2 * p_curr);
					}
						
					/* Roeland */
					if (MMlsvs.mkey == 2) {
											
						p1 = pow((MMlsvs.arg3 - p_curr) / MMlsvs.arg3, MMlsvs.arg6) *
						 	 pow((MMlsvs.arg4 - MMlsvs.arg5) / (T_curr - MMlsvs.arg5), MMlsvs.arg7);
						
						D.Area[i][j].mu[l] = MMlsvs.arg2 * pow((MMlsvs.arg1 / MMlsvs.arg2), p1);
					}

					/* Johari & Whalley */
					if (MMlsvs.mkey == 3) {

						p1 = exp(MMlsvs.arg3 * MMlsvs.arg2 / (MMlsvs.arg2 - p_curr));
						
						D.Area[i][j].mu[l] = MMlsvs.arg1 * p1;
					}
					
					/* Doolitle */
					if (MMlsvs.mkey == 4) {

						p1 = 1 + MMlsvs.arg5 * (T_curr - MMlsvs.arg2);	/* Vinf/VinfR */
						p2 = rhoR / D.Area[i][j].rho[l];				/* V/VR = rhoR/rho */
						p3 = p1 / (p2 - MMlsvs.arg4 * p1);

						D.Area[i][j].mu[l] = MMlsvs.arg1 * exp(MMlsvs.arg3 * MMlsvs.arg4 * (p3 - 1 / (1 - MMlsvs.arg4)));
						}

					/* improved Yasutomi */
					if (MMlsvs.mkey == 5) {

							p1 = MMlsvs.arg2 + MMlsvs.arg3 * log(1 + MMlsvs.arg4 * p_curr); /* TG(p) */
							p2 = pow((1 + MMlsvs.arg5 * p_curr), MMlsvs.arg6);				/* F(p) */
							p3 = -2.303 * MMlsvs.arg7 * (T_curr - p1) * p2;
							p4 = MMlsvs.arg8 + (T_curr - p1) * p2;

							D.Area[i][j].mu[l] = MMlsvs.arg1 * exp(p3 / p4);
						}
					}
				
				/* -- THERMAL CONDUCTIVITY -- */
					{
						
					/* constant */
					if (MMthc.mkey == 0) D.Area[i][j].th[l] = MMthc.arg1;
						
					/* pressure-depended */
					if (MMthc.mkey == 1) {
						D.Area[i][j].th[l] = MMthc.arg1 * (1 + MMthc.arg2 * p_curr / (1 + MMthc.arg3 * p_curr));
					}
					
					/* power-law */
					if (MMthc.mkey == 2) {

						p1 = rhoR / D.Area[i][j].rho[l]; /* V/VR = rhoR/rho */
						p2 = T_curr / TR;
						p3 = p1 * (1 + (-0.101) * p2 * pow(p1, 3));

						D.Area[i][j].th[l] = (MMthc.arg1 + MMthc.arg2 * pow(p3, -1.0 * MMthc.arg3));
					}
					
					}

				/* -- HEAT CAPACITY -- */
					{
						
					/* constant */
					if (MMhcp.mkey == 0) {

							D.Area[i][j].hc[l] = MMhcp.arg8;
					}
					
					/* linear */
					if (MMhcp.mkey == 1) {

						p1 = (T_curr / MMhcp.arg7) * pow((rhoR / D.Area[i][j].rho[l]), -4.0); /* V/VR = rhoR/rho */
						D.Area[i][j].hc[l] = (MMhcp.arg1 + MMhcp.arg2 * p1);

					}

					/* pressure-depended */

					if (MMhcp.mkey == 2) {

						p1 = MMhcp.arg3 * (1 + MMhcp.arg4 * p_curr + MMhcp.arg5 * p_curr * p_curr); /* beta(p) */
						p2 = 1 + MMhcp.arg1 * p_curr / (1 + MMhcp.arg2 * p_curr);
						p3 = 1 + p1 * (T_curr - MMhcp.arg6);

						D.Area[i][j].hc[l] = MMhcp.arg7 * p3 * p2;
						}
					}
						
			} /* l loop */


			/* -- SHEAR STRESS -- */
			{
				/* Newton behavior, eta=mu(p,T) */
				if (MMgnvs.mkey == 0) {

					/* pressure gradients calculations */
					if (i == (Nx - 1)) dpdx = 0; /* dpdx is zero at the outlet */
					else dpdx = ph * (D.Area[i + 1][j].p - D.Area[i][j].p) / (dx * a);

					if (j == (Ny - 1)) dpdy = 0; /* dpdy is zero at the outlet */
					else dpdy = ph * (D.Area[i][j + 1].p - D.Area[i][j].p) / (dy * b);

					for (l = 0; l < Nz; ++l) { /* fill vectors for integrate procedure */

						p1 = 1 / D.Area[i][j].mu[l];
						intf1.push_back(p1);		/* 1/eta */
					}
					
					p1 = calc_integrals(intf1, z1_curr, dz, Nz, 1, q); /* integrate 1/eta by dz/ */
					p2 = calc_integrals(intf1, z1_curr, dz, Nz, 2, q); /* integrate z/eta by dz/ */

					rtau_zx0 = ((vel2[0] - vel1[0]) - dpdx * p2) / p1;
					rtau_zy0 = ((vel2[1] - vel1[1]) - dpdy * p2) / p1;

					rtau_mod = mvec(rtau_zx0, rtau_zy0);
					check_LSS(D.SRR[0], D.SRR[1], p_curr, rtau_mod);
					D.Area[i][j].tau[0] = rtau_mod; /* shear stress on surface of body 1*/

					/* calculate tau in each z-point */
					for (l = 0; l < Nz; ++l) {

						z_curr = z1_curr + l * dz; /* current z-coordinate, m */

						p1 = rtau_zx0 + z_curr * dpdx;
						p2 = rtau_zy0 + z_curr * dpdy;
						p3 = mvec(p1, p2);

						D.Area[i][j].tau[l] = p3;
						D.Area[i][j].eta[l] = D.Area[i][j].mu[l];
						
					}
					intf1.clear();
				}

				/* Ree-Eyring */
				if (MMgnvs.mkey == 1) {

					tau0 = MMgnvs.arg1;
					
					/* calculate tau_zx0, tau_zy0 in (i,j)-point */
					e1 = 1; e2 = 1; 
					toler = 1e-4; 
					cnt = 0; cnt_max = 15; /* max iteration counter value */
					
					tau_zx0.push_back(0.0); tau_zy0.push_back(0.0); /* initial values */
						
					/* pressure gradients calculations */ 
					if (i == (Nx - 1)) dpdx = 0; /* dpdx is zero at the outlet */
						else dpdx = ph * (D.Area[i + 1][j].p - D.Area[i][j].p) / (dx * a);

					if (j == (Ny - 1)) dpdy = 0; /* dpdy is zero at the outlet */
						else dpdy = ph * (D.Area[i][j + 1].p - D.Area[i][j].p) / (dy * b);

					/* Newton-Raphson procedure */
					while ((e1 > toler) || (e2 > toler))  {
						
						if (cnt > cnt_max) break; /* avoid for endless loop */
						
						for (l = 0; l < Nz; ++l) {

							z_curr = z1_curr + l * dz; /* current z-coordinate, m */

							taux_curr = tau_zx0[cnt] + z_curr * dpdx;
							tauy_curr = tau_zy0[cnt] + z_curr * dpdy;
							tau_mod = mvec(taux_curr, tauy_curr);

							p1 = tau0 * pow(tauy_curr, 2) * sinh(tau_mod / tau0) +
								tau_mod * pow(taux_curr, 2) * cosh(tau_mod / tau0);

							p2 = tau0 * pow(taux_curr, 2) * sinh(tau_mod / tau0) +
								tau_mod * pow(tauy_curr, 2) * cosh(tau_mod / tau0);

							p3 = tau_mod * cosh(tau_mod / tau0) - tau0 * sinh(tau_mod / tau0);

							if (tau_mod == 0.0) p4 = 0.0; else p4 = 1 / (D.Area[i][j].mu[l] * pow(tau_mod, 3));

							if (tau_mod == 0.0) p5 = 0.0; else p5 = taux_curr * tauy_curr * p4;

							if (tau_mod == 0.0) p6 = 0.0; else p6 = taux_curr * tau0 * sinh(tau_mod / tau0) / (D.Area[i][j].mu[l] * tau_mod);
							if (tau_mod == 0.0) p7 = 0.0; else p7 = tauy_curr * tau0 * sinh(tau_mod / tau0) / (D.Area[i][j].mu[l] * tau_mod);

							intf1.push_back(p6);		/* function in integral for Fx */
							intf2.push_back(p7);		/* function in integral for Fy */
							intf3.push_back(p4 * p1);	/* function in integral for dFx/dtau_zx0 */
							intf4.push_back(p4 * p2);	/* function in integral for dFy/dtau_zy0 */
							intf5.push_back(p5 * p3);	/* function in integral for dFx/dtau_zy0 = dFy/dtau_zx0 */
						
						} /* l loop */
						
						b1  = calc_integrals(intf1, z1_curr, dz, Nz, 1, q) + (vel1[0] - vel2[0]); /* Fx */
						b2  = calc_integrals(intf2, z1_curr, dz, Nz, 1, q) + (vel1[1] - vel2[1]); /* Fy */
						
						a11 = calc_integrals(intf3, z1_curr, dz, Nz, 1, q); /* dFx/dtau_zx0 */
						a22 = calc_integrals(intf4, z1_curr, dz, Nz, 1, q); /* dFy/dtau_zy0 */
						a12 = calc_integrals(intf5, z1_curr, dz, Nz, 1, q); /* dFx/dtau_zy0 = dFy/dtau_zx0 */
						
						if (((a22 == 0) || (a11 == 0)) && (a12 == 0)) {
							p1 = 0.0;
							p2 = 0.0;
						}
						else {
							p1 = (-b1 * a22 + b2 * a12) / (a22 * a11 - a12 * a12);
							p2 = (-b2 * a11 + b1 * a12) / (a22 * a11 - a12 * a12);
						}
						
						double detM = a22 * a11 - a12 * a12; /* determinant of system matrix */

						/* new solution for tau_zx0, tau_zy0 */
						
							tau_zx0.push_back(p1);
							tau_zy0.push_back(p2);
						
							e1 = abs(p1 - tau_zx0[cnt]);
							e2 = abs(p2 - tau_zy0[cnt]);
						
							cnt = cnt + 1;
						
						/* clear vectors */
						intf1.clear(); intf2.clear(); intf3.clear(); intf4.clear(); intf5.clear();
						
					} /* while */

					/* calculate shear stress and general viscosity along z in (i,j)-point */
					rtau_zx0 = tau_zx0[cnt];
					rtau_zy0 = tau_zy0[cnt];
					
					rtau_mod = mvec(rtau_zx0, rtau_zy0);
					
					/* check LSS (according Habchi, 2018) */
					check_LSS(D.SRR[0], D.SRR[1], p_curr, rtau_mod);
					D.Area[i][j].tau[0] = rtau_mod; /* shear stress on surface of body 1*/
														
					for (l = 0; l < Nz; ++l) {
						
						z_curr = z1_curr + l * dz; /* current z-coordinate, m */
						
						p1 = rtau_zx0 + z_curr * dpdx;
						p2 = rtau_zy0 + z_curr * dpdy;
						p3 = 1 / D.Area[i][j].mu[l];
						p4 = mvec(p1, p2);
						if (p4 == 0) p5 = 0; else p5 = (tau0 / p4) * sinh(p4 / tau0);
												
						D.Area[i][j].tau[l] = p4 / ph;
						if (p5 == 0) D.Area[i][j].eta[l] = (1 / p3); 
						 else D.Area[i][j].eta[l] = 1 / (p3 * p5);
					
						/* clear vectors */
						tau_zx0.clear(); tau_zy0.clear();

					} /*l loop */

				} /* Ree-Eyring */

				/* Carreau */
				if (MMgnvs.mkey == 2) {

				}

				/* Carreau-Yasuda (dbl. modified) */
				if (MMgnvs.mkey == 3){
				
				}


			} /*shear stres group */

			//std::cout << "rho[" << i << "][" << j << "][0] = " << D.Area[i][j].rho[0] << '\n';

		} /* j loop */

	} /* i loop */
	
	/* and parameters for all models*/

}

void calc_densviscNB(Domain& D, MM_density MMrho, MM_ls_viscosity MMlsvs)
{

	int		i, j, jend, l, Nx, Ny, Nz;
	double	ph, rhoR, muR, etaR, TR;
	double	a, b, h, dx, dy, dz, z1_curr, z2_curr;
	double	T_curr, p_curr, rho_curr, d1, d2, z_curr;

	double p1, p2, p3, p4, p5, p6, p7;
	
	/* get domain data */
	a = D.a;
	b = D.b;
	Nx = D.Nx;
	Ny = D.Ny;
	dx = D.dx;
	dy = D.dy;
	rhoR = D.rhoR;

	//if (MMlsvs.mkey == 1) D.muR = MMlsvs.arg1;
	muR = D.muR;
	
	TR = D.T0;
	ph = D.ph;

	jend = 1; /* line contact is set by default */
	if (D.dkey == 2) jend = Ny;

	for (i = 0; i < Nx; ++i) {
		for (j = 0; j < jend; ++j) {

			Nz = D.Area[i][j].Nz;

			d1 = D.Area[i][j].d1;	/* m */
			d2 = D.Area[i][j].d2;	/* m */

			z1_curr = D.Area[i][j].z1 - d1;	/* m */
			z2_curr = D.Area[i][j].z2 + d2;	/* m */

			h = (z2_curr - z1_curr) + D.h0;	/* m */
			dz = h / (Nz - 1);		/* m */

			for (l = 0; l < Nz; ++l) {

				T_curr = D.Area[i][j].T[l] * TR;
				p_curr = D.Area[i][j].p * ph;

				// we store DIMENSION values for all parameters

				/* -- DENSITY -- */
				{
					/* this is a common terms for Tait and Murnaghan eqns */
					p1 = (1 + MMrho.arg4 * (T_curr - MMrho.arg2));
					p2 = 1 / (1 + MMrho.arg5);
					p5 = 1 / MMrho.arg5;

					switch (MMrho.mkey)
					{
					case 0: { // constant

						D.Area[i][j].rho[l] = MMrho.arg1;

						break;
					} //case 0

					case 1: {// Dowson & Higginson

						D.Area[i][j].rho[l] = MMrho.arg1 * (1 + MMrho.arg2 * p_curr / (1 + MMrho.arg3 * p_curr));

						break;
					} //case 1

					case 2: { // Tait

						p3 = 1 + (1 + MMrho.arg5) * p_curr / (MMrho.arg6 * exp(-1.0 * MMrho.arg3 * T_curr));

						if (p_curr == 0) D.Area[i][j].rho[l] = MMrho.arg1;
						else D.Area[i][j].rho[l] = MMrho.arg1 / (p1 * (1 - p2 * log(p3)));

						break;
					} //case 2

					case 3: { // Murnaghan

						p4 = 1 + MMrho.arg5 * p_curr / (MMrho.arg6 * exp(-1.0 * MMrho.arg4 * T_curr));

						if (p_curr == 0) D.Area[i][j].rho[l] = MMrho.arg1;
						else D.Area[i][j].rho[l] = (MMrho.arg1 / p1) * pow(p4, p5);

						break;
					} //case 3

					case 4: {// Dowson & Higginson var 1

						D.Area[i][j].rho[l] = MMrho.arg1 * (MMrho.arg2 + MMrho.arg3 * p_curr / (MMrho.arg2 + p_curr));

						break;
					} //case 1

					}
				}

				/* -- LS-VISCOSITY -- */
				{
					switch (MMlsvs.mkey)
					{
					case 0: { // constant

						D.Area[i][j].mu[l] = MMlsvs.arg1;
						break;
					} // case 0

					case 1: { // Barus

						D.Area[i][j].mu[l] = MMlsvs.arg1 * exp(MMlsvs.arg2 * p_curr);
						break;
					} // case 1

					case 2: { // Roeland

						p1 = pow((MMlsvs.arg3 - p_curr) / MMlsvs.arg3, MMlsvs.arg6) *
							pow((MMlsvs.arg4 - MMlsvs.arg5) / (T_curr - MMlsvs.arg5), MMlsvs.arg7);

						D.Area[i][j].mu[l] = MMlsvs.arg2 * pow((MMlsvs.arg1 / MMlsvs.arg2), p1);
						break;
					} // case 2

					case 3: {// Johari & Whalley

						p1 = exp(MMlsvs.arg3 * MMlsvs.arg2 / (MMlsvs.arg2 - p_curr));

						D.Area[i][j].mu[l] = MMlsvs.arg1 * p1;
						break;
					} // case 3

					case 4: {// Doolitle

						p1 = 1 + MMlsvs.arg5 * (T_curr - MMlsvs.arg2);	/* Vinf/VinfR */
						p2 = rhoR / D.Area[i][j].rho[l];				/* V/VR = rhoR/rho */
						p3 = p1 / (p2 - MMlsvs.arg4 * p1);

						D.Area[i][j].mu[l] = MMlsvs.arg1 * exp(MMlsvs.arg3 * MMlsvs.arg4 *
							(p3 - 1 / (1 - MMlsvs.arg4)));
						break;
					} // case 4

					case 5: { // improved Yasutomi

						p1 = MMlsvs.arg2 + MMlsvs.arg3 * log(1 + MMlsvs.arg4 * p_curr); /* TG(p) */
						p2 = pow((1 + MMlsvs.arg5 * p_curr), MMlsvs.arg6);				/* F(p) */
						p3 = -2.303 * MMlsvs.arg7 * (T_curr - p1) * p2;
						p4 = MMlsvs.arg8 + (T_curr - p1) * p2;

						D.Area[i][j].mu[l] = MMlsvs.arg1 * exp(p3 / p4);
						break;
					} // case 5


					} // switch
				}

				/* -- GEN-VISCOSITY (NEWTON BEHAVIOR) -- */
				D.Area[i][j].eta[l] = D.Area[i][j].mu[l];


				//std::cout << "rho[" << i << "][" << j << "][0] = " << D.Area[i][j].rho[0] << '\n';
			
			} // l loop
		} /* j loop */

	} /* i loop */

} // proc



void calc_fric_1D_old(Domain& D, Body &B1, Body &B2) {

	dbvector inveta, inveta_n;

	int		 i, j, l, ll, Nx, Ny, Nz;
	int		j_end;

	double	dx, dy, h, RX, z1, z2, dz, z_curr, res;
	double	eta_e, rho_e, eta_e1, rho_1, rho_2;
	double int1, int2, dpdx, us, tau_sum0, tau_sum1;
	double p_curr, rtauL;

	Nx = D.Nx;
	Ny = D.Ny;

	dx = D.dx;
	dy = D.dy;

	tau_sum0 = 0.0; tau_sum1 = 0.0;
	for (i = 1; i < Nx; ++i) {

		Nz = D.Area[i][0].Nz;
		h  = D.Area[i][0].h;
		// all parameters below is in dimension form
		z1 = (D.Area[i][0].z1 - D.Area[i][0].d1);	// start z-coord
		z2 = (D.Area[i][0].z2 + D.Area[i][0].d2 + D.h0);	// end z-coord
		dz = abs(z2 - z1) / (Nz - 1); // each [i,j]-point has its own dz
		p_curr = D.Area[i][0].p;

		//rtauL = calc_LSS(D.SRR[0], D.SRR[1], p_curr) * D.ph;
		rtauL = 1e9;
			
		// calculate eta, eta1 (dimension form)
		for (l = 0; l < Nz; ++l) 
			inveta.push_back(1 / D.Area[i][0].eta[l]);
			
		eta_e  = 1 / calc_integrals(inveta, z1, dz, Nz, 1, 0);
		eta_e1 = 1 / calc_integrals(inveta, z1, dz, Nz, 2, 0);
			
		// calculate velocity and shear stress in each z-position
		D.Area[i][0].u[0] = B1.vel[0]; // body surface velocity
		us = B2.vel[0] - B1.vel[0];

		// shear stress on 1st body's surface (Newtonian behavior)
		D.Area[i][0].tau[0] = D.ph * ((D.Area[i][0].p - D.Area[i - 1][0].p) / dx) * 
			(z1 - (eta_e / eta_e1)) + eta_e * us;

		/* check LSS */
		if (D.Area[i][0].tau[0] > rtauL) D.Area[i][0].tau[0] = rtauL;

		tau_sum0 = tau_sum0 + D.Area[i][0].tau[0];
		
		//std::cout << D.Area[i][0].tau[0] << "\n";

		inveta_n.push_back(inveta[0]);
					
		for (l = 1; l < Nz; ++l) {
				
			z_curr = z1 + l * dz;
			inveta_n.push_back(inveta[l]);
			int1 = calc_integrals(inveta_n, z1, dz, l + 1, 1, 0); // integral(1/eta, z=z1..z)
			int2 = calc_integrals(inveta_n, z1, dz, l + 1, 2, 0); // integral(z/eta, z=z1..z)
			
			dpdx = D.ph * ((D.Area[i][0].p - D.Area[i - 1][0].p) / dx); // pressure gradient (x-comp)
				
			D.Area[i][0].u[l] = D.Area[i][0].u[0] + dpdx * (int2 - (eta_e / eta_e1) * int1) +
			eta_e * us * int1;
				
			// Newtonian behavior
			D.Area[i][0].tau[l] = dpdx * (z_curr - (eta_e / eta_e1)) + eta_e * us;
			
			/* check LSS */
			if (D.Area[i][0].tau[l] > rtauL) D.Area[i][0].tau[l] = rtauL;

			// shear stress on 2nd body's surface (Newtonian behavior)
			if (l == Nz - 1) {
			
				D.Area[i][0].tau[l] = dpdx * (z2 - (eta_e / eta_e1)) + eta_e * us;
				/* check LSS */
				if (D.Area[i][0].tau[l] > rtauL) D.Area[i][0].tau[l] = rtauL;

				tau_sum1 = tau_sum1 + D.Area[i][0].tau[l];
							
			}
		} // l loop

		inveta.clear();
		inveta_n.clear();
		
	} // i loop

	// friction forces on bodies surfaces and their friction coefficients
	B1.Ft = abs(tau_sum0 * (D.a * dx));
	B2.Ft = abs(tau_sum1 * (D.a * dx));

	B1.fc = B1.Ft / D.W;
	B2.fc = B2.Ft / D.W;

} // proc

void calc_fric_1D(Domain& D, Body& B1, Body& B2) {

	dbvector inveta, inveta_n;

	int		 i, j, l, ll, Nx, Ny, Nz;
	int		j_end;

	double	dx, dy, h, RX, z1, z2, dz, z_curr, res;
	double	eta_e, rho_e, eta_e1, rho_1, rho_2;
	double int0, int1, int2, dpdx, us, tau0, tau_sum, tau_sum0, tau_sum1;
	double p_curr, rtauL, tau_curr, u1, u2;

	std::fstream extfile;

	Nx = D.Nx;
	Ny = D.Ny;

	dx = D.dx;
	dy = D.dy;

	u1 = B1.vel[0];
	u2 = B2.vel[0];
	us = B2.vel[0] - B1.vel[0];

	extfile.open("lc_tau0.txt", std::ios::out);

	tau_sum0 = 0.0; tau_sum1 = 0.0;
	for (i = 1; i < Nx; ++i) {

		Nz = D.Area[i][0].Nz;
		h = D.Area[i][0].h;
		z1 = 0.0;
		z2 = (D.Area[i][0].z2 + D.Area[i][0].d2 + D.h0);	// end z-coord
		dz = abs(z2 - z1) / (Nz - 1); // each [i,j]-point has its own dz
		
		p_curr = D.Area[i][0].p * D.ph; // [Pa]
		dpdx = D.ph * ((D.Area[i][0].p - D.Area[i - 1][0].p) / dx); // pressure gradient, [Pa]

		for (l = 0; l < Nz; ++l)
			inveta.push_back(1 / D.Area[i][0].eta[l]);

		eta_e  = calc_integrals(inveta, 0.0, dz, Nz, 1, 0);
		eta_e1 = calc_integrals(inveta, 0.0, dz, Nz, 2, 0);

		// shear stress on 1st body's surface (Newtonian behavior)
		tau0 = (u2 - u1 - dpdx * eta_e1) / eta_e;
		check_LSS(D.SRR[0], D.SRR[1], p_curr, tau0);
		D.Area[i][0].tau[0] = tau0;

		extfile << (-1.0 * D.kx) + i * D.dx << "\t" << tau0 << "\n";

		D.Area[i][0].u[0] = u1;

		inveta_n.push_back(inveta[0]);
		for (l = 1; l < Nz; ++l) {

			z_curr = z1 + l * dz;
			
			inveta_n.push_back(inveta[l]);
			int1 = calc_integrals(inveta_n, z1, dz, l + 1, 1, 0); // integral(1/eta, z=z1..z)
			int2 = calc_integrals(inveta_n, z1, dz, l + 1, 2, 0); // integral(z/eta, z=z1..z)

			D.Area[i][0].u[l] = D.Area[i][0].u[0] + dpdx * (int2 - (eta_e / eta_e1) * int1) +
				eta_e * us * int1;
						
			tau_curr = tau0 + z_curr * dpdx;
			check_LSS(D.SRR[0], D.SRR[1], p_curr, tau_curr);
			D.Area[i][0].tau[l] = tau_curr;

		} // l loop

		tau_sum0 = tau_sum0 + tau0;
		tau_sum1 = tau_sum1 + D.Area[i][0].tau[Nz - 1];

		inveta.clear();
	} // i loop

	// friction forces on bodies surfaces and their friction coefficients
	B1.Ft = abs(tau_sum0 * (D.a * dx));
	B2.Ft = abs(tau_sum1 * (D.a * dx));

	B1.fc = B1.Ft / D.W;
	B2.fc = B2.Ft / D.W;

	extfile.close();

} // proc


void calc_fric_2D_old(Domain& D, Body& B1, Body& B2)
{
	dbvector inveta, inveta_n;
	int		 i, j, l, ll, Nx, Ny, Nz;
	int		j_end;

	double	dx, dy, h, RX, z1, z2, dz, z_curr, res;
	double	eta_e, rho_e, eta_e1, rho_1, rho_2;
	double int1, int2, dpdx, dpdy, us, vs, tau_sum0, tau_sum1;
	double tauxz0, tauyz0, tauxz, tauyz;

	Nx = D.Nx;
	Ny = D.Ny;

	dx = D.dx;
	dy = D.dy;


	tau_sum0 = 0.0; tau_sum1 = 0.0;
	for (i = 1; i < Nx; ++i) {
		for (j = 1; j < Ny; ++j) {
			
			Nz = D.Area[i][j].Nz;
			h = D.Area[i][j].h;
		
			// all parameters and integrals below is in dimension form
			z1 = (D.Area[i][j].z1 - D.Area[i][j].d1);	// start z-coord
			z2 = (D.Area[i][j].z2 + D.Area[i][j].d2 + D.h0);	// end z-coord
			dz = abs(z2 - z1) / (Nz - 1); // each [i,j]-point has its own dz

			dpdx = D.ph * ((D.Area[i][j].p - D.Area[i - 1][j].p) / dx); // pressure gradient (x-comp)
			dpdy = D.ph * ((D.Area[i][j].p - D.Area[i][j - 1].p) / dy); // pressure gradient (y-comp)

			// calculate eta, eta1 (dimension form)
			for (l = 0; l < Nz; ++l)
				inveta.push_back(1 / D.Area[i][j].eta[l]);

			eta_e  = 1 / calc_integrals(inveta, z1, dz, Nz, 1, 0);
			eta_e1 = 1 / calc_integrals(inveta, z1, dz, Nz, 2, 0);

			// calculate velocity and shear stress in each z-position
			D.Area[i][j].u[0] = B1.vel[0]; // body surface velocity (x-comp)
			D.Area[i][j].v[0] = B1.vel[1]; // body surface velocity (y-comp)
			us = B2.vel[0] - B1.vel[0];
			vs = B2.vel[1] - B1.vel[1];

			// shear stress on 1st body's surface (Newtonian behavior)
			tauxz0 = dpdx * (z1 - (eta_e / eta_e1)) + eta_e * us;
			tauyz0 = dpdy * (z1 - (eta_e / eta_e1)) + eta_e * vs;

			D.Area[i][j].tau[0] = pow(tauxz0 * tauxz0 + tauyz0 * tauyz0,0.5);
			tau_sum0 = tau_sum0 + D.Area[i][j].tau[0];

			//std::cout << D.Area[i][0].tau[0] << "\n";

			inveta_n.push_back(inveta[0]);

			for (l = 1; l < Nz; ++l) {

			z_curr = z1 + l * dz;
			inveta_n.push_back(inveta[l]);
			int1 = calc_integrals(inveta_n, z1, dz, l + 1, 1, 0); // integral(1/eta, z=z1..z)
			int2 = calc_integrals(inveta_n, z1, dz, l + 1, 2, 0); // integral(z/eta, z=z1..z)

			
			D.Area[i][j].u[l] = D.Area[i][j].u[0] + dpdx * (int2 - (eta_e / eta_e1) * int1) +
				eta_e * us * int1;
			D.Area[i][j].v[l] = D.Area[i][j].v[0] + dpdy * (int2 - (eta_e / eta_e1) * int1) +
				eta_e * vs * int1;

			// Newtonian behavior
			tauxz = dpdx * (z_curr - (eta_e / eta_e1)) + eta_e * us;
			tauyz = dpdy * (z_curr - (eta_e / eta_e1)) + eta_e * vs;

			D.Area[i][j].tau[l] = pow(tauxz * tauxz + tauyz * tauyz, 0.5);
			
			
			// shear stress on 2nd body's surface (Newtonian behavior)
			if (l == Nz - 1) {

				tauxz = D.ph * ((D.Area[i][j].p - D.Area[i - 1][j].p) / dx) *
					(z2 - (eta_e / eta_e1)) + eta_e * us;
				tauyz = D.ph * ((D.Area[i][j].p - D.Area[i][j - 1].p) / dy) *
					(z2 - (eta_e / eta_e1)) + eta_e * vs;

				D.Area[i][j].tau[l] = pow(tauxz * tauxz + tauyz * tauyz, 0.5);
				tau_sum1 = tau_sum1 + D.Area[i][j].tau[l];

			}
		} // l loop

			inveta.clear();
			inveta_n.clear();
		} // j loop
	} // i loop

// friction forces on bodies surfaces and their friction coefficients
B1.Ft = abs(tau_sum0 * (D.a * dx) * (D.b * dy));
B2.Ft = abs(tau_sum1 * (D.a * dx) * (D.b * dy));

B1.fc = B1.Ft / D.W;
B2.fc = B2.Ft / D.W;

} // proc

void calc_fric_2D(Domain& D, Body& B1, Body& B2)
{
	dbvector inveta, inveta_n;
	int		 i, j, l, ll, Nx, Ny, Nz;
	int		j_end;

	double	dx, dy, h, RX, z1, z2, dz, z_curr, res;
	double	eta_e, rho_e, eta_e1, rho_1, rho_2;
	double  int1, int2, dpdx, dpdy, us, vs, tau_sum0, tau_sum1;
	double tauxz0, tauyz0, tauxz, tauyz;
	double u1, u2, v1, v2;
	double p_curr, tau_curr, tau0;

	std::fstream extfile;

	Nx = D.Nx;
	Ny = D.Ny;

	dx = D.dx;
	dy = D.dy;

	u1 = B1.vel[0];
	u2 = B2.vel[0];
	us = B2.vel[0] - B1.vel[0];
	
	v1 = B1.vel[1];
	v2 = B2.vel[1];
	vs = B2.vel[1] - B1.vel[1];

	extfile.open("pc_tau0.txt", std::ios::out);

	tau_sum0 = 0.0; tau_sum1 = 0.0;
	for (i = 1; i < Nx; ++i) {
		for (j = 1; j < Ny; ++j) {

			Nz = D.Area[i][j].Nz;
			h = D.Area[i][j].h;
			p_curr = D.Area[i][j].p * D.ph;

			// all parameters and integrals below is in dimension form
			//z1 = (D.Area[i][j].z1 - D.Area[i][j].d1);	// start z-coord
			z1 = 0.0; //since body 1 is rigid (for reduced geometry)
			z2 = (D.Area[i][j].z2 + D.Area[i][j].d2 + D.h0);	// end z-coord
			dz = abs(z2 - z1) / (Nz - 1); // each [i,j]-point has its own dz

			dpdx = D.ph * ((D.Area[i][j].p - D.Area[i - 1][j].p) / dx); // pressure gradient (x-comp)
			dpdy = D.ph * ((D.Area[i][j].p - D.Area[i][j - 1].p) / dy); // pressure gradient (y-comp)

			for (l = 0; l < Nz; ++l)
				inveta.push_back(1 / D.Area[i][j].eta[l]);

			eta_e  = 1 / calc_integrals(inveta, z1, dz, Nz, 1, 0);
			eta_e1 = 1 / calc_integrals(inveta, z1, dz, Nz, 2, 0);
	
			// shear stress on 1st body's surface (Newtonian behavior)
			tauxz0 = (u2 - u1 - dpdx * eta_e1) / eta_e;
			tauyz0 = (v2 - v1 - dpdy * eta_e1) / eta_e;
			tau0 = pow(tauxz0 * tauxz0 + tauyz0 * tauyz0, 0.5);
			
			/* check LSS (according Habchi, 2018) */
			check_LSS(D.SRR[0], D.SRR[1], p_curr, tau0);
			D.Area[i][j].tau[0] = tau0;

			extfile << (-1.0 * D.kx) + i * D.dx << "\t" << (-1.0 * D.ky) + j * D.dy << "\t" << tau0 << "\n";
			
			tau_sum0 = tau_sum0 + tau0;

			D.Area[i][j].u[0] = u1;
			D.Area[i][j].v[0] = v1;

			inveta_n.push_back(inveta[0]);
			for (l = 1; l < Nz; ++l) {

				z_curr = z1 + l * dz;
				
				inveta_n.push_back(inveta[l]);
				int1 = calc_integrals(inveta_n, z1, dz, l + 1, 1, 0); // integral(1/eta, z=z1..z)
				int2 = calc_integrals(inveta_n, z1, dz, l + 1, 2, 0); // integral(z/eta, z=z1..z)

				D.Area[i][j].u[l] = D.Area[i][j].u[0] + dpdx * (int2 - (eta_e / eta_e1) * int1) +
					eta_e * us * int1;
				D.Area[i][j].v[l] = D.Area[i][j].v[0] + dpdy * (int2 - (eta_e / eta_e1) * int1) +
					eta_e * vs * int1;

				tauxz = z_curr * dpdx + tauxz0;
				tauyz = z_curr * dpdy + tauyz0;
				tau_curr = pow(tauxz * tauxz + tauyz * tauyz, 0.5);
				
				/* check LSS (according Habchi, 2018) */
				check_LSS(D.SRR[0], D.SRR[1], p_curr, tau_curr);
				D.Area[i][j].tau[l] = tau_curr;

			} // l loop

			tau_sum1 = tau_sum1 + D.Area[i][j].tau[Nz-1];

			inveta.clear(); inveta_n.clear();
				
		} // j loop
	} // i loop

	extfile.close();

// friction forces on bodies surfaces and their friction coefficients
	B1.Ft = abs(tau_sum0 * (D.a * dx) * (D.b * dy));
	B2.Ft = abs(tau_sum1 * (D.a * dx) * (D.b * dy));

	B1.fc = B1.Ft / D.W;
	B2.fc = B2.Ft / D.W;

} // proc