#include "main.h"
#include "fftw3.h" /* FFTW3 library for Fourier transformations */


double IC_F(double x, double y)
{
/* local exoression in calc_IC */

	return(x + sqrt(x * x + y * y));

} /* proc */

double calc_IC_1D(double ki, double delta)
{
	/************************** INFLUENCE COEFFICIENTS CALCULATION ***************/
		// 1D domain, equidistant grid points, 0-order for shape function Ys

	double xki, a, res;
	
	xki = ki * delta;
	
	a = 0.5 * delta;
	
	res =(xki + a) * log(abs(xki + a)) - (xki - a) * log(abs(xki - a)) - 2 * delta;
	
	return (res);

} /* proc */

double calc_IC_2D(double ki, double lj, double delta)
{
/************************** INFLUENCE COEFFICIENTS CALCULATION ***************/
	// 2D domain, equidistant grid points, 0-order for shape function Ys
	
	double xki, ylj, a, b;
	double p1, p2, p3, p4, res;

	xki = ki * delta;
	ylj = lj * delta;

	a = 0.5 * delta;
	b = 0.5 * delta;

	p1 = (xki + a) * log(IC_F(ylj + b, xki + a) / IC_F(ylj - b, xki + a));
	p2 = (ylj + b) * log(IC_F(xki + a, ylj + b) / IC_F(xki - a, ylj + b));
	p3 = (xki - a) * log(IC_F(ylj - b, xki - a) / IC_F(ylj + b, xki - a));
	p4 = (ylj - b) * log(IC_F(xki - a, ylj - b) / IC_F(xki + a, ylj - b));

	res = p1 + p2 + p3 + p4;

	return (res);

} /* proc */

double calc_roughness()
{

	return(0);
} /* proc */

void calc_deformDS(Domain& D, Body B1, Body B2)
{
	int		k, i, Nx;
	double	summ, dx, C1, C2, ph, a, Estar;

	Nx = D.Nx;
	dx = D.dx;
	ph = D.ph;
	a = D.a;

	for (k = 0; k < Nx; ++k) {
		summ = 0.0;
		for (i = 0; i < Nx; ++i) 
			summ = summ + calc_IC_1D(k - i, dx) * D.Area[i][0].p;
				
		/* store deformations to domain */
		//C1 = 2 / (pi * B1.YMe);
		//C2 = 2 / (pi * B2.YMe);

		/* equvivalent (for reduced case) YM */
		Estar = pow(0.5 * ((1 - B1.PR * B1.PR) / B1.YM + (1 - B2.PR * B2.PR) / B2.YM), -1);
		C2 = 4 / (pi * Estar);
		
		/* new deformations and film thickness in dimension form */
		/* we store the ABSOLUTE VALUE of deformations, without sign! */
		//D.Area[k][0].d1 = abs(a * C1 * ph * summ); // m 
		D.Area[k][0].d1 = 0.0;
		D.Area[k][0].d2 = abs(a * C2 * ph * summ); // m 
	
	} // k loop

} //proc

void calc_deform(Domain& D, Body B1, Body B2)
{
	/************************** ELASTIC DEFORMATIONS CALCULATION ***************/

	fftw_complex* cpex, * cwaic;
	fftw_complex* cpex_cup, * cwaic_cup, * cmult, * cdspl;

	fftw_plan	ddft1, ddft2, idft;

	int		i, j, Nx, Ny, flag10, flag11, strt,dkey;
	double	p_curr, waic_curr, dx, dy, a, b, kx, ky;
	double	shift01, shift10, shift11;
	double	Re_pcup, Im_pcup, Re_wcup, Im_wcup;
	double	C1, C2, ph, Estar;

	Nx = D.Nx;
	Ny = D.Ny;
	dx = D.dx;
	dy = D.dy;
	a  = D.a;
	b  = D.b;
	kx = D.kx;
	ky = D.ky;
	ph = D.ph;
	dkey = D.dkey;

	switch (dkey) {

	case 1: {
			
		/* dynamic arrays initialization */
		cpex  = new fftw_complex[2 * Nx - 1];
		cwaic = new fftw_complex[2 * Nx - 1];

		for (i = 0; i < (2 * Nx - 1); ++i) {

			/* write pressure (zero-padding in extended part) */
			if (i < Nx) p_curr = D.Area[i][0].p;
			else p_curr = 0;

			cpex[i][0] = p_curr;
			cpex[i][1] = 0.0; /* imaginary part */

			/* calc current WAIC-matrix element and write it */
			if (i < Nx) {
				cwaic[i][0] = calc_IC_1D(-i, dx);
				cwaic[i][1] = 0.0; /* imaginary part */
			}
			else {
				cwaic[i][0] = calc_IC_1D(-2 * Nx + 1 + i, dx);
				cwaic[i][1] = 0.0; /* imaginary part */
			}
			//std::cout << "cwaic[" << i << "] = " << cwaic[i][0] << "\n";
			//std::cout << "cpex[" << i << "] = " << cpex[i][0] << "\n";
		} //i loop

		//std::cout << "\n\n";

		/* make DDFT */
		cpex_cup  = new fftw_complex[2 * Nx - 1];
		cwaic_cup = new fftw_complex[2 * Nx - 1];
		cmult	  = new fftw_complex[2 * Nx - 1];
		cdspl	  = new fftw_complex[2 * Nx - 1];

		ddft1 = fftw_plan_dft_1d((2 * Nx - 1), cpex, cpex_cup, FFTW_FORWARD, FFTW_ESTIMATE);
		fftw_execute(ddft1);

		ddft2 = fftw_plan_dft_1d((2 * Nx - 1), cwaic, cwaic_cup, FFTW_FORWARD, FFTW_ESTIMATE);
		fftw_execute(ddft2);

		/* normalization and complex multiplication */
		double norm = 1.0 / pow(2 * Nx - 1,0.5); // symmetric normalization
		//double norm = 1.0;
		for (i = 0; i < (2 * Nx - 1); ++i) {

			Re_pcup = cpex_cup[i][0] * norm;
			Im_pcup = cpex_cup[i][1] * norm;
			Re_wcup = cwaic_cup[i][0] * norm;
			Im_wcup = cwaic_cup[i][1] * norm;

			//std::cout << "cpex_cup[" << i << "] = " << cpex_cup[i][0] << "\n";
			//std::cout << "cwaic_cup[" << i << "] = " << cwaic_cup[i][0] << "\n";
			
			cmult[i][0] = (Re_pcup * Re_wcup - Im_pcup * Im_wcup);
			cmult[i][1] = (Re_pcup * Im_wcup + Im_pcup * Re_wcup);
		}

		/* make IDFT */
		idft = fftw_plan_dft_1d((2 * Nx - 1), cmult, cdspl, FFTW_BACKWARD, FFTW_ESTIMATE);
		fftw_execute(idft);

		/* store deformations to domain */
				
		//Estar = pow(0.5 * ((1 - B1.PR * B1.PR) / B1.YM + (1 - B2.PR * B2.PR) / B2.YM), -1);
		Estar = D.rYm;
		C2 = 4 / (pi * Estar);
				
		for (i = 0; i < Nx; ++i) {
			
			/* new deformations and film thickness in dimension form */
			/* we store the ABSOLUTE VALUE of deformations, without sign! */
			//D.Area[i][0].d1 = abs(a * C1 * ph * cdspl[i][0]); /* real part, m */
			D.Area[i][0].d1 = 0.0;
			D.Area[i][0].d2 = abs(a * C2 * ph * cdspl[i][0]); /* real part, m */
		
			//std::cout << "d1[" << i << "] = " << D.Area[i][0].d1 << "\n";
		}

		/* free memory */
		delete[] cpex_cup;
		delete[] cwaic_cup;
		delete[] cmult;
		delete[] cdspl;
		delete[] cpex;
		delete[] cwaic;
		
		break;
	} //case
	case 2: { /* point contact */

		/* dynamic arrays initialization */
		cpex = new fftw_complex[(2 * Nx) * (2 * Ny)];
		cwaic = new fftw_complex[(2 * Nx) * (2 * Ny)];

		shift11 = 0; shift10 = 0;
		flag11 = 0; flag10 = 0;

		for (i = 0; i < (2 * Nx); ++i) {

			shift01 = 0;
			for (j = 0; j < (2 * Ny); ++j) {

				/* write pressure (zero-padding in extended part) */
				if ((i < Nx) && (j < Ny)) p_curr = D.Area[i][j].p;
				else p_curr = 0;

				cpex[(2 * Ny) * i + j][0] = p_curr;
				cpex[(2 * Ny) * i + j][1] = 0.0; /* imaginary part*/

				/* calc current WAIC-matrix element and write it */
				if ((i < Nx) && (j < Ny)) {
					waic_curr = calc_IC_2D(i, j, dx);
				}
				else if ((i < Nx) && (Ny <= j)) {
					waic_curr = calc_IC_2D(i, -Nx + shift01, dx);
					shift01 = shift01 + 1;
				}
				else if ((Nx <= i) && (j < Ny)) {
					waic_curr = calc_IC_2D(j, -Ny + shift10, dx);
					flag10 = 1;
				}
				else if ((Nx <= i) && (Ny <= j)) {
					waic_curr = calc_IC_2D(-Nx + shift01, -Ny + shift11, dx);
					flag11 = 1;
					shift01 = shift01 + 1;
				}

				cwaic[(2 * Ny) * i + j][0] = waic_curr;
				cwaic[(2 * Ny) * i + j][1] = 0.0; /* imaginary part*/

			} /* j loop */

			if (flag11 == 1) shift11 = shift11 + 1;
			if (flag10 == 1) shift10 = shift10 + 1;

		} /* i loop */

		/* make DDFT */
		cpex_cup  = new fftw_complex[(2 * Nx) * (2 * Ny)];
		cwaic_cup = new fftw_complex[(2 * Nx) * (2 * Ny)];
		cmult	  = new fftw_complex[(2 * Nx) * (2 * Ny)];
		cdspl	  = new fftw_complex[(2 * Nx) * (2 * Ny)];

		ddft1 = fftw_plan_dft_2d((2 * Nx), (2 * Ny), cpex, cpex_cup, FFTW_FORWARD, FFTW_ESTIMATE);
		fftw_execute(ddft1);

		ddft2 = fftw_plan_dft_2d((2 * Nx), (2 * Ny), cwaic, cwaic_cup, FFTW_FORWARD, FFTW_ESTIMATE);
		fftw_execute(ddft2);

		/* normalization and complex multiplication */
		double norm = 1.0 / pow(2 * Nx * 2 * Ny, 0.5); // symmetric normalization
		for (i = 0; i < (2 * Nx) * (2 * Ny); ++i) {

			Re_pcup = cpex_cup[i][0] * norm;
			Im_pcup = cpex_cup[i][1] * norm;
			Re_wcup = cwaic_cup[i][0] * norm;
			Im_wcup = cwaic_cup[i][1] * norm;

			cmult[i][0] = (Re_pcup * Re_wcup - Im_pcup * Im_wcup);
			cmult[i][1] = (Re_pcup * Im_wcup + Im_pcup * Re_wcup);
			
			//std::cout << "cpex_cup[" << i << "] = " << cpex_cup[i][0] << "\n";
			//std::cout << "cwaic_cup[" << i << "] = " << cwaic_cup[i][0] << "\n";
		}

		/* make IDFT */
		idft = fftw_plan_dft_2d((2 * Nx), (2 * Ny), cmult, cdspl, FFTW_BACKWARD, FFTW_ESTIMATE);
		fftw_execute(idft);

		/* store deformations to domain */
			
		//Estar = pow(0.5 * ((1 - B1.PR * B1.PR) / B1.YM + (1 - B2.PR * B2.PR) / B2.YM), -1);
		Estar = D.rYm;
		C2 = 2 / (pi * Estar);

		strt = 0;
		for (i = 0; i < Nx; ++i) {
			for (j = 0; j < Ny; ++j) {

				/* new deformations and film thickness in dimension form */
				//D.Area[i][j].d1 = abs(a * C1 * ph * cdspl[strt + j][0]); /* real part, m */
				D.Area[i][j].d1 = 0.0;
				D.Area[i][j].d2 = abs(a * C2 * ph * cdspl[strt + j][0]); /* real part, m */

				if (j == (Ny - 1)) strt = strt + 2 * Ny;
			}
		}

		/* free memory */
		delete[] cpex_cup;
		delete[] cwaic_cup;
		delete[] cmult;
		delete[] cdspl;
		delete[] cpex;
		delete[] cwaic;

		break;
	} //case

	} //switch

} /* proc */