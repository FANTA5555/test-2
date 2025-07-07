#ifndef MAIN_H
#define MAIN_H

#include <iostream>
#include <fstream>
#include <cmath> 
#include <vector>
#include <string>  

double const pi = acos(-1);			/* Pi number */

typedef std::vector<double>		dbvector;
typedef std::string				string;


struct Point {
	/************* Domian point (i,j) with TEHL-model data ****************/

	int		 Nz;				/* point number along z axis (z-resolution)	  */
	double   z1_n;				/* z coordinate of 1st (bottom) body		  */
	double   z2_n;				/* z coordinate of 2nd (top) body			  */
	double	 z1;				/* z coordinate of 1st body (reduced geom)	  */
	double	 z2;				/* z coordinate of 2nd body (reduced geom)	  */
	double	 h;					/* film thickness, m						  */
	double	 d1;				/* surface deformation for 1st body, m		  */
	double	 d2;				/* surface deformation for 2nd body, m		  */
	double	 p;					/* pressure									  */
	double	 pold;				/* pressure at previous iteration			  */
	double	 dp;				/* pressure addition (uses in point contact)  */
	dbvector T;					/* temperature								  */
	dbvector Told;				/* temperature at previous iteration		  */
	dbvector rho;				/* density									  */
	dbvector mu;				/* low-shear viscosity						  */
	dbvector eta;				/* generalized viscosity					  */
	dbvector tau;				/* shear stress								  */
	dbvector th;				/* thermal conductivity						  */
	dbvector hc;				/* heat capacity							  */
	dbvector u;					/* x component of lubricant velocity		  */
	dbvector v;					/* y component of lubricant velocity		  */
};

struct Drefine{
	/************* Refinement data for TEHL domains ****************/
	double	x1;
	double	x2;
	double	y1;
	double	y2;
	int		val;
};

struct MM_density{

	int		mkey;
	double	arg1;
	double	arg2;
	double	arg3;
	double	arg4;
	double	arg5;
	double	arg6;
	string	IT1;
	string	narg1;
	string	narg2;
	string	narg3;
	string	narg4;
	string	narg5;
	string	narg6;

};

struct MM_ls_viscosity {

	int		mkey;
	double	arg1;
	double	arg2;
	double	arg3;
	double	arg4;
	double	arg5;
	double	arg6;
	double	arg7;
	double	arg8;
	string	IT1;
	string	narg1;
	string	narg2;
	string	narg3;
	string	narg4;
	string	narg5;
	string	narg6;
	string	narg7;
	string	narg8;
};

struct MM_gn_viscosity {

	int		mkey;
	double	arg1;
	double	arg2;
	double	arg3;
	double	arg4;
	double	arg5;
	double	arg6;
	double	arg7;
	string	IT1;
	string	narg1;
	string	narg2;
	string	narg3;
	string	narg4;
	string	narg5;
	string	narg6;
	string	narg7;
	
};

struct MM_th_conductivity {

	int		mkey;
	double	arg1;
	double	arg2;
	double	arg3;
	double	arg4; /* TR */
	double	arg5; /* k (const) */
	string	IT1;
	string	narg1;
	string	narg2;
	string	narg3;
	string	narg4;
	string	narg5;

};

struct MM_ht_capacity {

	int		mkey;
	double	arg1;
	double	arg2;
	double	arg3;
	double	arg4;
	double	arg5;
	double	arg6;
	double	arg7; /* TR */
	double	arg8; /* Cp (const) */
	string	IT1;
	string	narg1;
	string	narg2;
	string	narg3;
	string	narg4;
	string	narg5;
	string	narg6;
	string	narg7;
	
};


struct Domain {
	/********************* Domain for TEHL-model **************************/

	int			dkey;			/* domain type:								  */
								/*				1 - 1D (line contact)		  */	
								/*				2 - 2D (point contact)		  */
	int			rsx;
	int			rsy;
	int			Nx;				/* number of grid points in x-direction		  */
	int			Ny;				/* number of grid points in y-direction		  */
	int			Nz0;			/* number of grid points in z-direction		  */
	Drefine		RF;				/* refine data for subarea in domain		  */
	double		dx;				/* grid resolution in x-direction			  */
	double		dy;				/* grid resolution in y-direction			  */
	double		a;				/* a-semiaxis of contact ellipse, m			  */
	double		b;				/* b-semiaxis of contact ellipse, m			  */
	double		kx;				/* scale factor for area in x-direction		  */
	double		ky;				/* scale factor for area in y-direction		  */
	double		rRX;			/* RX for reduced geometry, m				  */
	double		rRY;			/* RY for reduced geometry, m				  */
	double		rYm;			/* Young modulus for material (reduced geom.) */
	double		W;				/* normal contact force, N (or N/m)			  */
	double		ph;				/* maximum of herzian pressure, Pa			  */
	double		delta;			/* maximum closure between bodies, m		  */
	double		h0;				/* inital gap between bodies, m				  */		
	double		rhoR;			/* ref. value for density, kg / m^3			  */
	double		muR;			/* ref. value for low-shear viscosity, Pa s   */
	double		etaR;			/* ref. value for generalized viscosity, Pa s */
	double		KR;				/* ref. value for th¡£conductivity, W/£¨m K£©  */
	double		CpR;			/* ref. value for heat capacity, J / (kg K)   */
	double		T0;				/* environment temperature, K				  */
	double		SRR[2];			/* slide-to-roll ratio						  */
	double		W_HD;			/* dimensionless load prameter				  */
	double		G_HD;			/* dimensionless geometry prameter			  */
	double		U_HD;			/* dimensionless velocity prameter			  */
	double		Hc;				/* central film thickness (empirical expr.)	  */
	double		Hm;				/* minimum film thickness (empirical expr.)	  */
	MM_density	RHO;
	MM_ls_viscosity	LSV;
	MM_gn_viscosity GNV;
	MM_th_conductivity THC;
	MM_ht_capacity HCP;
	string		pth;			/* path to the output folder				  */
	string		IT1;			/* additional inf. about domain (optional)	  */
	Point** Area;				/* domain points							  */
};


struct Body {
	/**** Geometry and physical parameters for body in dimension form ***/
	int			Num;				/* body number: (1 - bottom; 2 - top)		*/
	double		Rx;					/* radii of curvature (x-direction), m		*/
	double		Ry;					/* radii of curvature (y-direction), m		*/
	double		YM;					/* Young modulus, Pa						*/
	double		YMe;				/* effective YM, (YM/(1-PR*PR)), Pa			*/
	double		PR;					/* Poisson ratio							*/
	double		K;					/* thermal conductivity coeff., W/(m K)		*/
	double		Cp;					/* heat capacity, J/(kg K)					*/
	double		vel[2];				/* X,Y comp. of contact area velocity, m/s	*/
	double		Ts;					/* surface body temperature, K				*/
	double		Ft;					/* friction force on surface, N				*/
	double		fc;					/* friction coefficient						*/
	double		IT1;				/* additional item (for 2nd body - h0)		*/
	string		IT2;				/* additional inf. about body (optional)    */
};


struct LData {

	double** P;
	double** H;
	double H0;
};

struct MGData {

	LData* level;

};





void prsmodif(Domain& D);

void init_from_file(Domain& D, int kval);

void domain_initialize(int dkey, Body B1, Body B2, int rsx, int rsy, double a, double b, double kx, double ky, double W, double ph, double delta,
	double h0, double rhoR, double muR, double KR, double CpR, double T0, int pkey, string IT1, Domain& D);

void domain_finalize(Domain &D);

void elbody_initialize(int N, double Rx, double Ry, double YM, double PR, double K, double Cp,
	double vel[2], double Ts, double IT1, string IT2, Body& B);

void herzian_parameters(Body B1, Body B2, int kval, double W, double len, double& a, double& b, double& ph, double& delta);
/* void calc_body_geom(Body B, int i, int j, double dx, double dy, double a, double b, double kx, double ky, double h0); */

void set_density_model(MM_density& RHO, string fname);

void set_ls_viscosity_model(MM_ls_viscosity& LSV, string fname);

void set_gn_viscosity_model(MM_gn_viscosity& GNV, string fname);

void set_thermal_conductivity_model(MM_th_conductivity& THC, string fname);

void set_heat_capacity_model(MM_ht_capacity& HCP, string fname);

#endif