#ifndef SOLVERS_H
#define SOLVERS_H

void dumpvar1D(Domain D, int kval, int dtype, int kout);

void dumpvar2D(Domain D, int kval, int src, int kout);

void calc_triplets(int n, int m, double* A, int*& Ti, int*& Tj, double*& Tx);

void calc_rhoeps(int i, int j, Domain D, Body B1, Body B2, double eps_curr, double rhox_curr, double rhoy_curr);


void calc_coef1D(int i, Domain D, Body B1, Body B2, double* eps, double* rhox, double* H,
	double& cA, double& cB, double& cE, double& cF, double& cJ, int kval);

void calc_coef2D(int i, int j, Domain D, Body B1, Body B2, double** eps, double** rhox, double** rhoy, double** H,
	double& cA, double& cB, double& cC, double& cD, double& cE, double& cF, int kval);

double dFdP(int k, int l, int i, int j, double** rhox, double** rhoy, double dx, double dy, double theta);


void calc_pressure_field_1D_old(Domain& D, Body B1, Body B2);
void calc_pressure_field_1D(Domain& D, Body B1, Body B2);
void calc_pressure_field_1D_backup(Domain& D, Body B1, Body B2);

void calc_pressure_field_2D_old(Domain& D, Body B1, Body B2);

void calc_pressure_field_2D_LNRLX(Domain& D, Body B1, Body B2);

void calc_pressure_field_2D(Domain& D, Body B1, Body B2);

void calc_pressure_field_2D_FMG(Domain& D, Body B1, Body B2, int kmax);

void prnpress(Domain D);


#endif
