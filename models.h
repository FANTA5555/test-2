#ifndef MODELS_H
#define MODELS_H

double mvec(double x, double y);

void check_LSS(double SRRx, double SRRy, double prs, double& tau);

double calc_integrals(dbvector f, double a, double dz, int n, int kval, int q);

void calc_gen_viscosity();

void calc_densviscNB(Domain& D, MM_density MMrho, MM_ls_viscosity MMlsvs);




void calc_physical_fields(Domain &D, MM_density MMrho, MM_ls_viscosity MMlsvs, MM_gn_viscosity MMgnvs,
	MM_th_conductivity MMthc, MM_ht_capacity MMhcp, double vel1[2], double vel2[2]);

void calc_fric_1D_old(Domain& D, Body &B1, Body &B2);
void calc_fric_1D(Domain& D, Body& B1, Body& B2);

void calc_fric_2D_old(Domain& D, Body& B1, Body& B2);

void calc_fric_2D(Domain& D, Body& B1, Body& B2);

#endif