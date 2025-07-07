#ifndef DEFORM_H
#define DEFORM_H

double IC_F(double x, double y);

double calc_IC_1D(double ki, double delta);

double calc_IC_2D(double ki, double lj, double delta);

double calc_roughness();

void calc_deformDS(Domain& D, Body B1, Body B2);

void calc_deform(Domain &D, Body B1, Body B2);


#endif