#pragma once
void ctof(int nc, double uc[], int nf, double uf[]);
void ftoc(int nf, double uf[], double rf[], int nc, double uc[],
	double rc[]);
void gauss_seidel(int n, double r[], double u[], double &dif_l1);
int i4_log_2(int i);
int i4_power(int i, int j);
void monogrid_poisson_1d(int n, double a, double b, double ua, double ub,
	double force(double x), double exact(double x), int *it_num,
	double u[]);
void monogrid_poisson_fem_1d(int n, double a, double b, double ua, double ub,
	double force(double x), double exact(double x), int *it_num,
	double u[]);
void multigrid_poisson_1d(int n, double a, double b, double ua, double ub,
	double force(double x), double exact(double x), int *it_num,
	double u[]);
double r8_max(double x, double y);
double *equidist_new(int n, double a, double b);