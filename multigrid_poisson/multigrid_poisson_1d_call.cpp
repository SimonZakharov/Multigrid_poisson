#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <fstream>

#include "multigrid_poisson_1d.hpp"

void test01_mono();
void test01_multi();
void test02_mono();
void test02_multi();
double exact1(double x);
double force1(double x);
double exact2(double x);
double force2(double x);

using namespace std;

ofstream output;

int main(int argc, char *argv[])
{
	output.open("output.txt");
	output << endl;
	output << "MULTIGRID_POISSON_1D_PRB:\n";

	test01_mono();
	test01_multi();
	//test02_mono();
	//test02_multi();
	
	output << "\n";
	output << "MULTIGRID_POISSON_1D_PRB:\n";
	output << "  Normal end of execution.\n";
	output << "\n";

	output.close();
	//system("pause");
	return 0;
}

void test01_mono()
{
	double a;
	double b;
	double difmax;
	int i;
	int it_num;
	int k;
	int n;
	double *u;
	double ua;
	double ub;
	double *x;

	if (!output.is_open())
	{
		printf("MULTIGRID_POISSON_1D - Fatal error!\n");
		printf("Failed to open the output file.\n");
		exit(1);
	}

	output << "\n";
	output << "TEST01_MONO\n";
	output << "  MONOGRID_POISSON_1D solves a 1D Poisson BVP\n";
	output << "  using the Gauss-Seidel method.\n";

	a = 0.0;
	b = 1.0;
	ua = 0.0;
	ub = 0.0;

	output << "\n";
	output << "  -u''(x) = 1, for 0 < x < 1\n";
	output << "  u(0) = u(1) = 0.\n";
	output << "  Solution is u(x) = ( -x^2 + x ) / 2\n";

	for (k = 5; k <= 5; k++)
	{
		n = i4_power(2, k);

		u = (double *)malloc((n + 1) * sizeof(double));
		x = equidist_new(n + 1, a, b);

		output << "\n";
		output << "  Mesh index K = " << k << endl;
		output << "  Number of intervals N=2^K = " << n << endl;
		output << "  Number of nodes = 2^K+1 = " << n + 1 << endl;

		monogrid_poisson_1d(n, a, b, ua, ub, force1, exact1, &it_num, u);

		output << "\n";
		output << "     I        X(I)      U(I)         U Exact(X(I))\n";
		output << "\n";
		for (i = 0; i < n + 1; i++)
		{
			output << i << "\t\t" << x[i] << "\t\t" << u[i] << "\t\t" << exact1(x[i]) << endl;
		}

		output << "\n";

		difmax = 0.0;
		for (i = 0; i < n + 1; i++)
		{
			difmax = r8_max(difmax, fabs(u[i] - exact1(x[i])));
		}
		output << "  Maximum error = " << difmax << endl;
		output << "  Number of iterations = " << it_num << endl;

		free(u);
		free(x);
	}
	
	n = i4_power(2, k);

	u = (double *)malloc((n + 1) * sizeof(double));
	x = equidist_new(n + 1, a, b);

	output << "\n";
	output << "TEST01_MONO\n";
	output << "  MONOGRID_POISSON_1D solves a 1D Poisson BVP\n";
	output << "  using the finite elements method.\n";

	monogrid_poisson_fem_1d(n, a, b, ua, ub, force1, exact1, &it_num, u);

	output << "\n";
	output << "  Mesh index K = " << k << endl;
	output << "  Number of intervals N=2^K = " << n << endl;
	output << "  Number of nodes = 2^K+1 = " << n + 1 << endl;
	
	output << "\n";
	output << "     I        X(I)      U(I)         U Exact(X(I))         Error\n";
	output << "\n";
	for (i = 0; i < n + 1; i++)
	{
		output <<  i << "\t\t" << x[i] << "\t\t" << u[i] << "\t\t" << exact1(x[i]) << "\t\t" << abs(u[i] - exact1(x[i])) << endl;
	}

	output << "\n";

	difmax = 0.0;
	for (i = 0; i < n + 1; i++)
	{
		difmax = r8_max(difmax, fabs(u[i] - exact1(x[i])));
	}
	output << "  Maximum error = " << difmax << endl;
	output << "  Number of iterations = " << it_num << endl;

	free(u);
	free(x);

	return;
}

void test01_multi()
{
	double a;
	double b;
	double difmax;
	int i;
	int it_num;
	int k;
	int n;
	double *u;
	double ua;
	double ub;
	double *x;

	if (!output.is_open())
	{
		output << "MULTIGRID_POISSON_1D - Fatal error!\n";
		output << "Failed to open the output file.\n";
		exit(1);
	}
	output << "\n";
	output << "TEST01_MULTI\n";
	output << "  MULTIGRID_POISSON_1D solves a 1D Poisson BVP\n";
	output << "  using the multigrid method.\n";

	a = 0.0;
	b = 1.0;
	ua = 0.0;
	ub = 0.0;

	output << "\n";
	output << "  -u''(x) = 1, for 0 < x < 1\n";
	output << "  u(0) = u(1) = 0.\n";
	output << "  Solution is u(x) = ( -x^2 + x ) / 2\n";

	for (k = 5; k <= 5; k++)
	{
		n = i4_power(2, k);

		u = (double *)malloc((n + 1) * sizeof(double));
		x = equidist_new(n + 1, a, b);

		output << "\n";
		output << "  Mesh index K = " << k << endl;
		output << "  Number of intervals N=2^K = " << n << endl;
		output << "  Number of nodes = 2^K+1 = " << n + 1 << endl;

		multigrid_poisson_1d(n, a, b, ua, ub, force1, exact1, &it_num, u);

		output << "\n";
		output << "     I        X(I)      U(I)         U Exact(X(I))\n";
		output << "\n";
		for (i = 0; i < n + 1; i++)
		{
			output << i << "\t\t" << x[i] << "\t\t" << u[i] << "\t\t" << exact1(x[i]) << "\t\t" << abs(u[i] - exact1(x[i])) << endl;
		}

		output << "\n";

		difmax = 0.0;
		for (i = 0; i < n + 1; i++)
		{
			difmax = r8_max(difmax, fabs(u[i] - exact1(x[i])));
		}
		output << "  Maximum error = " << difmax << endl;
		output << "  Number of iterations = " << it_num << endl;

		free(u);
		free(x);
	}
	return;
}

double exact1(double x)
{
	double value;
	value = 0.5 * (-x * x + x);
	return value;
}

double force1(double x)
{
	double value;
	value = 1.0;
	return value;
}

void test02_mono()
{
	double a;
	double b;
	double difmax;
	int i;
	int it_num;
	int k;
	int n;
	double *u;
	double ua;
	double ub;
	double *x;

	printf("\n");
	printf("TEST02_MONO\n");
	printf("  MONOGRID_POISSON_1D solves a 1D Poisson BVP\n");
	printf("  using the Gauss-Seidel method.\n");

	a = 0.0;
	b = 1.0;
	ua = 0.0;
	ub = 0.0;

	printf("\n");
	printf("  -u''(x) = - x * (x+3) * exp(x), for 0 < x < 1\n");
	printf("  u(0) = u(1) = 0.\n");
	printf("  Solution is u(x) = x * (x-1) * exp(x)\n");

	for (k = 5; k <= 5; k++)
	{
		n = i4_power(2, k);

		u = (double *)malloc((n + 1) * sizeof(double));
		x = equidist_new(n + 1, a, b);

		printf("\n");
		printf("  Mesh index K = %d\n", k);
		printf("  Number of intervals N=2^K = %d\n", n);
		printf("  Number of nodes = 2^K+1 =   %d\n", n + 1);

		monogrid_poisson_1d(n, a, b, ua, ub, force2, exact2, &it_num, u);

		printf("\n");
		printf("     I        X(I)      U(I)         U Exact(X(I))\n");
		printf("\n");
		for (i = 0; i < n + 1; i++)
		{
			printf("  %4d  %10f  %14g  %14g\n", i, x[i], u[i], exact2(x[i]));
		}

		printf("\n");

		difmax = 0.0;
		for (i = 0; i < n + 1; i++)
		{
			difmax = r8_max(difmax, fabs(u[i] - exact2(x[i])));
		}
		printf("  Maximum error = %g\n", difmax);
		printf("  Number of iterations = %d\n", it_num);

		free(u);
		free(x);
	}
	return;
}

void test02_multi()
{
	double a;
	double b;
	double difmax;
	int i;
	int it_num;
	int k;
	int n;
	double *u;
	double ua;
	double ub;
	double *x;

	printf("\n");
	printf("TEST02_MULTI\n");
	printf("  MULTIGRID_POISSON_1D solves a 1D Poisson BVP\n");
	printf("  using the multigrid method.\n");

	a = 0.0;
	b = 1.0;
	ua = 0.0;
	ub = 0.0;

	printf("\n");
	printf("  -u''(x) = - x * (x+3) * exp(x), for 0 < x < 1\n");
	printf("  u(0) = u(1) = 0.\n");
	printf("  Solution is u(x) = x * (x-1) * exp(x)\n");

	for (k = 5; k <= 5; k++)
	{
		n = i4_power(2, k);

		u = (double *)malloc((n + 1) * sizeof(double));
		x = equidist_new(n + 1, a, b);

		printf("\n");
		printf("  Mesh index K = %d\n", k);
		printf("  Number of intervals N=2^K = %d\n", n);
		printf("  Number of nodes = 2^K+1 =   %d\n", n + 1);

		multigrid_poisson_1d(n, a, b, ua, ub, force2, exact2, &it_num, u);

		printf("\n");
		printf("     I        X(I)      U(I)         U Exact(X(I))\n");
		printf("\n");
		for (i = 0; i < n + 1; i++)
		{
			printf("  %4d  %10f  %14g  %14g\n", i, x[i], u[i], exact2(x[i]));
		}

		printf("\n");

		difmax = 0.0;
		for (i = 0; i < n + 1; i++)
		{
			difmax = r8_max(difmax, fabs(u[i] - exact2(x[i])));
		}
		printf("  Maximum error = %g\n", difmax);
		printf("  Number of iterations = %d\n", it_num);

		free(u);
		free(x);
	}
	return;
}

double exact2(double x)
{
	double value;
	value = x * (x - 1.0) * exp(x);
	return value;
}

double force2(double x)
{
	double value;
	value = -x * (x + 3.0) * exp(x);
	return value;
}