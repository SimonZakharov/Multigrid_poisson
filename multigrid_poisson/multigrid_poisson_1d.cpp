# include <cstdlib>
# include <iostream>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

//	����������� ���������� - ����� ��� ���������� ������ ����� multigrid_poisson_1d_call.cpp
#include "multigrid_poisson_1d.hpp"

//	������� ������ �� ������ �����
//	from Coarser TO Finer grid 
void ctof(int nc, double uc[], int nf, double uf[])
{
	int ic;
	int iff;

	for (ic = 0; ic < nc; ic++)
	{
		iff = 2 * ic;
		uf[iff] = uf[iff] + uc[ic];
	}

	for (ic = 0; ic < nc - 1; ic++)
	{
		iff = 2 * ic + 1;
		uf[iff] = uf[iff] + 0.5 * (uc[ic] + uc[ic + 1]);
	}

	return;
}

//	�������� ������� � ������ �� ������ �����
//	from Finer TO Coarser grid
void ftoc(int nf, double uf[], double rf[], int nc, double uc[], double rc[])
{
	int ic;
	int iff;

	for (ic = 0; ic < nc; ic++)
	{
		uc[ic] = 0.0;
	}

	rc[0] = 0.0;
	for (ic = 1; ic < nc - 1; ic++)
	{
		iff = 2 * ic;
		rc[ic] = 4.0 * (rf[iff] + uf[iff - 1] - 2.0 * uf[iff] + uf[iff + 1]);
	}
	rc[nc - 1] = 0.0;

	return;
}

//	��� ���� ��� ������������� ��������
//	����� ������-�������
void gauss_seidel(int n, double r[], double u[], double &dif_l1)
{
	int i;
	double u_old;

	dif_l1 = 0.0;

	for (i = 1; i < n - 1; i++)
	{
		u_old = u[i];
		u[i] = 0.5 * (u[i - 1] + u[i + 1] + r[i]);
		dif_l1 = dif_l1 + fabs(u[i] - u_old);
	}

	return;
}

//	����� ��������
//	����� ��� ������� ���������������� �������, ������� ���������� � ���
//	��� ������������� �������� �������� �������
void tma(int n, double *diag, double *up, double *down, double *f, double *u)
{
	double *d, *l;
	d = new double[n]; l = new double[n];
	//	������ ��������
	d[0] = -up[0] / diag[0];
	l[0] = f[0] / diag[0];
	for (int i = 1; i < n; i++)
	{
		if (i != n - 1)
			d[i] = -up[i] / (diag[i] + down[i - 1] * d[i - 1]);
		else d[i] = 0;
		l[i] = (f[i] - down[i - 1] * l[i - 1]) / (diag[i] + down[i - 1] * d[i - 1]);
	}
	//	�������� ��������
	u[n - 1] = l[n - 1];
	for (int i = n - 2; i >= 0; i--)
	{
		u[i] = d[i] * u[i + 1] + l[i];
	}
	delete[]d;
	delete[]l;
}

//	��� ������� ���������� ����� ����� ��������� log_2(i)
//	������ ��� i4_log_2(i) + 1 - ��� ���������� �������� ���� � ������ ����� i
//	��������:
//        I  I4_LOG_10
//    -----  --------
//        0    0
//        1    0
//        2    1
//        3    1
//        4    2
//        5    2
//        7    2
//        8    3
//        9    3
//     1000    9
//     1024   10
int i4_log_2(int i)
{
	int i_abs;
	int two_pow;
	int value;

	if (i == 0)
	{
		value = 0;
	}
	else
	{
		value = 0;
		two_pow = 2;

		i_abs = abs(i);

		while (two_pow <= i_abs)
		{
			value = value + 1;
			two_pow = two_pow * 2;
		}
	}

	return value;
}


//	���������� �������� i � ������� j
int i4_power(int i, int j)
{
	int k;
	int value;

	if (j < 0)
	{
		if (i == 1)
		{
			value = 1;
		}
		else if (i == 0)
		{
			cerr << "\n";
			cerr << "I4_POWER - Fatal error!\n";
			cerr << "  I^J requested, with I = 0 and J negative.\n";
			exit(1);
		}
		else
		{
			value = 0;
		}
	}
	else if (j == 0)
	{
		if (i == 0)
		{
			cerr << "\n";
			cerr << "I4_POWER - Fatal error!\n";
			cerr << "  I^J requested, with I = 0 and J = 0.\n";
			exit(1);
		}
		else
		{
			value = 1;
		}
	}
	else if (j == 1)
	{
		value = i;
	}
	else
	{
		value = 1;
		for (k = 1; k <= j; k++)
		{
			value = value * i;
		}
	}
	return value;
}

//	����� ���������� ������� ���������� ������ �������� ������� �������� ���������
//	� �������� �������� ���������� ����� ������-�������
//	�������, ������� �������, ��� ��������� � ������������� �������
//	��� �������� �� ���� (���������):
//	n - ���������� �������������
//	a, b - ������ � ����� �������
//	ua, ub - ��������� ������� ������� ���� �� ������ �������
//	force - ��� �������, ������� � ������ �����
//	exact - ��� �������, ������������ ������ �������
//	��� ��������� (���������):
//	it_num - ���������� ��������
//	u - �������
void monogrid_poisson_1d(int n, double a, double b, double ua, double ub, double force(double x), double exact(double x), int *it_num, double u[])
{
	double d1;
	double h;
	int i;
	double *r;
	double tol;
	double *x;
	//
	//  �������������
	//
	tol = 0.0001;
	//
	//  ����������� ����
	//
	x = equidist_new(n + 1, a, b);
	//
	//  ��������� ������ �����
	//
	r = new double[n + 1];

	r[0] = ua;
	h = (b - a) / (double)(n);
	for (i = 1; i < n; i++)
	{
		r[i] = h * h * force(x[i]);
	}
	r[n] = ub;

	for (i = 0; i <= n; i++)
	{
		u[i] = 0.0;
	}
	*it_num = 0;
	//
	//  �������� ������ ������-�������
	//
	for (; ; )
	{
		*it_num = *it_num + 1;

		gauss_seidel(n + 1, r, u, d1);

		if (d1 <= tol)
		{
			break;
		}
	}

	delete[] r;
	delete[] x;

	return;
}

//	����� ���������� ������� ���������� ������ �������� ������� �������� ���������
//	� �������� �������� ������������ ����� ��������
//	�������, ������� �������, ��� ��������� � ������������� �������
//	��� �������� �� ���� (���������):
//	n - ���������� �������������
//	a, b - ������ � ����� �������
//	ua, ub - ��������� ������� ������� ���� �� ������ �������
//	force - ��� �������, ������� � ������ �����
//	exact - ��� �������, ������������ ������ �������
//	��� ��������� (���������):
//	it_num - ���������� ��������
//	u - �������
void monogrid_poisson_fem_1d(int n, double a, double b, double ua, double ub, double force(double x), double exact(double x), int *it_num, double u[])
{
	int i;
	double *x;
	double h;

	double g[2][2]; //	���������� ������� ���������
	double el_r[2];	//	���������� ������ ������ �����
	//	������� �������
	double *diag = new double[n];	for (int i = 0; i < n; i++) diag[i] = 0;
	double *up = new double[n];		for (int i = 0; i < n; i++) up[i] = 0;
	double *down = new double[n];	for (int i = 0; i < n; i++) down[i] = 0;
	//	���������� ������ ������ �����
	double *r = new double[n + 1];		fill(&r[0], &r[n + 1], 0.0);
	//	����������� �����
	x = equidist_new(n + 1, a, b);
	h = (b - a) / (double)n;
	//	���� ������ �� ���������
	for (i = 0; i < n; i++)
	{
		//	��������� ���������� ������� ���������
		g[0][0] = g[1][1] = 1.0 / h;
		g[0][1] = g[1][0] = -1.0 / h;
		//	��������� ���������� ������ ������ �����
		el_r[0] = h * (2 * force(x[i]) + force(x[i + 1])) / 6;
		el_r[1] = h * (force(x[i]) + 2 * force(x[i + 1])) / 6;
		//	�������� ���������� ������� �������
		diag[i] += g[0][0];
		diag[i + 1] += g[1][1];
		up[i] += g[0][1];
		down[i] += g[1][0];
		//	�������� ���������� ������ ������ �����
		r[i] += el_r[0]; r[i + 1] += el_r[1];
	}
	//	����� ��������� �������
	diag[0] = 1; up[0] = 0;
	diag[n] = 1; up[n] = 0; down[n - 1] = 0;
	r[0] = ua; r[n] = ub;
	tma(n + 1, diag, up, down, r, u);
	return;
}

//	������� � �������� ���� ������
//	����� ���������� ������� ���������� ������ ��������, ����������� � �������������� �������������� ������
//	���������� ������:
//	-U''(X) = F(X), A < X < B
//	��������� ������� ������� ����: U(A) = UA; U(B) = UB
//	��� �������� �� ���� (���������):
//	n - ����� �������������� (������ ���� �������� 2)
//	a, b - ����� ������ � ����� �������
//	ua, ub - ��������� ������� �� ������ (������� ����)
//	f - ��� �������, ������� � ������ ����� ���������
//	exact - ��� �������, ����������� ������ �������
//	��� ��������� (���������):
//	it_num - ���������� ��������
//	u - �������
void multigrid_poisson_1d(int n, double a, double b, double ua, double ub, double f(double x), double exact(double x), int *it_num, double u[])
{
	double d0;
	double d1;
	double h;
	int i;
	int it;
	int j;
	int k;
	int l;
	int ll;
	int m;
	int nl;
	double *r;
	double s;
	double tol;
	double utol;
	double *uu;
	double *x;
	//
	//  ���������, ���������� �� ������ (�������� �� N �������� 2)
	//
	k = i4_log_2(n);

	if (n != i4_power(2, k))
	{
		cout << "\n";
		cout << "MULTIGRID_POISSON_1D - Fatal error!\n";
		cout << "  N is not a power of 2.\n";
		exit(1);
	}

	nl = n + n + k - 2;
	//
	//  �������������
	//
	it = 4;
	*it_num = 0;
	tol = 0.0001;
	utol = 0.7;
	m = n;
	//
	//  ����������� �����.
	//
	x = equidist_new(n + 1, a, b);
	//
	//  ��������� ������ �����.
	//
	r = new double[nl];
	r[0] = ua;
	h = (b - a) / (double)(n);
	for (i = 1; i < n; i++)
	{
		r[i] = h * h * f(x[i]);
	}
	r[n] = ub;

	uu = new double[nl];

	for (i = 0; i < nl; i++)
	{
		uu[i] = 0.0;
	}
	//
	//  L ����� �� ������� ���������
	//  LL ����� �� ������������� ���������
	//
	l = 0;
	ll = n - 1;
	//
	//  ����� �������� ������ ������-�������
	//
	d1 = 0.0;
	j = 0;

	for (; ; )
	{
		d0 = d1;
		j = j + 1;
		gauss_seidel(n + 1, r + l, uu + l, d1);
		*it_num = *it_num + 1;
		//
		//  �� ������ ������ ���������� �� ����� 4 ��������
		//
		if (j < it)
		{
			continue;
		}
		//
		//	�������� ����������, ����������� ������������������, ��������� �� ����� ������ ����� - ����� �� �����
		//
		else if (d1 < tol && n == m)
		{
			break;
		}
		//
		//	�������� ����������; ���������� ������������������; ��������� �� ������ �����
		//
		else if (d1 < tol)
		{
			ctof(n + 1, uu + l, n + n + 1, uu + (l - 1 - n - n));

			n = n + n;
			ll = l - 2;
			l = l - 1 - n;
			j = 0;
		}
		//
		//	�������� ����������; ���������� ���������; �������� ������� �� ������ �����
		//
		else if (utol * d0 <= d1 && 2 < n)
		{
			ftoc(n + 1, uu + l, r + l, (n / 2) + 1, uu + (l + n + 1), r + (l + n + 1));

			n = n / 2;
			l = ll + 2;
			ll = ll + n + 1;
			j = 0;
		}
	}

	for (i = 0; i < n + 1; i++)
	{
		u[i] = uu[i];
	}
	delete[] r;
	delete[] uu;
	delete[] x;

	return;
}

//	����� ���������� ���������� �� ���� ������������ ����� ���� R8 (double)
double r8_max(double x, double y)
{
	double value;

	if (y < x)
	{
		value = x;
	}
	else
	{
		value = y;
	}
	return value;
}

//	�������� �������:
//	������� �� a_first �� a_last ����������� n ��������������� ������� �� n - 1 ��������
//	����� ��� n ����� ����� ������������ ��� ���� ������ �����
double *equidist_new(int n, double a_first, double a_last)
{
	double *a;
	int i;

	a = new double[n];

	if (n == 1)
	{
		a[0] = (a_first + a_last) / 2.0;
	}
	else
	{
		for (i = 0; i < n; i++)
		{
			a[i] = ((double)(n - 1 - i) * a_first
				+ (double)(i)* a_last)
				/ (double)(n - 1);
		}
	}
	return a;
}