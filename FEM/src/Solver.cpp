#include<vector>
#include<iostream>
#include <iterator>
using namespace std;
class LOS_
{
private:
	/*vector<double> X, di, al;
	vector < int32_t> b;
	vector<int> ia, ja;*/
	int n, m;
	int maxiter = 1000;
	double nev = 0, eps = 1e-6;
	double* x0 = new double[n] {1};

	double DotProduct(double* x, double* y, int n)
	{
		double res = 0;
		for (int i = 0; i < n; i++)
		{
			res += x[i] * y[i];
		}
		return res;
	}

	void MatrixMult(int* ia, int* ja, int n, double* al, double* di, double* x, double* b)
	{
		for (int i = 0; i < n; i++)
		{
			b[i] = x[i] * di[i];
			for (int k = ia[i]-1; k < ia[i + 1]-1; k++)
			{
				int j = ja[k]-1;
				b[i] += al[k] * x[j];
				b[j] += al[k] * x[i];
			}
		}
	}

	int MSG(int* ia, int* ja, int n, double* al, double* di, double* x, double* b, int maxiter, double eps)
	{
		double bnorm = sqrt(DotProduct(b, b, n));
		double* r = new double[n];
		double* p = new double[n];
		double* q = new double[n];
		MatrixMult(ia, ja, n, al, di, x, r);
		for (int i = 0; i < n; i++)
		{
			r[i] = b[i] - r[i];
			p[i] = r[i];
		}
		int k = 0;
		double alpha, betta, rnorm = sqrt(DotProduct(r, r, n));
		while (k<maxiter && rnorm / bnorm>eps)
		{
			MatrixMult(ia, ja, n, al, di, p, q);
			alpha = DotProduct(r, r, n) / DotProduct(q, p, n);
			betta = 1 / DotProduct(r, r, n);
			for (int i = 0; i < n; i++)
			{
				x[i] += alpha * p[i];
				r[i] -= alpha * q[i];
			}
			rnorm = sqrt(DotProduct(r, r, n));
			betta *= rnorm * rnorm;
			for (int i = 0; i < n; i++)
			{
				p[i] = r[i] + betta * p[i];
			}
			k++;
		}
		printf_s("relative residual is: %e\n", rnorm / bnorm);
		return k;
	}

	int LOS(int* ia, int* ja, int n, double* al, double* di, double* x, double* b, int maxiter, double eps)
	{
		double bnorm = sqrt(DotProduct(b, b, n));
		double* r = new double[n];
		double* p = new double[n];
		double* z = new double[n];
		double* Ar = new double[n];
		MatrixMult(ia, ja, n, al, di, x, r);
		for (int i = 0; i < n; i++)
		{
			r[i] = b[i] - r[i];
			z[i] = r[i];
		}
		MatrixMult(ia, ja, n, al, di, z, p);
		int k = 0;
		double alpha, betta, rnorm = sqrt(DotProduct(r, r, n));
		while (k<maxiter && rnorm / bnorm>eps)
		{
			alpha = DotProduct(p, r, n) / DotProduct(p, p, n);
			for (int i = 0; i < n; i++)
			{
				x[i] += alpha * z[i];
				r[i] -= alpha * p[i];
			}
			MatrixMult(ia, ja, n, al, di, r, Ar);
			betta = -DotProduct(p, Ar, n) / DotProduct(p, p, n);
			rnorm = sqrt(DotProduct(r, r, n));
			for (int i = 0; i < n; i++)
			{
				z[i] = r[i] + betta * z[i];
				p[i] = Ar[i] + betta * p[i];
			}
			k++;
		}
		printf_s("relative residual is: %e\n", rnorm / bnorm);
		return k;
	}
public:
	LOS_(vector<int> _ig, vector<int> _jg, vector<double> _gg, vector<double> _diag, vector<double> _d, uint32_t _N)
	{
		n = _N;
		int* ia = &_ig[0];
		int* ja = &_jg[0];
		double* di = &_diag[0];
		double* al = &_gg[0];
		double* x = new double[n];
		
		double* b = &_d[0];
		int m = _ig[_ig.size()- 1] - 1;
		MSG(ia, ja, n, al, di, x0, b, maxiter, eps);
	}

	vector<double> get_q()
	{
		vector<double> output(x0, x0+sizeof(x0));
		return output;
	}

};