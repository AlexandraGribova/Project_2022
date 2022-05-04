#include<vector>
#include<iostream>
using namespace std;
class LOS_
{
private:
	vector<double> x0;
	int n, m;
	int maxiter = 1000;
	double nev = 0, eps = 1e-6;

	double DotProduct(vector<double> x, vector<double> y, int n)
	{
		double res = 0;
		for (int i = 0; i < n; i++)
		{
			res += x[i] * y[i];
		}
		return res;
	}

	void MatrixMult(vector<int> ia, vector<int> ja, int n, vector<double> al, vector<double> di, vector<double> x, vector<double>& b)
	{
		for (int i = 0; i < n; i++)
		{
			b[i] = x[i] * di[i];
			for (int k = ia[i] - 1; k < ia[i + 1] - 1; k++)
			{
				int j = ja[k] - 1;
				b[i] += al[k] * x[j];
				b[j] += al[k] * x[i];
			}
		}
	}

	int MSG(vector<int> ia, vector<int> ja, int n, vector<double> al, vector<double> di, vector<double>& x, vector<double> b, int maxiter, double eps)
	{
		double bnorm = sqrt(DotProduct(b, b, n));
		vector<double> r(n);
		vector<double> p(n);
		vector<double> q(n);


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

	/*int LOS(int* ia, int* ja, int n, double* al, double* di, double* x, double* b, int maxiter, double eps)
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
	}*/
public:
	LOS_(vector<int> _ig, vector<int> _jg, vector<double> _gg, vector<double> _diag, vector<double> _d, uint32_t n)
	{		
		x0.resize(n, 1);
		MSG(_ig, _jg, n, _gg, _diag, x0, _d, maxiter, eps);
	}

	vector<double> get_q()
	{
		return x0;
	}
	
};