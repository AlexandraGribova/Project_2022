#include<vector>
#include<iostream>
using namespace std;
class LOS_solver
{
private:
	vector<double> X, d, di, gg, pr;
	vector<int> ig, jg;
	vector<double> r, z, p, boof, boof1, L, U;
	int N;
	int maxiter = 100000;
	double nev = 0, eps = 1e-6;
	void r0_s()
	{
		int i, j, k;
		int i0, i1;
		for (i = 0; i < N; i++)
		{
			i0 = ig[i]-1;
			i1 = ig[i + 1]-1;
			r[i] = pr[i] - d[i] * X[i];//f-Ax0
			for (k = i0; k < i1; k++)
			{
				j = jg[k];
				r[i] -= gg[k] * X[j];
				r[j] -= gg[k] * X[i];
			}
		}
	}
	void L_1(vector <double> pr, vector <double>& r)
	{
		double s;
		int	i0, i1, k;
		for (int i = 0; i < N; i++) {
			s = pr[i];
			i0 = ig[i]-1;
			i1 = ig[i + 1]-1;
			for (k = i0; k < i1; k++)
				s -= L[k] * r[jg[k]];
			r[i] = s / di[i];
		}
	}
	void U_1(vector <double> pr, vector <double>& z)
	{
		int	i0, i1, j, k;
		z = pr;
		for (int i = N - 1; i >= 0; i--)
		{
			i0 = ig[i]-1;
			i1 = ig[i + 1]-1;
			z[i] = z[i] / di[i];
			for (k = i0; k < i1; k++)
			{
				j = jg[k];
				z[j] -= U[k] * z[i];
			}
		}
	}
	void p0()
	{
		int i, k, j, i0, i1;
		for (i = 0; i < N; i++)
		{
			p[i] = d[i] * z[i];
			i0 = ig[i]-1;
			i1 = ig[i + 1]-1;
			for (k = i0; k < i1; k++)
			{
				j = jg[k];
				p[i] += gg[k] * z[j];
				p[j] += gg[k] * z[i];
			}
		}
	}
	double scalar_mult(vector <double> v1, vector <double> v2, int size)
	{
		double s = 0;
		for (int i = 0; i < size; i++)
			s += v1[i] * v2[i];
		return s;
	}

	void X_k(double a, vector <double> z)
	{
		for (int i = 0; i < N; i++)
			X[i] += a * z[i];
	}
	void R_k(double a, vector <double> p)
	{
		for (int i = 0; i < N; i++)
			r[i] -= a * p[i];
	}
	double Norm(vector<double> X)
	{
		double norma = 0;
		for (int i = 0; i < N; i++)
		{
			norma += X[i] * X[i];
		}
		norma = sqrt(norma);
		return norma;
	}

	void AVec(vector <double>& x, vector <double>& y)
	{
		for (int i = 0; i < N; i++) {
			y[i] = d[i] * x[i];
			int i0 = ig[i]-1;
			int i1 = ig[i + 1]-1;
			for (int k = i0; k < i1; k++) {
				int j = jg[k];
				y[i] += gg[k] * x[j];
				y[j] += gg[k] * x[i];
			}
		}
	}
	void Z_k(vector <double> dat, double b, vector <double> d, vector <double>& z)
	{
		for (int i = 0; i < N; i++)
			z[i] = dat[i] + b * d[i];
	}
	void vec_DI(vector <double>vec, vector <double>& res)
	{
		for (int i = 0; i < N; i++)
			res[i] = vec[i] / d[i];
	}
	double _nev()
	{

		int i, j, k;
		int i0, i1;
		for (i = 0; i < N; i++)
		{
			i0 = ig[i]-1;
			i1 = ig[i + 1]-1;
			boof1[i] = pr[i] - d[i] * X[i];//f-Ax0
			for (k = i0; k < i1; k++)
			{
				j = jg[k];
				boof1[i] -= gg[k] * X[j];
				boof1[j] -= gg[k] * X[i];
			}
		}
		return(Norm(boof1) / Norm(pr));
	}
public:
	LOS_solver(vector<int> _ig, vector<int> _jg, vector<double> _gg, vector<double> _diag, uint32_t _N)
	{
		ig = _ig;
		jg = _jg;
		gg = _gg;
		N = _N;
		d = _diag;
		X.resize(N);
		for (int k = 0; k < N; k++)X[k] = 1;
		pr.resize(N);
		di.resize(N);
		p.resize(N);
		r.resize(N);
		z.resize(N);
		boof.resize(N);
		boof1.resize(N);
		int size = ig[N];
		L.resize(ig[N]);
		U.resize(ig[N]);
		LUS_factorisation();
		LOS_DI();
		for (int i = 0; i < N; i++) {
			cout.setf(ios::scientific);
			cout << X[i] << endl;
		}

	}
	void LOS_LU()
	{
		double a, b, nr, nf;
		r0_s();//r=f-A*x
		L_1(r, r);//r=(L^-1)*r
		U_1(r, z);//z=(U^-1)*r
		p0();//p=Az
		L_1(p, p);//p=(L^-1)*p
		nr = Norm(r);
		nf = Norm(pr);
		for (int i = 0; i < maxiter; i++)
		{
			if (nr <= eps * nf)
				break;
			a = scalar_mult(p, r, N) / scalar_mult(p, p, N);//a=(p,r)/(p,p)  (3.35)
			X_k(a, z);//x=x(k-1)+a*z  (3.36)
			R_k(a, p);//r=r(k-1)+p*a  (3.37)
			U_1(r, boof);//boof=(U^-1)*r
			AVec(boof, boof1);//boof1=A*boof Все сделать по образцу этой функции
			L_1(boof1, boof);//boof=(L^-1)*boof1
			b = -scalar_mult(p, boof, N) / scalar_mult(p, p, N);//b=-(p,boof)/(p,p)  (3.38)
			U_1(r, boof1);//boof1=(U^-1)*r
			Z_k(boof1, b, z, z);//z=boof1+b*z(k-1)
			Z_k(boof, b, p, p);//p=длинное выражение+b*p  (3.40)
			nr = Norm(r);
			cout.setf(ios::scientific);
			cout << i << "  Условие малости относительной невязки:" << nr / nf << endl;
		}
		nev = _nev();//nev=|f-Ax|/|f|
		cout << "  Относительная невязка:" << nev << endl;
	}

	void LOS_DI()
	{
		double a, b, nr, nf;
		r0_s();//r=f-A*x
		vec_DI(r, r);//r=r/di
		z = r;
		p0();//p=Az
		vec_DI(p, p);//p=p/di
		nr = Norm(r);
		nf = Norm(pr);
		for (int i = 0; i < maxiter; i++)
		{
			if (nr <= eps * nf)
				break;
			a = scalar_mult(p, r, N) / scalar_mult(p, p, N);//a=(p,r)/(p,p)  (3.35)
			X_k(a, z);//x=x(k-1)+a*z  (3.36)
			R_k(a, p);//r=r(k-1)+p*a  (3.37)
			AVec(r, boof);//boof=A*r 
			vec_DI(boof, boof);
			b = -scalar_mult(p, boof, N) / scalar_mult(p, p, N);//b=-(p,boof)/(p,p)  (3.38)
			Z_k(r, b, z, z);//z=r+b*z(k-1)
			Z_k(boof, b, p, p);//p=boof+b*p  (3.40)
			nr = Norm(r);
			cout.setf(ios::scientific);
			cout << i << "  Условие малости относительной невязки:" << nr / nf << endl;
		}
		nev = _nev();//nev=|f-Ax|/|f|
		cout << "  Относительная невязка:" << nev << endl;
	}

	void LOS()
	{
		double a, b, nr, nf;
		r0_s();//r=f-A*x
		z = r;//z0=r0
		AVec(z, p);//p0=A*z0
		nr = Norm(r);
		nf = Norm(pr);
		for (int i = 0; i < maxiter; i++)
		{
			if (nr <= eps * nf)
				break;
			a = scalar_mult(p, r, N) / scalar_mult(p, p, N);//a=(p,r)/(p,p)  (3.35)
			X_k(a, z);//x=x(k-1)+a*z  (3.36)
			R_k(a, p);//r=r(k-1)+p*a  (3.37)
			AVec(r, boof);//boof=A*r 
			b = -scalar_mult(p, boof, N) / scalar_mult(p, p, N);//b=-(p,boof)/(p,p)  (3.38)
			Z_k(r, b, z, z);//z=r+b*z(k-1)
			Z_k(boof, b, p, p);//p=boof+b*p  (3.40)
			nr = Norm(r);
			cout.setf(ios::scientific);
			if (i % 1000 == 0)
				cout << i << " Условие малости относительной невязки:" << nr / nf << endl;

		}
		nev = _nev();//nev=|f-Ax|/|f|
		cout << "  Относительная невязка:" << nev << endl;
	}

	void LUS_factorisation()
	{
		int i, j, k;
		int i0, i1, ki, kj;
		double suml, sumu, sumd;

		for (i = 0; i < N; i++)
		{
			i0 = ig[i]-1;
			i1 = ig[i + 1]-1;
			sumd = 0;
			for (k = i0; k < i1; k++)
			{
				j = jg[k];
				ki = i0;
				kj = ig[j]-1;//что такое kj - соответсвующий номер элемента для домножения
				suml = sumu = 0;
				while (ki < k)
					if (jg[kj] == jg[ki]) {
						sumu += U[ki] * L[kj];
						suml += U[kj] * L[ki];
						kj++;
						ki++;
					}
					else
						if (jg[kj] < jg[ki])
							kj++;
						else
							ki++;
				U[k] = (gg[k] - sumu) / di[j];
				L[k] = (gg[k] - suml) / di[j];
				sumd += U[k] * L[k];
			}
			di[i] = sqrt(d[i] - sumd);
		}

	}
};