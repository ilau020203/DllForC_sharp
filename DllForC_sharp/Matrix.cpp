#include "pch.h"
#include <ctime>
#include <stdexcept>
#include "Matrix.h"

using namespace std;

int main() { return 0; }
Matrix::Matrix(int n)
{
	if (n < 1)
		throw invalid_argument("Order must be an integer greater than 1.");
	N = n;
	cols = new double[n];
	rows = new double[n];
	cols[0] = 150 + 2 * n;
	rows[0] = cols[0];
	for (int i = 1; i < N; i++)
	{
		cols[i] = 3 * n * i * 17 % 100 - 50;
		rows[i] = 4 * n * i * 21 % 100 - 50;
	}
}
Matrix::Matrix(int n, double* cols_, double* rows_)
{
	N = n;
	if (n < 1)
		throw invalid_argument("Order must be an integer greater than 1.");
	cols = new double[n];
	rows = new double[n];
	if (cols[0] != rows[0])
		throw invalid_argument("First element must be equal.");
	if (cols[0] == 0)
		throw invalid_argument("First element must not be 0.");
	for (int i = 0; i < N; i++)
	{
		cols[i] = cols_[i];
		rows[i] = rows_[i];
	}
}
Matrix::Matrix(const Matrix& other)
{
	N = other.N;
	cols = new double[N];
	rows = new double[N];
	for (int i = 0; i < N; i++)
	{
		cols[i] = other.cols[i];
		rows[i] = other.rows[i];
	}
}
Matrix::~Matrix()
{
	delete[] rows;
	delete[] cols;
}
Matrix& Matrix::operator=(const Matrix& other)
{
	if (this == &other)
		return *this;
	delete[] cols;
	delete[] rows;
	N = other.N;
	cols = new double[N];
	rows = new double[N];
	for (int i = 0; i < N; i++)
	{
		cols[i] = other.cols[i];
		rows[i] = other.rows[i];
	}
	return *this;
}
void Matrix::Solve(double* B, double* x)
{
	int m = N;
	double* f = new double[m];
	double* b = new double[m];

	int i, j, n;
	double rk, sk, tk, Gk, Fk, r__;
	f[0] = 1 / rows[0];
	b[0] = f[0];
	x[0] = f[0] * B[0];
	for (n = 1; n < m; ++n)
	{
		Fk = 0;
		Gk = 0;
		for (j = 0; j < n; ++j)
		{
			Fk += (double)rows[n - j] * f[j];
			Gk += (double)cols[j + 1] * b[j];
		}
		rk = (double)1 / (1 - Gk * Fk);
		sk = (double)(-Fk * rk);
		tk = (double)(-Gk * rk);
		Fk = f[0];
		Gk = b[0];
		f[0] = Fk * rk;
		b[0] = Fk * tk;
		for (j = 1; j <= n - 1; ++j)
		{
			Fk = f[j];
			f[j] = (double)(Fk * rk + Gk * sk);
			r__ = b[j];
			b[j] = (double)(Fk * tk + Gk * rk);
			Gk = r__;
		}
		f[n] = (double)Gk * sk;
		b[n] = (double)Gk * rk;
		r__ = B[n];
		for (i = 0; i <= n - 1; i++)
		{
			r__ -= (double)rows[n - i] * x[i];
		}
		for (i = 0; i <= n - 1; i++)
		{
			x[i] += (double)b[i] * r__;
		}
		x[n] = b[n] * r__;
	}
	delete[] f;
	delete[] b;
}

void solutionAuto(int n, int repeat, double* B, double* x, double* duration)
{
	Matrix matrix(n);
	int start = clock();
	for (int i = 0; i < repeat; i++)
	{
		matrix.Solve(B, x);
	}
	int stop = clock();
	*duration = (double)(stop - start) / CLOCKS_PER_SEC * 1000;
}
void solution(int n, int repeat, double* B, double* x, double* duration, double* cols, double* rows)
{
	Matrix matrix(n, cols, rows);
	int start = clock();
	for (int i = 0; i < repeat; i++)
	{
		matrix.Solve(B, x);
	}
	int stop = clock();
	*duration = (double)(stop - start) / CLOCKS_PER_SEC * 1000;
}