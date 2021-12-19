#pragma once

#ifdef MATRIXSOLVER_EXPORTS
#define MATRIXSOLVER_API __declspec(dllexport)
#else
#define MATRIXSOLVER_API __declspec(dllimport)
#endif

class Matrix
{
private:
	double* cols;
	double* rows;
	int N;

public:
	Matrix(int);
	Matrix(int, double*, double*);
	Matrix(const Matrix&);
	~Matrix();
	Matrix& operator=(const Matrix&);
	void Solve(double*, double*);
};


extern "C" MATRIXSOLVER_API void solutionAuto(int n, int repeat, double* B, double* x, double* duration);
extern "C" MATRIXSOLVER_API void solution(int n, int repeat, double* B, double* x, double* duration, double* cols, double* rows);