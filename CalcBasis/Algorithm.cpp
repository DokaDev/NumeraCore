#include "pch.h"
#include <vector>
using namespace std;

extern "C" __declspec(dllexport) double Derivative(double (*f)(double), double x, double h) {
	return (f(x + h) - f(x - h)) / (2 * h);
}

extern "C" __declspec(dllexport) double TrapezoidalRule(double (*f)(double), double a, double b, int n) {
	double h = (b - a) / n;
	double sum = 0.5 * (f(a) + f(b));
	for(int i = 1; i < n; i++) {
		sum += f(a + i * h);
	}
	return h * sum;
}

extern "C" __declspec(dllexport) void GaussianElimination(vector<vector<double>>& mat) {
	int n = mat.size();
	for(int i = 0; i < n; i++) {
		// Pivot
		for(int j = i + 1; j < n; j++) {
			double factor = mat[j][i] / mat[i][i];
			for(int k = i; k <= n; k++) {
				mat[j][k] -= factor * mat[i][k];
			}
		}
	}

	// Back substitution
	vector<double> x(n);
	for(int i = n - 1; i >= 0; i--) {
		x[i] = mat[i][n];
		for(int j = i + 1; j < n; j++) {
			x[i] -= mat[i][j] * x[j];
		}
		x[i] = x[i] / mat[i][i];
	}
}

extern "C" __declspec(dllexport) double GradientDescent(double (*f)(double), double (*df)(double), double initial, double learningRate, int iterations) {
	double x = initial;
	for(int i = 0; i < iterations; i++) {
		x -= learningRate * df(x);
	}
	return x;
}