#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double phiFunction(double a, double b, int m, double x);
double f(double a, double b, double x);

int main(){
	int n = 100;
	double a = 2.0, b = 6.0;
	double h = (b - a) / double(n);
	double *c = NULL;
	c = new double[n];
	for (int m = 1; m < n; ++m){
		c[m] = 0.0;
		for (int j = 1; j < n; j++){
			c[m] += (1.0 / double(n)) * f(a, b, a + double(j) * h) * phiFunction(a, b, m, a + double(j) * h);
		}
	}
	double temp = 0.0;
	for (int j = 0; j < n; j++){
		temp = 0.0;
		for (int m = 1; m < n; ++m){
			temp += c[m] * phiFunction(a, b, m, a + (double(j) + 0.5) * h);
		}
		printf("%f (True) - %f (Solved) = %f\n", f(a, b, a + (double(j) + 0.5) * h), temp, f(a, b, a + (double(j) + 0.5) * h) - temp);
	}
	return 0;
}

double phiFunction(double a, double b, int m, double x){
	return sqrt(2.0) * sin(M_PI * double(m) * (x - a) / (b - a)); 
}

double f(double a, double b, double x){
	return (x - a) * (x - b) * cos(x);
}