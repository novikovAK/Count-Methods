#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double f(double x, double y);
double integrate(FILE* out, int n, double a, double b, double true_int, double h, double p);

int main(){
	FILE* out = fopen("out.txt", "w");
	int n = 0, max_n = 0;
	double a = 0.0, b = 1.0, true_int = 0.0, h = 0.0, p = 0.0, max_p = 0.0;
	true_int = 0.5 * (exp(1.0) - 1.0);
	for (n = 100; n <= 1000; n += 100){
		h = (b - a) / double(n);
		p = integrate(out, n, a, b, true_int, h, p);
		if (fabs(p) > max_p){
			max_p = fabs(p);
			max_n = n;
		}
	}
	fprintf(out, "\nmax p = %f for n = %d\n", max_p, max_n);
	fclose(out);
	return 0;
}

double integrate(FILE* out, int n, double a, double b, double true_int, double h, double p){
	double solve_int = 0.0;
	for (int i = 0; i < n; i++){
		for (int j = 0; j < n; j++){
			solve_int += 1.0 / 6.0 * (f((double(j) + 0.5) * h, double(i) * h) + f((double(j) + 1.0) * h, (double(i) + 0.5) * h) + f((double(j) + 0.5) * h, (double(i) + 0.5) * h)) * pow(h, 2);
			solve_int += 1.0 / 6.0 * ( f(double(j) * h, (double(i) + 0.5) * h) + f((double(j) + 0.5) * h, (double(i) + 1.0) * h) + f((double(j) + 0.5) * h, (double(i) + 0.5) * h)) * pow(h, 2);
		}
	}
	p = log(fabs(true_int - solve_int)) / log(double(h));
	fprintf(out, "n = %d, delta = %1.0e = (h = %f) ^ %f\n", n, fabs(true_int - solve_int), h, p);
	return p;
}

double f(double x, double y){
	return y * exp(x);
}