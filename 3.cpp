#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>

void Integral(double a, double b, int N);
void print(FILE* out, double delta, double estimate);
double max_dkf(double a, double b, int index, int k);
double dkf(double x, int index, int k);

double M = 1.0;
int functionsCount = 4;

int main(){
	double a = 1.0, b = 1.2;
	int N = 11;
	Integral(a, b, N);
	return 0;
}

void Integral(double a, double b, int N){
	FILE* out = fopen("out.txt", "w");
	double squadro = 0.0, true_int = 0.0, delta = 0.0, estimate = 0.0, x_0 = 0.0, x_m = 0.0, x_p = 0.0, ak = 0.0, bk = 0.0;
	std::vector <double> true_int_vector(functionsCount);
	std::vector <double> squadro_vector(functionsCount);
	for (int i = 0; i < functionsCount; ++i){
		true_int_vector.at(i) = 0.0;
		squadro_vector.at(i) = 0.0;
	}
	fprintf(out, "a = %f, b = %f, N = %d\n\n", a, b, N);
	for (int k = 0; k < N; ++k){
		ak = a + ((b - a) / double(N)) * double(k);
		bk = a + ((b - a) / double(N)) * double(k + 1);
		x_0 = 0.5 * (ak + bk), x_m = 0.5 * (ak + bk) - 0.5 * (bk - ak) * sqrt(3.0 / 5.0), x_p = 0.5 * (ak + bk) + 0.5 * (bk - ak) * sqrt(3.0 / 5.0);
		fprintf(out, "ak = %f, bk = %f\n\n", ak, bk);
		for (int index = 1; index <= functionsCount; ++index){
			switch (index){
				case 1: {
					fprintf(out, "f(x) = exp(x)\n"); 
					break;
				}
				case 2: {
					fprintf(out, "f(x) = sin(x)\n");
					break;
				}
				case 3: {
					fprintf(out, "f(x) = x^10\n");
					break;
				}
				case 4: {
					fprintf(out, "f(x) = sin^2(%.0f * x)\n", M);
					break;
				}
				default: exit (-1);
			}
			true_int = dkf(bk, index, -1) - dkf(ak, index, -1);
			true_int_vector.at(index - 1) += true_int;
			fprintf(out, "true integral = %.20f\n\n", true_int);
			// Poligons for the left point formula
			squadro = (bk - ak) * dkf(ak, index, 0);
			fprintf(out, "Poligons for the left point formula: %.20f\n", squadro);
			delta = fabs(squadro - true_int);
			estimate = 0.5 * pow(bk - ak, 2) * max_dkf(ak, bk, index, 1);
			print(out, delta, estimate);
			// Poligons for the central point formula
			squadro = (bk - ak) * dkf(x_0, index, 0);
			fprintf(out, "Poligons for the central point formula: %.20f\n", squadro);
			delta = fabs(squadro - true_int);
			estimate = (1.0 / 24.0) * pow(bk - ak, 3) * max_dkf(ak, bk, index, 2);
			print(out, delta, estimate);
			// Trapezes formula
			squadro = 0.5 * (bk - ak) * (dkf(ak, index, 0) + dkf(bk, index, 0));
			fprintf(out, "Trapezes formula: %.20f\n", squadro);
			delta = fabs(squadro - true_int);
			estimate = (1.0 / 12.0) * pow(bk - ak, 3) * max_dkf(ak, bk, index, 2);
			print(out, delta, estimate);
			// Simpson's formula
			squadro = (1.0 / 6.0) * (bk - ak) * (dkf(ak, index, 0) + 4.0 * dkf(x_0, index, 0) + dkf(bk, index, 0));
			fprintf(out, "Simpson's formula: %.20f\n", squadro);
			delta = fabs(squadro - true_int);
			estimate = (1.0 / 2880.0) * pow(bk - ak, 5) * max_dkf(ak, bk, index, 4);
			print(out, delta, estimate);
			// Gauss's three nodes formula
			squadro = (1.0 / 18.0) * (bk - ak) * (5.0 * dkf(x_m, index, 0) + 8.0 * dkf(x_0, index, 0) + 5.0 * dkf(x_p, index, 0));
			squadro_vector.at(index - 1) += squadro;
			fprintf(out, "Gauss's three nodes formula: %.20f\n", squadro);
			delta = fabs(squadro - true_int);
			estimate = (1.0 / 737280.0) * pow(bk - ak, 7) * max_dkf(ak, bk, index, 6);
			print(out, delta, estimate);
		}
	}
	if (N > 1){
		fprintf(out, "\n");
		for (int index = 1; index <= functionsCount; ++index){
			switch (index){
				case 1: {
					fprintf(out, "f(x) = exp(x)\n"); 
					break;
				}
				case 2: {
					fprintf(out, "f(x) = sin(x)\n");
					break;
				}
				case 3: {
					fprintf(out, "f(x) = x^10\n");
					break;
				}
				case 4: {
					fprintf(out, "f(x) = sin^2(%.0f * x)\n", M);
					break;
				}
				default: exit (-1);
			}
			fprintf(out, "Integral for [%f, %f]:\nGauss's three nodes formula: %.20f\nTrue integral: %.20f\ndelta = %.20f\n\n", a, b, squadro_vector.at(index - 1), true_int_vector.at(index - 1), fabs(squadro_vector.at(index - 1) - true_int_vector.at(index - 1)));
		}
		if (N > 10){
			double del = fabs(squadro_vector.at(3) - true_int_vector.at(3));
			printf("del = %f\n", del);
			double R = log(1.0 / del);
			printf("R = %f\n", R);
			printf("int(R) = %d int(log(N)) = %d\n", int(R), int(log(N)));
			double lnC = double(int(R) % int(log(N)));
			printf("lnC = %f\n", lnC);
			double C = 1.0 / exp((lnC));
			printf("C = %f\n", C);
			double p = (R - lnC) / log(N);
			printf("p = %f\n", p);
			double acc = C / exp(p * log(N));
			printf("acc = %f\n", acc);
			fprintf(out, "%.30f (delta) ~ %.10f (C) / %d (N) ^ %.10f (p) = %.30f\naccuracy = %.30f\n", del, C, N, p, acc, fabs(del - acc));
		}
	}
	fclose(out);
}

void print(FILE* out, double delta, double estimate){
	if (delta < estimate){
		fprintf(out, "%.20f (delta) <  %.20f (estimate):\t\tOK\n\n", delta, estimate);
	} else {
		fprintf(out, "%.20f (delta) >  %.20f (estimate):\t\tERROR\n\n", delta, estimate);
	}
}

double max_dkf(double a, double b, int index, int k){
	int n = 100;
	double maximum = 0.0;
	for (int i = 0; i < n; i++){
		if (fabs(dkf(a + (double(i) / double(n)) * (b - a), index, k)) > maximum){
			maximum = fabs(dkf(a + (double(i) / double(n)) * (b - a), index, k));
		}
	}
	return maximum;
}

double dkf(double x, int index, int k){
	switch (index){
		case 1: {
			switch (k){
				case -1: return exp(x);
				case 0: return exp(x);
				case 1: return exp(x);
				case 2: return exp(x);
				case 4: return exp(x);
				case 6: return exp(x);
				default: exit (-1);
			}
		}
		case 2: {
			switch (k){
				case -1: return -1.0 * cos(x);
				case 0: return sin(x);
				case 1: return cos(x);
				case 2: return -1.0 * sin(x);
				case 4: return sin(x);
				case 6: return -1.0 * sin(x);
				default: exit (-1);
			}
		}
		case 3: {
			switch (k){
				case -1: return pow(x, 11) / 11.0;
				case 0: return pow(x, 10);
				case 1: return 10.0 * pow(x, 9);
				case 2: return 90.0 * pow(x, 8);
				case 4: return 5040.0 * pow(x, 6);
				case 6: return 151200.0 * pow(x, 4);
				default: exit (-1);
			}
		}
		case 4: {
			switch (k){
				case -1: return 0.5 * x - (sin(2.0 * M * x)) / (4.0 * M);
				case 0: return pow(sin(M * x), 2);
				case 1: return 2.0 * M * sin(M * x) * cos(M * x);
				case 2: return 2.0 * pow(M, 2) * cos(2.0 * M * x);
				case 4: return -8.0 * pow(M, 4) * cos(2.0 * M * x);
				case 6: return 32.0 * pow(M, 6) * cos(2.0 * M * x);
				default: exit (-1);
			}
		}
		default: exit (-1);
	}
}
