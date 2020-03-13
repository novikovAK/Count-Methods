//plot "out.txt" using 1:2 with lines title "f(x)", "out.txt" using 1:3 with lines title "L(x)", "out.txt" using 1:4 with lines title "G(x)"
//plot "out.txt" using 1:5 with lines title "f(x) - L(x)", "out.txt" using 1:6 with lines title "f(x) - G(x)", "out.txt" using 1:7 with lines title "L(x) - G(x)"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "gauss.cpp"

double f(double x);
void create_uniform_set(int N, double a, double h, double* &X, double* &FX);
void create_cheb_set(int N, double a, double b, double h, double* &X, double* &FX);
void create_lagrange_polinomians(FILE* out, int N, double h, double* &X, double* &FX, double* c);
void create_matrix_of_coeff(double** &MX, double* &B, int N, double* X, double* FX);
double f_LES(double* c, int N, double x);

int main(){
	int N = 20, a = -1.0, b = 1.0;
	double h = (b - a) / double(N - 1);
	double* X = new double[N];
	double* FX = new double[N];
	double** MX = new double*[N];
	double* B = new double[N];
	for (int i = 0; i < N; ++i){
		MX[i] = new double[N];
	}
	FILE *out = fopen("out.txt", "w");
	if (out == NULL){
		printf("File opening error: out\n");
		exit (-1);
	}
	//create_uniform_set(N, a, h, X, FX);
	create_cheb_set(N, a, b, h, X, FX);
	create_matrix_of_coeff(MX, B, N, X, FX);
	double* c = NULL;
	c = solve_LES_by_Gauss_with_squared_matrix(N, MX, B);
	printf("LES Solution:\n");
	for (int l = 0; l < N; ++l){
		if (l < N - 1){
			printf("%.2e*x^%d +\n", c[l], l);
		} else {
			printf("%.2e*x^%d\n", c[l], l);
		}
	}
	create_lagrange_polinomians(out, N, h, X, FX, c);
	fclose(out);
	return 0;
}

void create_lagrange_polinomians(FILE* out, int N, double h, double* &X, double* &FX, double* c){
	double L = 0.0;
	double* PHI = new double[N];
	double* Y = new double[3];
	double f_LES_v = 0.0;
	Y[0] = 0.0;
	Y[1] = 1.0 / 3.0;
	Y[2] = 2.0 / 3.0;
	for (int l = 0; l < N; ++l){
		for (int k = 0; k < 3; ++k){
			L = 0.0;
			for (int i = 0; i < N; ++i){
				PHI[i] = 1.0;
				for (int j = 0; j < N; ++j){
					if (j != i){
						if (l < N - 1){
							PHI[i] *= (X[l] + Y[k] * (X[l + 1] - X[l]) - X[j]);
						}
						if (fabs((X[i] - X[j])) < pow(10, -18)){
							printf("ALARM_1\n");
						}
						PHI[i] /= (X[i] - X[j]);
					}
				}
				L += FX[i] * PHI[i];
			}
			if (l < N - 1){
				f_LES_v = f_LES(c, N, X[l] + Y[k] * (X[l+1] - X[l]));
				fprintf(out, "%.5f %.5f %.5f %.5f %.10f %.10f %.10f\n", X[l] + Y[k] * (X[l+1] - X[l]), f(X[l] + Y[k] * (X[l+1] - X[l])), L, f_LES_v, f(X[l] + Y[k] * (X[l+1] - X[l])) - L, f(X[l] + Y[k] * (X[l+1] - X[l])) - f_LES_v, L - f_LES_v);
			}
			/*
			if (X[l] + Y[k] * h < 0){
				fprintf(out, "x = %.5f\t\tf-L = %.2e\t\tf-f_LES = %.2e\n", X[l] + Y[k] * h, f(X[l] + Y[k] * h) - L, f(X[l] + Y[k] * h) - f_LES_v);
			} else {
				fprintf(out, "x = +%.5f\t\tf-L = %.2e\t\tf-f_LES = %.2e\n", X[l] + Y[k] * h, f(X[l] + Y[k] * h) - L, f(X[l] + Y[k] * h) - f_LES_v);
			}
			*/
		}
	}
}

void create_cheb_set(int N, double a, double b, double h, double* &X, double* &FX){
	FILE* in = fopen("in.txt", "w");
	if (in == NULL){
		printf("File opening error: in\n");
		exit (-1);
	}
	fprintf(in, "%d\n", N);
	for (int i = 1; i <= N; ++i){
		X[i - 1] = (a + b) / 2.0 + ((b - a) / 2.0) * cos(((2 * i - 1) * M_PI) / (2 * N));
		FX[i - 1] = f(X[i - 1]);
		fprintf(in, "%f %f\n", X[i - 1], FX[i - 1]);
	}
	fclose(in);	
}

void create_uniform_set(int N, double a, double h, double* &X, double* &FX){
	FILE* in = fopen("in.txt", "w");
	if (in == NULL){
		printf("File opening error: in\n");
		exit (-1);
	}
	fprintf(in, "%d\n", N);
	for (int i = 0; i < N; ++i){
		X[i] = a + i*h;
		FX[i] = f(X[i]);
		fprintf(in, "%f %f\n", X[i], FX[i]);
	}
	fclose(in);
}

void create_matrix_of_coeff(double** &MX, double* &B, int N, double* X, double* FX){
	for (int i = 0; i < N; ++i){
		//printf("X[%d] = %f\n", i, X[i]);
		for (int j = 0; j < N; ++j){
			MX[i][j] = pow(X[i], j);
			//printf("MX[%d][%d] = %f\n", i, j, MX[i][j]);
		}
		B[i] = FX[i];
		//printf("B[%d] = %f\n", i, B[i]);
	}
}

double f(double x){
	return 1.0 / (25.0 * x * x + 1.0);
	//return pow(x, 49);
}

double f_LES(double* c, int N, double x){
	double temp = 0.0;
	for (int i = 0; i < N; ++i){
		temp += c[i] * pow(x, i);
	}
	return temp;
}
