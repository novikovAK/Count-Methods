//plot "out.txt" using 1:2 with lines title "y(x)", "out.txt" using 1:3 with lines title "F(x)"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "gauss.cpp"
#include "qr.cpp"

void create_uniform_set(int M, double a, double h, double* &X, double* &Y);
void create_matrix_of_coeff(double** &A, int M, int N, double* X);
double f(double x);
void check(FILE* out, int M, int N, double* X, double* c);

int main(){
	int M = 10, N = 10, a = -1.0, b = 1.0;
	double h = (b - a) / double(M), temp = 0.0;
	double* X = new double[M + 1];
	double* Y = new double[M + 1];
	double** A = new double*[M + 1];
	double* B = new double[M + 1];
	for (int i = 0; i <= M; ++i){
		A[i] = new double[N + 1];
	}
	FILE *out = fopen("out.txt", "w");
	if (out == NULL){
		printf("File opening error: out\n");
		exit (-1);
	}
	create_uniform_set(M, a, h, X, Y);
	create_matrix_of_coeff(A, M, N, X);
	double** Q = new double*[M + 1];
	double** R = new double*[M + 1];
	for (int i = 0; i <= M; ++i){
		Q[i] = new double[M + 1];
		R[i] = new double[N + 1];
	}
	getQRdecomposition(M + 1, N + 1, Q, R, A);
	for (int i = 0; i <= M; ++i){
		temp = 0.0;
		for (int j = 0; j <= M; ++j){
			temp += Q[i][j] * Y[j];
		}
		B[i] = temp;
	}
	double** G = new double*[N + 1];
	double* BG = new double[N + 1];
	for (int i = 0; i <= N; ++i){
		G[i] = new double[N + 1];
	}
	for (int i = 0; i <= N; ++i){
		for (int j = 0; j <= N; ++j){
			G[i][j] = R[i][j];
		}
		BG[i] = B[i];
	}
	double* c = solve_LES_by_Gauss_with_squared_matrix(N + 1, G, BG);
	printf("LES Solution:\n");
	for (int l = 0; l <= N; ++l){
		if (l == 0){
			//printf("%.2e +\n", c[l]);
			printf("%.16f +\n", c[l]);
		} else if (l < N){
			//printf("%.2e*x^%d +\n", c[l], l);
			printf("%.16f * x^%d +\n", c[l], l);
		} else {
			//printf("%.2e*x^%d\n", c[l], l);
			printf("%.16f * x^%d\n", c[l], l);
		}
	}
	check(out, M, N, X, c);
	fclose(out);
	return 0;
}

void create_uniform_set(int M, double a, double h, double* &X, double* &Y){
	FILE* in = fopen("in.txt", "w");
	if (in == NULL){
		printf("File opening error: in\n");
		exit (-1);
	}
	fprintf(in, "%d\n", M);
	for (int i = 0; i <= M; ++i){
		X[i] = a + i*h;
		Y[i] = f(X[i]);
		fprintf(in, "%f %f\n", X[i], Y[i]);
	}
	fclose(in);
}

void create_matrix_of_coeff(double** &A, int M, int N, double* X){
	for (int i = 0; i <= M; ++i){
		//printf("X[%d] = %f\n", i, X[i]);
		for (int j = 0; j <= N; ++j){
			A[i][j] = pow(X[i], j);
			//printf("A[%d][%d] = %f\n", i, j, A[i][j]);
		}
	}
}

double f(double x){
	return 1.0 / (25.0 * x * x + 1.0);
	//return pow(x, 10) + 1.0;
}

void check(FILE* out, int M, int N, double* X, double* c){
	double* H = new double[3];
	H[0] = 0.0;
	H[1] = 1.0 / 3.0;
	H[2] = 2.0 / 3.0;
	double temp = 0.0;
	for (int l = 0; l < M; ++l){
		for (int k = 0; k < 3; ++k){
			temp = 0.0;
			for (int i = 0; i <= N; ++i){
				temp += c[i] * pow(X[l] + H[k] * (X[l + 1] - X[l]), i);
			}
			fprintf(out, "%.5f %.5f %.5f %.20f\n", X[l] + H[k] * (X[l + 1] - X[l]), f(X[l] + H[k] * (X[l + 1] - X[l])), temp, f(X[l] + H[k] * (X[l + 1] - X[l])) - temp);
		}
	}
}
