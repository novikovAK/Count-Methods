#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void getQRdecomposition(int M, int N, double** &Q, double** &R, double** A);
double sign(double x);
double scalar_multiplication(double* v, int N, double* u, int M);
double* normalize(int N, double* v);
double** matrixes_mulriplication(double** A, int M, int N, double** B, int P, int Q);
double** vector_times_vector_T(double* v, int N, double* u, int M);
double deltaIJint(int x, int y);
double deltaIJdouble(double x, double y);
/*
int main(){
	FILE* in = fopen("in.txt", "r");
	int N = 0, M = 0;
	if (fscanf(in, "%d%d", &M, &N) != 2){
		printf("Error (main): cannot read N or M\n");
		return 0;
	}
	double** A = new double*[M]; // M x N - matrix
	double** A_initial = new double*[M]; // M x N - matrix
	double** Q = new double*[M]; // M x M - matrix
	double** R = new double*[M]; // M x N - matrix
	for (int i = 0; i < M; ++i){
		A_initial[i] = new double[N];
		A[i] = new double[N];
		R[i] = new double[N];
		Q[i] = new double[M];
	}
	for (int i = 0; i < M; ++i){
		for (int j = 0; j < N; ++j){
			if (fscanf(in, "%lf", &A[i][j]) != 1){
				printf("Error (main): cannot read A[%d][%d]\n", i, j);
				return 0;
			}
			A_initial[i][j] = A[i][j];
		}
	}
	printf("A:\n");
	for (int i = 0; i < M; ++i){
		for (int j = 0; j < N; ++j){
			printf("%f\t", A[i][j]);
		}
		printf("\n");
	}
	getQRdecomposition(M, N, Q, R, A);
	printf("Q:\n");
	for (int i = 0; i < M; ++i){
		for (int j = 0; j < M; ++j){
			printf("%f\t", Q[i][j]);
		}
		printf("\n");
	}
	printf("R:\n");
	for (int i = 0; i < M; ++i){
		for (int j = 0; j < N; ++j){
			printf("%f\t", R[i][j]);
		}
		printf("\n");
	}
	printf("Q^T * R:\n");
	double temp = 0.0;
	for (int i = 0; i < M; ++i){
		for (int j = 0; j < N; ++j){
			temp = 0.0;
			for (int k = 0; k < M; ++k){
				temp += Q[k][i] * R[k][j];
			}
			printf("%f\t", temp);
		}
		printf("\n");
	}
	printf("Q * A:\n");
	for (int i = 0; i < M; ++i){
		for (int j = 0; j < N; ++j){
			temp = 0.0;
			for (int k = 0; k < M; ++k){
				temp += Q[i][k] * A_initial[k][j];
			}
			printf("%f\t", temp);
		}
		printf("\n");
	}
	return 0;
}
*/
void getQRdecomposition(int M, int N, double** &Q, double** &R, double** A){
	if (N > M){
		printf("Error(getQRdecomposition): This function is only for M x N matrixes where M >= N. M = %d, N = %d\n", M, N);
		exit (-1);
	}
	double alpha = 0.0;
	double* w = new double[M];
	double** U = new double*[M];
	for (int i = 0; i < M; ++i){
		U[i] = new double[M];
	}
	for (int i = 0; i < M; ++i){
		for (int j = 0; j < M; ++j){
			Q[i][j] = deltaIJint(i, j);
		}
	}
	for (int i = 0; i < N; ++i){
		alpha = 0.0;
		for (int k = 0; k < M; ++k){
			w[k] = 0.0;
		}
		for (int k = i; k < M; ++k){
			alpha += A[k][i] * A[k][i];
		}
		alpha = sqrt(alpha);
		w[i] = A[i][i] - alpha;
		for (int k = i + 1; k < M; ++k){
			w[k] = A[k][i];
		}
		normalize(M, w);
		U = vector_times_vector_T(w, M, w, M);
		for (int k = 0; k < M; ++k){
			for (int l = 0; l < M; ++l){
				U[k][l] = deltaIJint(k, l) - 2.0 * U[k][l];
			}
		}
		Q = matrixes_mulriplication(U, M, M, Q, M, M);
		A = matrixes_mulriplication(U, M, M, A, M, N);
	}
	for (int i = 0; i < M; ++i){
		for (int j = 0; j < N; ++j){
			R[i][j] = A[i][j];
		}
	}
}

double deltaIJint(int x, int y){
	if (x == y){
		return 1.0;
	} else {
		return 0.0;
	}
}

double deltaIJdouble(double x, double y){
	if (fabs(x-y) < pow(10, -20)){
		return 1.0;
	} else {
		return 0.0;
	}
}

double sign(double x){
	if (x < 0){
		return -1.0;
	} else if (x > 0){
		return 1.0;
	} else {
		return 0.0;
	}
}

double* normalize(int N, double* v){
	double length = sqrt(scalar_multiplication(v, N, v, N));
	if (length < pow(10, -20)){
		for (int i = 0; i < N; ++i){
			v[i] = 0.0;
		}
	} else {
		for (int i = 0; i < N; ++i){
			v[i] /= length;
		}
	}
	return v;
}

double scalar_multiplication(double* v, int N, double* u, int M){
	if (N != M){
		printf("Error (scalar_multiplication): Sizes of vectors aren't equal\n");
		exit (-1);
	}
	double temp = 0.0;
	for (int i = 0; i < N; ++i){
		temp += u[i] * v[i];
	}
	return temp;
}

double** matrixes_mulriplication(double** A, int M, int N, double** B, int P, int Q){
	if (N != P){
		printf("Error(matrixes_mulriplication): Sizes of matrixes aren't equal\n");
		exit (-1);
	}
	double temp = 0.0;
	double** C = new double*[M];
	for (int i = 0; i < M; ++i){
		C[i] = new double[P];
	}
	for (int i = 0; i < M; ++i){
		for (int j = 0; j < P; ++j){
			temp = 0.0;
			for (int s = 0; s < N; ++s){
				temp += A[i][s] * B[s][j];
			}
			C[i][j] = temp;
		}
	}
	return C;
}

double** vector_times_vector_T(double* v, int N, double* u, int M){
	if (N != M){
		printf("Error(vector_times_vector_T): Sizes of vectors aren't equal\n");
		exit (-1);
	}
	double** C = new double*[N];
	for (int i = 0; i < N; ++i){
		C[i] = new double[N];
	}
	for (int i = 0; i < N; ++i){
		for (int j = 0; j < N; ++j){
			C[i][j] = v[i] * u[j];
		}
	}
	return C;
}
