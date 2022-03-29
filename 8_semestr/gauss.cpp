#include <stdio.h>
#include <stdlib.h>
#include <math.h>


double* solve_LES_by_Gauss_with_squared_matrix(int N, double** A, double* b);

/*
int main(){
	FILE* in = fopen("in.txt", "r");
	int N = 0;
	fscanf(in, "%d", &N);
	double** A = new double*[N];
	double* b = new double[N];
	double* x = new double[N];
	for (int i = 0; i < N; ++i){
		A[i] = new double[N];
	}
	for (int i = 0; i < N; ++i){
		for (int j = 0; j < N; ++j){
			fscanf(in, "%lf", &A[i][j]);
		}
		fscanf(in, "%lf", &b[i]);
	}
	x = solve_LES_by_Gauss_with_squared_matrix(N, A, b);
	for (int i = 0; i < N; ++i){
		printf("%f\n", x[i]);
	}
	return 0;
}
*/ 

double* solve_LES_by_Gauss_with_squared_matrix(int N, double** A, double* b){
	int i = 0, j = 0, r = 0;
	double temp = 0.0;
	while (i < N && j < N){
		temp = A[i][j];
		for (int k = 0; k < N; ++k){
			A[i][k] /= temp;
		}
		b[i] /= temp;
		for (int l = 0; l < N; ++l){
			if (l != i){
				temp = A[l][j];
				for (int k = 0; k < N; ++k){
					A[l][k] -= temp * A[i][k];
				}
				b[l] -= temp * b[i];
			}
		}
		i++;
		j++;
	}
	return b;
}
