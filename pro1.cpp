#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

void initTrueMatrix(double** &MatrixTrue);
void initFilterMatrix(double** &MatrixFilter);
double f(double x, double y);
void printMatrix(double** matrix);

int M = 10, N = 10;
double AX = 1.0, BX = 4.0, AY = 0.5, BY = 2.0, hX = (BX - AX) / double(N), hY = (BY - AY) / double(M);  

int main(){
	srand(time(NULL));
	double** MatrixTrue = new double*[M + 1];
	double** MatrixFilter = new double*[M + 1];
	for (int i = 0; i <= M; i++){
		MatrixTrue[i] = new double[N + 1];
		MatrixFilter[i] = new double[N + 1];
	}
	initTrueMatrix(MatrixTrue);
	initFilterMatrix(MatrixFilter);
	printMatrix(MatrixTrue);
	printMatrix(MatrixFilter);
	return 0;
}

void initTrueMatrix(double** &MatrixTrue){
	for (int i = 0; i <= M; i++){
		for (int j = 0; j <= N; j++){
			MatrixTrue[i][j] = f(AX + j * hX, AY + i * hY);
		}
	}
}

void initFilterMatrix(double** &MatrixFilter){
	for (int i = 0; i <= M; i++){
		for (int j = 0; j <= N; j++){
			MatrixFilter[i][j] = rand() % 2;
		}
	}
}

double f(double x, double y){
	return (x - AX) * (x - BX) * (log(cos(y) + 2.0) + tan(log(y)));
}

void printMatrix(double** matrix){
	for (int i = 0; i <= M; i++){
		printf("| %f ", matrix[i][0]);
		for (int j = 1; j < N; j++){
			printf("%f ", matrix[i][j]);
		}
		printf("%f |\n", matrix[i][N]);
	}
}