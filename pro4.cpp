#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

void calculateSpeedOfConvergence(FILE* out, int index, int n, int m, double h_x, double h_y, double** MatrixTrue, double** MatrixRecovery);
double** splainCubedRecover(int n, int m, double a_y, double h_y, double** MatrixTrue);
void initTrueMatrix(int n, int m, double a_x, double b_x, double h_x, double a_y, double b_y, double h_y, double** &MatrixTrue);
double f(double a_x, double b_x, double x, double y);
double* solve_LES_by_Gauss_with_squared_matrix(int N, double** A, double* b);
void printMatrixes(FILE* out, int n, int m, double a_x, double h_x, double a_y, double h_y, double** MatrixTrue, double** MatrixRecovery);

int main(){
	srand(time(NULL));
	FILE* out = fopen("out.txt", "w");
	if (out == NULL){
		printf("out.txt file opening error\n");
		return -1;
	}
	int n = 0, m = 0;
	double a_x = 2.0, b_x = 6.0, a_y = 1.0, b_y = 4.0, h_x = 0.0, h_y = 0.0; 
	double** MatrixTrue = NULL;
	double** MatrixRecovery = NULL;
	int** MatrixFilter = NULL;
	for (n = 10; n <= 10; n += 10){
		for (m = 10; m <= 10; m += 10){
			h_x = (b_x - a_x) / double(n);
			h_y = (b_y - a_y) / double(m);
			MatrixTrue = new double*[m + 1];
			MatrixFilter = new int*[m + 1];
			for (int i = 0; i <= m; i++){
				MatrixTrue[i] = new double[n + 1];
				MatrixFilter[i] = new int[n + 1];
			}
			initTrueMatrix(n, m, a_x, b_x, h_x, a_y, b_y, h_y, MatrixTrue);
			initFilterMatrix(n, m, MatrixFilter);
			MatrixRecovery = splainCubedRecover(n, m, a_y, h_y, MatrixTrue);
			//MatrixRecovery = furieRecover(n, m, a_x, b_x, h_x, MatrixTrue, MatrixFilter);
			printMatrixes(out, n, m, a_x, h_x, a_y, h_y, MatrixTrue, MatrixRecovery);
			calculateSpeedOfConvergence(out, 0, n, m, h_x, h_y, MatrixTrue, MatrixRecovery);
		}
	}
	fclose(out);
	return 0;
}

void calculateSpeedOfConvergence(FILE* out, int index, int n, int m, double h_x, double h_y, double** MatrixTrue, double** MatrixRecovery){
	double sum = 0.0, p = 0.0;
	int count = 0;
	for (int i = 0; i <= m; i++){
		for (int j = 0; j <= n; j++){
			count++;
			sum += fabs(MatrixTrue[i][j] - MatrixRecovery[i][j]);		
		}
	}
	sum /= double(count);
	fprintf(out, "n = %d\tm = %d\t", n, m);
	if (sum > pow(10, -15)){
		if (index == 0){
			p = log(sum) / log (h_y);
		} else {
			p = log(sum) / log (h_x);
		}
		fprintf(out, "p = %.10f\n", p);
	} else {
		fprintf(out, "The sum of deltas is less than 10 ^ -15\n");
	}
}

double** splainCubedRecover(int n, int m, double a_y, double h_y, double** MatrixTrue){
	int N = 0, tempInt = 0;
	double *X = NULL, *Y = NULL, *M = NULL, *P = NULL;
	double **C = NULL, **MatrixRecovery;
	MatrixRecovery = new double*[m + 1];
	for (int i = 0; i <= m; ++i){
		MatrixRecovery[i] = new double[n + 1];
	}
	for (int k = 0; k <= n; k++){
		// Make the available points count zero
		N = 0;
		// Figure out how many points are available
		for (int s = 0; s <= m; s++){
			if (MatrixFilter[s][k] == 1){
				MatrixRecovery[s][k] = MatrixTrue[s][k];
				N++;
			}
		}
		// Create X-vector contained available points and Y-vector contained the function value in these available points
		X = new double[N];
		Y = new double[N];
		tempInt = 0;
		for (int s = 0; s <= m; s++){
			if (MatrixFilter[s][k] == 1){
				Y[tempInt] = MatrixTrue[s][k];
				X[tempInt] = a_y + s * h_y;
				tempInt++;
			}
		}
		// Create C-matrix for cubed splaines
		C = new double*[N - 2];
		for (int i = 0; i < N -2; i++){
			C[i] = new double[N - 2];
		}
		// Fll C-matrix for cubed splaines like in the report
		for (int i = 0; i < N -2; i++){
			for (int j = 0; j < N - 2; j++){
				if (i == j){
					C[i][j] = (X[i + 2] - X[i]) / 3.0;
				} else if (j == i + 1){
					C[i][j] = (X[j + 1] - X[j]) / 6.0;
				} else if (i == j + 1){
					C[i][j] = (X[i + 1] - X[i]) / 6.0;
				} else {
					C[i][j] = 0.0;
				}
			}
		}
		// Create P-vector for the right part of Linear Equations System CM = P
		P = new double[N - 2];
		for (int i = 0; i < N -2; i++){
			P[i] = (Y[i + 2] - Y[i + 1]) / (X[i + 2] - X[i + 1]) - (Y[i + 1] - Y[i]) / (X[i + 1] - X[i]);
		}
		// M - is the vector of the second derivatives for the solution of Linear Equations System CM = P, but M is Nx1 vector becouse M[0] = M[N - 1] = 0 but P and the solution of Linear Equations System CM = P are (N - 2)x1 vector. So we put the solution of Linear Equations System CM = P into P 
		P = solve_LES_by_Gauss_with_squared_matrix(N - 2, C, P);
		// Initialise M like Nx1 vector
		M = new double[N];
		// Add the extreme values
		M[0] = 0.0;
		M[N - 1] = 0.0;
		// Copy P into M
		for (int i = 1; i < N - 1; ++i){
			M[i] = P[i - 1];
		}
		// Recover uknown values using splaines
		tempInt = 0;
		for (int s = 0; s <= m; s++){
			if (MatrixFilter[s][k] == 1){
				++tempInt;
			} else {
				MatrixRecovery[s][k] = M[tempInt - 1] * pow(X[tempInt] - (a_y + s * h_y), 3) / (6.0 * (X[tempInt] - X[tempInt - 1])) + M[tempInt] * pow((a_y + s * h_y) - X[tempInt - 1], 3) / (6.0 * (X[tempInt] - X[tempInt - 1])) + (Y[tempInt - 1] - (M[tempInt - 1] * pow(X[tempInt] - X[tempInt - 1], 2)) / 6.0) * (X[tempInt] - (a_y + s * h_y)) / (X[tempInt] - X[tempInt - 1]) + (Y[tempInt] - (M[tempInt] * pow(X[tempInt] - X[tempInt - 1], 2)) / 6.0) * (a_y + s * h_y - X[tempInt - 1]) / (X[tempInt] - X[tempInt - 1]);
			}
		}
	}
	return MatrixRecovery;
}

void initTrueMatrix(int n, int m, double a_x, double b_x, double h_x, double a_y, double b_y, double h_y, double** &MatrixTrue){
	for (int i = 0; i <= m; i++){
		for (int j = 0; j <= n; j++){
			MatrixTrue[i][j] = f(a_x, b_x, a_x + j * h_x, a_y + i * h_y);
		}
	}
}

double f(double a_x, double b_x, double x, double y){
	return (x - a_x) * (x - b_x) * (cos(x) - log(y * y));
}

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

void printMatrixes(FILE* out, int n, int m, double a_x, double h_x, double a_y, double h_y, double** MatrixTrue, int** MatrixFilter, double** MatrixRecovery){
	fprintf(out, "(X, Y): \n");
	for (int i = 0; i <= m; i++){
		for (int j = 0; j <= n; j++){
			fprintf(out, "(%.3f, %.3f)\t", a_x + j * h_x, a_y + i * h_y) ;
		}
		fprintf(out, "\n");
	}
	fprintf(out, "\n");
	fprintf(out, "MatrixTrue: \n");
	for (int i = 0; i <= m; i++){
		fprintf(out, "|\t%3.6f\t", MatrixTrue[i][0]);
		for (int j = 1; j < n; j++){
			fprintf(out, "%3.6f\t", MatrixTrue[i][j]);
		}
		fprintf(out, "%3.6f\t|\n", MatrixTrue[i][n]);
	}
	fprintf(out, "\n");
	fprintf(out, "MatrixFilter: \n");
	for (int i = 0; i <= m; i++){
		fprintf(out, "|\t%d\t", MatrixFilter[i][0]);
		for (int j = 1; j < n; j++){
			fprintf(out, "%d\t", MatrixFilter[i][j]);
		}
		fprintf(out, "%d\t|\n", MatrixFilter[i][n]);
	}
	fprintf(out, "\n");
	fprintf(out, "MatrixRecovery: \n");
	for (int i = 0; i <= m; i++){
		fprintf(out, "|\t%3.6f\t", MatrixRecovery[i][0]);
		for (int j = 1; j < n; j++){
			fprintf(out, "%3.6f\t", MatrixRecovery[i][j]);
		}
		fprintf(out, "%3.6f\t|\n", MatrixRecovery[i][n]);
	}
	fprintf(out, "\n");
	fprintf(out, "MatrixDelta: \n");
	for (int i = 0; i <= m; i++){
		fprintf(out, "|\t%3.8f\t", fabs(MatrixTrue[i][0] - MatrixRecovery[i][0]));
		for (int j = 1; j < n; j++){
			fprintf(out, "%3.8f\t", fabs(MatrixTrue[i][j] - MatrixRecovery[i][j]));
		}
		fprintf(out, "%3.8f\t|\n", fabs(MatrixTrue[i][n] - MatrixRecovery[i][n]));
	}
	fprintf(out, "\n");
}