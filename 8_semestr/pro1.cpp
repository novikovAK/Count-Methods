#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

void calculateSpeedOfConvergence(FILE* out, int n, int m, double a_x, double b_x, double a_y, double b_y, double** MatrixTrue, int** MatrixFilter, double** MatrixRecovery);
double** furieRecover(int n, int m, double a_x, double b_x, double h_x, double** MatrixTrue, int** MatrixFilter);
double recoverValueScale(double below, double above, double t, double Y_below, double Y_above, double Y_q, double X_below, double X_above, double q);
double phiFunction(double a_x, double b_x, int p, double x);
double** splainHorCubedRecover(int n, int m, double a_x, double h_x, double** MatrixTrue, int** MatrixFilter);
double** splainVerCubedRecover(int n, int m, double a_y, double h_y, double** MatrixTrue, int** MatrixFilter);
void initTrueMatrix(int n, int m, double a_x, double b_x, double h_x, double a_y, double b_y, double h_y, double** &MatrixTrue);
void initFilterMatrix(int n, int m, int** &MatrixFilter);
double f(double a_x, double b_x, double x, double y);
double* solve_LES_by_Gauss_with_squared_matrix(int N, double** A, double* b);
void printMatrixes(FILE* out, int n, int m, double a_x, double h_x, double a_y, double h_y, double** MatrixTrue, int** MatrixFilter, double** MatrixRecovery);

int functionIndex = 1;

int main(){
	srand(time(NULL));
	FILE* out = fopen("out.txt", "w");
	if (out == NULL){
		printf("out.txt file opening error\n");
		return -1;
	}
	int n = 0, m = 0; // Количество разбиений по оси ОХ и ОУ соответственно
	double a_x = 0.0, b_x = 1.0, a_y = 0.0, b_y = 1.0, h_x = 0.0, h_y = 0.0; // Задаём область и шаг равномерного разбиения
	double** MatrixTrue = NULL; // Матрица истинных значений функции
	//double** MatrixRecoveryVer = NULL; // Матрица восстановленных значений функции
	double** MatrixRecoveryHor = NULL; // Матрица восстановленных значений функции
	int** MatrixFilter = NULL; // Матрица индикаторов
	for (functionIndex = 2; functionIndex <= 2; functionIndex++){
		//fprintf(out, "functionIndex: %d ", functionIndex);
		for (n = 200; n <= 500; n += 50){
			for (m = 3; m <= 3; m += 1){
				h_x = (b_x - a_x) / double(n); // Вычисляем шаг по оси ОХ
				h_y = (b_y - a_y) / double(m); // Вычисляем шаг по оси ОУ
				MatrixTrue = new double*[m + 1]; // Инициализируем матрицу истинных значений функции
				MatrixFilter = new int*[m + 1]; // Инициализируем матрицу индикаторов
				for (int i = 0; i <= m; i++){
					MatrixTrue[i] = new double[n + 1];
					MatrixFilter[i] = new int[n + 1];
				}
				initTrueMatrix(n, m, a_x, b_x, h_x, a_y, b_y, h_y, MatrixTrue); // Заполняем матрицу истинных значений функции
				initFilterMatrix(n, m, MatrixFilter); // Заполняем матрицу индикаторов
				MatrixRecoveryHor = splainHorCubedRecover(n, m, a_x, h_x, MatrixTrue, MatrixFilter); // Заполняем матрицу восстановленных значений функции методом сплайнов
				//MatrixRecoveryVer = splainVerCubedRecover(n, m, a_y, h_y, MatrixTrue, MatrixFilter); // Заполняем матрицу восстановленных значений функции методом сплайнов
				//for (int i = 0; i <= m; i++){
				//	for (int j = 0; j <= n; j++){
				//		MatrixRecoveryHor[i][j] = 0.5 * (MatrixRecoveryVer[i][j] + MatrixRecoveryHor[i][j]);
				//	}
				//}
				//MatrixRecovery = furieRecover(n, m, a_x, b_x, h_x, MatrixTrue, MatrixFilter); // Заполняем матрицу восстановленных значений функции методом Фурье
				//printMatrixes(out, n, m, a_x, h_x, a_y, h_y, MatrixTrue, MatrixFilter, MatrixRecoveryHor); // Печатем все матрицы
				//fprintf(out, "Hor: "); 
				//calculateSpeedOfConvergence(out, n, m, a_x, b_x, a_y, b_y, MatrixTrue, MatrixFilter, MatrixRecoveryHor);
				//printMatrixes(out, n, m, a_x, h_x, a_y, h_y, MatrixTrue, MatrixFilter, MatrixRecoveryVer); 
				//fprintf(out, "Ver: ");
				//calculateSpeedOfConvergence(out, n, m, a_x, b_x, a_y, b_y, MatrixTrue, MatrixFilter, MatrixRecoveryVer); // Вычисляем степень p для отчётности
				//for (int i = 0; i <= m; i++){
				//	for (int j = 0; j <= n; j++){
				//		MatrixRecoveryHor[i][j] = 0.5 * (MatrixRecoveryVer[i][j] + MatrixRecoveryHor[i][j]);
				//	}
				//}
				//fprintf(out, "Mid: ");
				calculateSpeedOfConvergence(out, n, m, a_x, b_x, a_y, b_y, MatrixTrue, MatrixFilter, MatrixRecoveryHor);
			}
		}
	}
	fclose(out);
	return 0;
}

void calculateSpeedOfConvergence(FILE* out, int n, int m, double a_x, double b_x, double a_y, double b_y, double** MatrixTrue, int** MatrixFilter, double** MatrixRecovery){
	double sum = 0.0, p = 0.0, h = (b_x - a_x) / double(n);
	//double h = ((b_x - a_x) * (b_y - a_y)) / double(m * n);
	int count = 0;
	for (int i = 0; i <= m; i++){
		for (int j = 0; j <= n; j++){
			if (MatrixFilter[i][j] == 0){
				count++;
				sum += fabs(MatrixTrue[i][j] - MatrixRecovery[i][j]);
			}
		}
	}
	sum /= double(count);
	fprintf(out, "n = %d\tm = %d\t", n, m);
	if (sum > pow(10, -15)){
		p = log(sum) / log (h);
		fprintf(out, "p = %.10f\n", p);
	} else {
		fprintf(out, "The sum of deltas is less than 10 ^ -15\n");
	}
}

double** furieRecover(int n, int m, double a_x, double b_x, double h_x, double** MatrixTrue, int** MatrixFilter){
	double temp = 0.0, below = 0.0, above = 0.0, t = 0.0, q = 0.0, X_below = 0.0, X_above = 0.0, Y_below = 0.0, Y_above = 0.0, Y_q = 0.0;
	int N = 0, tempInt = 0, count = 0;
	double *c = NULL, *X = NULL, *Y = NULL;
	double **MatrixRecovery = new double*[m + 1];
	for (int s = 0; s <= m; ++s){
		MatrixRecovery[s] = new double[n + 1];
	}
	for (int s = 0; s <= m; ++s){
		N = 0;
		for (int j = 0; j <= n; ++j){
			if (MatrixFilter[s][j] == 1){
				MatrixRecovery[s][j] = MatrixTrue[s][j];
				N++;
			}
		}
		X = new double[N];
		Y = new double[N];
		tempInt = 0;
		for (int j = 0; j <= n; j++){
			if (MatrixFilter[s][j] == 1){
				Y[tempInt] = MatrixTrue[s][j];
				X[tempInt] = a_x + double(tempInt) * h_x;
				tempInt++;
			}
		}
		/*
		printf("X Y:\n");
		for (int i = 0; i < N; i++){
			printf("%f %f\n", X[i], Y[i]);
		}
		printf("\n");
		*/
		c = new double[N];
		double H = (X[N - 1] - X[0]) / double(N - 1);
		for (int p = 1; p < N - 1; ++p){
			c[p] = 0.0;
			for (int j = 1; j < N - 1; j++){
				c[p] += (1.0 / double(N - 1)) * Y[j] * phiFunction(X[0], X[N - 1], p, X[0] + double(j) * H);
			}
		}
		/*
		temp = 0.0;
		for (int j = 0; j < N; j++){
			temp = 0.0;
			for (int p = 1; p < N - 1; ++p){
				temp += c[p] * phiFunction(X[0], X[N - 1], p, X[0] + double(j) * H);
			}
			printf("True: %f ", Y[j]);
			printf("Solved: %f\n", temp);
		}
		printf("\n");
		*/
		below = 0.0;
		above = 0.0; 
		t = 0.0; 
		q = 0.0; 
		X_below = 0.0; 
		X_above = 0.0; 
		Y_below = 0.0; 
		Y_above = 0.0; 
		Y_q = 0.0;
		tempInt = 0;
		count = 0;
		for (int i = 0; i <= n; i++){
			if (MatrixFilter[s][i] == 1){
				below = a_x + double(i) * h_x;
				tempInt++;
			} else {
				X_below = X[tempInt - 1];
				X_above = X[tempInt];
				Y_below = Y[tempInt - 1];
				Y_above = Y[tempInt];
				count = 0;
				for (int k = i + 1; k <= n; ++k){
					if (MatrixFilter[s][k] == 1){
						count = k - i + 1;
						above = a_x + k * h_x;
						break;
					}
				}
				for (int j = 1; j < count; ++j){
					t = below + h_x * j;
					q = X_below + ((X_above - X_below) / double(count)) * j;
					Y_q = 0.0;
					for (int p = 1; p < N - 1; ++p){
						Y_q += c[p] * phiFunction(X[0], X[N - 1], p, q);
					}
					//printf("below = %f above = %f t = %f\nX_below = %f X_above = %f q = %f\nY_below = %f Y_above = %f Y_q = %f\n", below, above, t, X_below, X_above, q, Y_below, Y_above, Y_q);
					MatrixRecovery[s][i + j - 1] = recoverValueScale(below, above, t, Y_below, Y_above, Y_q, X_below, X_above, q);
					//printf("delta[%d][%d] = %f\n", s, i + j - 1, MatrixRecovery[s][i + j - 1] - MatrixTrue[s][i + j - 1]);
				}
				i += (count - 2);
			}
		}
	}
	return MatrixRecovery;
}

double recoverValueScale(double below, double above, double t, double Y_below, double Y_above, double Y_q, double X_below, double X_above, double q){
	return Y_below - ((below - t) * (Y_below - Y_q) * (X_below - X_above)) / ((below - above) * (X_below - q));
}

double phiFunction(double a_x, double b_x, int p, double x){
	return sqrt(2.0) * sin(M_PI * double(p) * (x - a_x) / (b_x - a_x)); 
}

double** splainHorCubedRecover(int n, int m, double a_x, double h_x, double** MatrixTrue, int** MatrixFilter){
	int N = 0, tempInt = 0;
	double *X = NULL, *Y = NULL, *M = NULL, *P = NULL;
	double **C = NULL, **MatrixRecovery;
	MatrixRecovery = new double*[m + 1];
	for (int i = 0; i <= m; ++i){
		MatrixRecovery[i] = new double[n + 1];
	}
	for (int s = 0; s <= m; s++){
		// Make the available points count zero
		N = 0;
		// Figure out how many points are available
		for (int k = 0; k <= n; k++){
			if (MatrixFilter[s][k] == 1){
				MatrixRecovery[s][k] = MatrixTrue[s][k];
				N++;
			}
		}
		// Create X-vector contained available points and Y-vector contained the function value in these available points
		X = new double[N];
		Y = new double[N];
		tempInt = 0;
		for (int k = 0; k <= n; k++){
			if (MatrixFilter[s][k] == 1){
				Y[tempInt] = MatrixTrue[s][k];
				X[tempInt] = a_x + k * h_x;
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
		for (int k = 0; k <= n; k++){
			if (MatrixFilter[s][k] == 1){
				++tempInt;
			} else {
				MatrixRecovery[s][k] = M[tempInt - 1] * pow(X[tempInt] - (a_x + k * h_x), 3) / (6.0 * (X[tempInt] - X[tempInt - 1])) + M[tempInt] * pow((a_x + k * h_x) - X[tempInt - 1], 3) / (6.0 * (X[tempInt] - X[tempInt - 1])) + (Y[tempInt - 1] - (M[tempInt - 1] * pow(X[tempInt] - X[tempInt - 1], 2)) / 6.0) * (X[tempInt] - (a_x + k * h_x)) / (X[tempInt] - X[tempInt - 1]) + (Y[tempInt] - (M[tempInt] * pow(X[tempInt] - X[tempInt - 1], 2)) / 6.0) * (a_x + k * h_x - X[tempInt - 1]) / (X[tempInt] - X[tempInt - 1]);
			}
		}
	}
	return MatrixRecovery;
}

double** splainVerCubedRecover(int n, int m, double a_y, double h_y, double** MatrixTrue, int** MatrixFilter){
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

void initFilterMatrix(int n, int m, int** &MatrixFilter){	
	for (int i = 0; i <= m; i++){
		for (int j = 0; j <= n; j++){
			if (rand() % 2 == 0){
				MatrixFilter[i][j] = 0;
			} else {
				MatrixFilter[i][j] = 1;
			}
			/*
			if (i == 0 || j == 0 || i == m || j == n){
				MatrixFilter[i][j] = 1;
			} else {
				if (rand() % 2 == 0){
					MatrixFilter[i][j] = 0;
				} else {
					MatrixFilter[i][j] = 1;
				}
			}
			*/
		}
	}
}

double f(double a_x, double b_x, double x, double y){
	switch (functionIndex){
		case 0: return (x * x - 1.0) * (y * y - 1.0);
		case 1: return (x * x - 1.0) * (y * y - 1.0) * exp(-1.0 * x * y);
		case 2: return sin(M_PI * y) + sin(2.0 * M_PI * y) + sin(50.0 * M_PI * y);
		case 3: return (x * x - 1.0) * (y * y - 1.0) * atan(y) * exp(x);
		case 4: return (x * x - 1.0) * (y * y - 1.0) * (x * x - 0.01) * (y * y - 0.01) * exp(x * y);
		case 5: return (x * x - 1.0) * (cos(x) - log(y * y + 3.0));
		default: exit (-1); 
	}
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

/*

*/