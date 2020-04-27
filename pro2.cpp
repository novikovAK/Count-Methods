#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

double recoverValueScale(double below, double above, double t, double Y_below, double Y_above, double Y_q, double X_below, double X_above, double q);
double phiFunction(double a, double b, int m, double x);
double f(double a, double b, double x);

int main(){
	srand(time(NULL));
	int n = 300;
	double a = 2.0, b = 6.0;
	double h = (b - a) / double(n);
	/*
	printf("x f(x):\n");
	for (int i = 0; i <= n; ++i){
		printf("%f %f\n", a + double(i) * h, f(a, b, a + double(i) * h));
	}
	printf("\n");
	*/
	double *c = NULL, *R = new double[n + 1];
	int *I = new int[n + 1];
	I[0] = 1;
	I[n] = 1;
	for (int i = 1; i < n; ++i){
		if (rand() % 10 == 0){
			I[i] = 0;
		} else {
			I[i] = 1;
		}
	}
	/*
	printf("I: ");
	for (int i = 0; i <= n; ++i){
		printf("%d ", I[i]);	
	}
	printf("\n\n");
	*/
	int N = 0, tempInt = 0;
	double *X = NULL, *Y = NULL;
	for (int s = 0; s <= n; s++){
		if (I[s] == 1){
			R[s] = f(a, b, a + double(s) * h);
			N++;
		}
	}
	X = new double[N];
	Y = new double[N];
	tempInt = 0;
	for (int s = 0; s <= n; s++){
		if (I[s] == 1){
			Y[tempInt] = f(a, b, a + double(s) * h);
			X[tempInt] = a + double(tempInt) * h;
			tempInt++;
		}
	}
	/*
	printf("X Y:\n");
	for (int s = 0; s < N; s++){
		printf("%f %f\n", X[s], Y[s]);
	}
	printf("\n");
	*/
	c = new double[N];
	double H = (X[N - 1] - X[0]) / double(N - 1);
	for (int m = 1; m < N - 1; ++m){
		c[m] = 0.0;
		for (int j = 1; j < N - 1; j++){
			c[m] += (1.0 / double(N - 1)) * Y[j] * phiFunction(X[0], X[N - 1], m, X[0] + double(j) * H);
		}
	}
	/*
	double temp = 0.0;
	for (int j = 0; j < N; j++){
		temp = 0.0;
		for (int m = 1; m < N - 1; ++m){
			temp += c[m] * phiFunction(X[0], X[N - 1], m, X[0] + double(j) * H);
		}
		printf("True: %f ", Y[j]);
		printf("Solved: %f\n", temp);
	}
	printf("\n");
	*/
	double below = 0.0, above = 0.0, t = 0.0, q = 0.0, X_below = 0.0, X_above = 0.0, Y_below = 0.0, Y_above = 0.0, Y_q = 0.0;
	tempInt = 0;
	int count = 0;
	for (int i = 0; i <= n; i++){
		if (I[i] == 1){
			below = a + double(i) * h;
			tempInt++;
		} else {
			X_below = X[tempInt - 1];
			X_above = X[tempInt];
			Y_below = Y[tempInt - 1];
			Y_above = Y[tempInt];
			count = 0;
			for (int k = i + 1; k <= n; ++k){
				if (I[k] == 1){
					count = k - i + 1;
					above = a + k * h;
					break;
				}
			}
			for (int j = 1; j < count; ++j){
				t = below + h * j;
				q = X_below + ((X_above - X_below) / double(count)) * j;
				Y_q = 0.0;
				for (int m = 1; m < N - 1; ++m){
					Y_q += c[m] * phiFunction(X[0], X[N - 1], m, q);
				}
				//printf("below = %f above = %f t = %f\nX_below = %f X_above = %f q = %f\nY_below = %f Y_above = %f Y_q = %f\n", below, above, t, X_below, X_above, q, Y_below, Y_above, Y_q);
				R[i + j - 1] = recoverValueScale(below, above, t, Y_below, Y_above, Y_q, X_below, X_above, q);
				//printf("R[%d] = %f\n\n", i + j - 1, R[i + j - 1]);
			}
			i += (count - 2);
		}
	}
	//printf("\n");
	for (int i = 0; i <= n; ++i){
		if (I[i] == 0){
			printf("del[%d]: %f\n", i, f(a, b, a + i * h) - R[i]);
		}
	}
	return 0;
}

double recoverValueScale(double below, double above, double t, double Y_below, double Y_above, double Y_q, double X_below, double X_above, double q){
	return Y_below - ((below - t) * (Y_below - Y_q) * (X_below - X_above)) / ((below - above) * (X_below - q));
}

double phiFunction(double a, double b, int m, double x){
	return sqrt(2.0) * sin(M_PI * double(m) * (x - a) / (b - a)); 
}

double f(double a, double b, double x){
	return (x - a) * (x - b) * cos(x);
}