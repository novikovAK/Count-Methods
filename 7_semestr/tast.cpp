#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

void f(double** &d2r, double &G, double* &m, double** &r);
double length_rij(double** &r, int i, int j);

int main(){
	srand( time( NULL ) );
	double** r = new double*[3];
	for (int i = 0; i < 3; ++i){
		r[i] = new double[3];
	}
	double** d2r = new double*[3];
	for (int i = 0; i < 3; ++i){
		d2r[i] = new double[3];
	}
	double* m = new double[3];
	double G = 6.67 * pow(10, -11);
	
	m[0] = 1.0;
	m[1] = 0.5;
	m[2] = 0.75;
	for (int i = 0; i < 3; ++i){
		for (int j = 0; j < 3; ++j){
			r[i][j] = rand() % 10 + 1;
			printf("%f\t", r[i][j]);
		}
		printf("\n");
	}
	f(d2r, G, m, r);
	for (int i = 0; i < 3; ++i){
		for (int j = 0; j < 3; ++j){
			printf("%f\t", d2r[i][j]);
		}
		printf("\n");
	}
	return 0;
}

void f(double** &d2r, double &G, double* &m, double** &r){
	for (int k = 0; k < 3; ++k){
		for (int i = 0; i < 3; ++i){
			d2r[k][i] = 0.0;
		}
	}
	for (int k = 0; k < 3; ++k){
		d2r[k][0] = ((-1) * G * m[1] * (r[k][0] - r[k][1])) / pow(length_rij(r, 0, 1), 3) - ((-1) * G * m[2] * (r[k][0] - r[k][2])) / pow(length_rij(r, 0, 2), 3);
		d2r[k][1] = ((-1) * G * m[0] * (r[k][1] - r[k][0])) / pow(length_rij(r, 1, 0), 3) - ((-1) * G * m[2] * (r[k][1] - r[k][2])) / pow(length_rij(r, 1, 2), 3);
		d2r[k][2] = ((-1) * G * m[0] * (r[k][2] - r[k][0])) / pow(length_rij(r, 2, 0), 3) - ((-1) * G * m[1] * (r[k][2] - r[k][1])) / pow(length_rij(r, 2, 1), 3);
	}
}

double length_rij(double** &r, int i, int j){
	double length = sqrt((r[0][i] - r[0][j]) * (r[0][i] - r[0][j]) + (r[1][i] - r[1][j]) * (r[1][i] - r[1][j]) + (r[2][i] - r[2][j])*(r[2][i] - r[2][j]));
	printf("length = %f\n", length);
	return length;
}
