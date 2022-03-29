// to plot y solve print:  plot "plot.txt" using 1:2 title "solve y" with lines
// to plot y true and y solve print:  plot "plot.txt" using 1:2 title "true y" with lines, "plot.txt" using 1:3 title "solve y" with lines

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void solve(FILE* out, FILE* plot, int N, double h, double* alpha, double* beta, double* a, double* b, double* c, double* y, double norma, double p);
double f_f(double h, int N, int k);
double k_f(double h, int k, double half);
double p_f(double h, int k);
double y_true(double h, int k);

bool out_print = true, plot_print = true; //file print markers 

int main(){
	int N = -1;
	double h = 0.0, norma = 0.0, p = 0.0;
	FILE* out = fopen("out.txt", "w");
	if (out == NULL){
		printf("Incorrect file out.txt\n");
		return 0;
	}
	FILE* plot = fopen("plot.txt", "w");
	if (plot == NULL){
		printf("Incorrect file plot.txt\n");
		return 0;
	}
	double *alpha = NULL;
	double *beta = NULL;
	double *a = NULL;
	double *b = NULL;
	double *c = NULL;
	double *y = NULL;
	printf("Set N:\n");
	if (scanf("%d", &N) != 1){
		printf("Incorrect N\n");
		return 0;
	}
	solve(out, plot, N, h, alpha, beta, a, b, c, y, norma, p);
	fclose(out);
	fclose(plot);
	return 0;
}

void solve(FILE* out, FILE* plot, int N, double h, double* alpha, double* beta, double* a, double* b, double* c, double* y, double norma, double p){
	h = 1.0 / double(N);
	if (out_print){
		fprintf(out, "N = %d  h = %f\n\n", N, h);
	}
	alpha = new double[N + 1];
	beta = new double[N + 1];
	a = new double[N + 1];
	b = new double[N + 1];
	c = new double[N + 1];
	y = new double[N + 1];
	for (int j = 0; j <= N; ++j){
		if (j == 0){
			a[j] = 0.0;
			b[j] = -1.0 / h;
			c[j] = -1.0 * (h * p_f(h, j)) / (2 * k_f(h, 0, 0)) - 1.0 / h;
			if (out_print){
				fprintf(out, "b[%d] = %.5f  c[%d] = %.5f\n", j, b[j], j, c[j]);
			}
		} else if (j == N){
			a[j] = 0.0;
			b[j] = 0.0;
			c[j] = 1.0;
			if (out_print){
				fprintf(out, "a[%d] = %.5f  c[%d] = %.5f\n", j, a[j], j, c[j]);
			}
		} else {
			a[j] = (1.0 / (h * h)) * k_f(h, j, -0.5);
			b[j] = (1.0 / (h * h)) * k_f(h, j, 0.5);
			c[j] = (1.0 / (h * h)) * (k_f(h, j, 0.5) + k_f(h, j, -0.5)) + p_f(h, j);
			if (out_print){
				fprintf(out, "a[%d] = %.5f  b[%d] = %.5f  c[%d] = %.5f\n", j, a[j], j, b[j], j, c[j]);
			}
		}
	}
	if (out_print){
		fprintf(out, "\n");
	}
	for (int k = 1; k < N; ++k){
		if (k == 1){
			alpha[k] = b[0] / c[0];
			beta[k] = f_f(h, N, 0) / c[0];
			if (out_print){
				fprintf(out, "alpha[%d] = %.5f  beta[%d] = %.5f\n", k, alpha[k], k,  beta[k]);
			}
			alpha[k + 1] = (b[k]) / (c[k] - a[k] * alpha[k]);
			beta[k + 1] = (f_f(h, N, k) + a[k] * beta[k]) / (c[k] - a[k] * alpha[k]);
			if (out_print){
				fprintf(out, "alpha[%d] = %.5f  beta[%d] = %.5f\n", k + 1, alpha[k + 1], k + 1, beta[k + 1]);
			}
		} else {
			alpha[k + 1] = (b[k]) / (c[k] - a[k] * alpha[k]);
			beta[k + 1] = (f_f(h, N, k) + a[k] * beta[k]) / (c[k] - a[k] * alpha[k]);
			if (out_print){
				fprintf(out, "alpha[%d] = %.5f  beta[%d] = %.5f\n", k + 1, alpha[k + 1], k + 1, beta[k + 1]);
			}
		}
	}
	if (out_print){
		fprintf(out, "\n");
	}
	for (int k = N; k >= 0; --k){
		if (k == N){
			y[k] = (f_f(h, N, k) + a[k] * beta[k]) / (c[k] - a[k] * alpha[k]);
		} else {
			y[k] = alpha[k + 1] * y[k + 1] + beta[k + 1];
		}
		if (out_print){
			fprintf(out, "y[%d] = %.5f\n", k, y[k]);
		}
	}
	if (plot_print){
		for (int k = 0; k <= N; ++k){
			//fprintf(plot, "%.10f %.10f\n", double(k) * h, y[k]); //for unknown true solve
			fprintf(plot, "%.10f %.10f %.10f\n", double(k) * h, y_true(h, k), y[k]); //for known true solve
		}
	}
	if (out_print){
		fprintf(out, "\n");
	}
}

double f_f(double h, int N, int k){
	if (k == 0){
 		return ((h * 9.0 * M_PI * M_PI) / (8.0 * k_f(h, 0, 0))) * cos(1.5 * M_PI * double(k) * h);
	} else if (k == N){
		return 0.0;
	} else {
		return ((-9.0 * M_PI * M_PI) / 4.0) * cos(1.5 * M_PI * double(k) * h);
	}
}

double k_f(double h, int k, double half){
	/*
	if (half == 0.5 || half == -0.5 || half == 0.0){
		return 1 + (double(k) + half) * h;
	} else {
		printf("Error in the k function\n");
		exit (-1);
	}
	return -1;
	*/
	return -1.0;
}

double p_f(double h, int k){
	//return 1.0 + 0.5 * sin(M_PI * double(k) * h);
	return 0.0;
}

double y_true(double h, int k){
	return cos(1.5 * M_PI * double(k) * h);
}
