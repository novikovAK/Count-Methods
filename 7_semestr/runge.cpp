// plot h: plot "2.txt" using 1:2 with lines title "h(x)"
// plot solve and true M = 1: plot "2.txt" using 1:3 with lines title "y(x) true", "2.txt" using 1:4 with lines title "y(x) solve"
// plot solve and true M = 2: plot "2.txt" using 1:3 with lines title "y(x) true 1", "2.txt" using 1:4 with lines title "y(x) solve 1", "2.txt" using 1:5 with lines title "y(x) true 2", "2.txt" using 1:6 with lines title "y(x) solve 2"  

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void step(double &x_curr, double &h, double* &f_out, double* &y_curr, double* &y_new, double* &y_new_2, double* &y_temp, double* &k1, double* &k2);
void step_auto(FILE* f_1, FILE* f_2, double &x_curr, double &h, double &NORM, double* &f_out, double* &y_curr, double* &y, double* &y_new, double* &y_new_2, double* &y_temp, double* &k1, double* &k2);
void f_function(double* &f_out, double x_curr, double* y_curr);
void y_true(double* &y, double &x_curr);
void norma(double &NORM, double h, double* y_new, double* y_new_2);
void print_solution(FILE* f_1, FILE* f_2, double &x_curr, double &h, double* &y_curr, double* &y);

double A = 1.0, EPS_m = pow(10, -12), EPS_M = pow(10, -10), H_MIN = pow(10, -16), alpha2 = 1.0, beta21 = 1.0, p1 = 0.5, p2 = 0.5, x = 2 * M_PI, degree = 0.0, min_degree = 10.0, max_degree = 0.0, min_degree_x = 0.0, max_degree_x = 0.0;
int M = 1, min_degree_i = 0, max_degree_i = 0;
bool index_h = false, usual_eiler = false, autostep_eiler = false, usual_runge = false, autostep_runge = true;

int main (){
	FILE* f_1 = fopen("1.txt", "w");
	FILE* f_2 = fopen("2.txt", "w");
	int N = 100;
	double h = x / double(N), x_curr = 0.0, NORM = 0.0;
	double* f_out = new double[M];
	double* y_curr = new double[M];
	double* y_new = new double[M];
	double* y_new_2 = new double[M];
	double* y_temp = new double[M];
	double* y = new double[M];
	double* y_0 = new double[M];
	double* k1 = new double[M];
	double* k2 = new double[M];
	y_0[0] = 1.0;
	y_0[1] = 0.0;
	for (int i = 0; i < M; ++i){
		f_out[i] = 0.0;
		y_curr[i] = y_0[i];
		y_new[i] = 0.0;
		y_new_2[i] = 0.0;
		y_temp[i] = 0.0;
		y[i] = 0.0;
	}
	print_solution(f_1, f_2, x_curr, h, y_curr, y);
	while (x_curr < x){
		step_auto(f_1, f_2, x_curr, h, NORM, f_out, y_curr, y, y_new, y_new_2, y_temp, k1, k2);
	}
	fclose (f_1);
	fclose (f_2);
	return 0;
}

void step(double &x_curr, double &h, double* &f_out, double* &y_curr, double* &y_new, double* &y_new_2, double* &y_temp, double* &k1, double* &k2){
	if (autostep_runge){
		//case 1
		f_function(f_out, x_curr, y_curr);
		for (int i = 0; i < M; ++i){
			k1[i] = h * f_out[i];
			y_new[i] = y_curr[i] + beta21 * k1[i];
		}
		f_function(f_out, x_curr + alpha2 * h, y_new);
		for (int i = 0; i < M; ++i){
			k2[i] = h * f_out[i];
			y_new[i] = y_curr[i] + p1 * k1[i] + p2 * k2[i];
		}
		//case 2 h = h / 2
		h *= 0.5;
		f_function(f_out, x_curr, y_curr);
		for (int i = 0; i < M; ++i){
			k1[i] = h * f_out[i];
			y_temp[i] = y_curr[i] + beta21 * k1[i];
		}
		f_function(f_out, x_curr + alpha2 * h, y_temp);
		for (int i = 0; i < M; ++i){
			k2[i] = h * f_out[i];
			y_temp[i] = y_curr[i] + p1 * k1[i] + p2 * k2[i];
		}
		x_curr += h;
		f_function(f_out, x_curr, y_temp);
		for (int i = 0; i < M; ++i){
			k1[i] = h * f_out[i];
			y_new_2[i] = y_temp[i] + beta21 * k1[i];
		}
		f_function(f_out, x_curr + alpha2 * h, y_new_2);
		for (int i = 0; i < M; ++i){
			k2[i] = h * f_out[i];
			y_new_2[i] = y_temp[i] + p1 * k1[i] + p2 * k2[i];
		}
		//restore initial values
		x_curr -= h;
		h *= 2.0;
	}
	if (usual_runge){
		f_function(f_out, x_curr, y_curr);
		for (int i = 0; i < M; ++i){
			k1[i] = h * f_out[i];
			y_new[i] = y_curr[i] + beta21 * k1[i];
		}
		f_function(f_out, x_curr + alpha2 * h, y_new);
		for (int i = 0; i < M; ++i){
			k2[i] = h * f_out[i];
			y_new[i] = y_curr[i] + p1 * k1[i] + p2 * k2[i];
		}
	}
	if (autostep_eiler){
		//case 1
		f_function(f_out, x_curr, y_curr);
		for (int i = 0; i < M; ++i){
			y_new[i] = y_curr[i] + h * f_out[i];
		}
		//case 2
		h *= 0.5;
		for (int i = 0; i < M; ++i){
			y_temp[i] = y_curr[i] + h * f_out[i];
		}
		x_curr += h;
		f_function(f_out, x_curr, y_temp);
		for (int i = 0; i < M; ++i){
			y_new_2[i] = y_temp[i] + h * f_out[i];
		}
		x_curr -= h;
		h *= 2.0;
	}
	if (usual_eiler){
		f_function(f_out, x_curr, y_curr);
		for (int i = 0; i < M; ++i){
			y_new[i] = y_curr[i] + h * f_out[i];
		}
	}
}

void step_auto(FILE* f_1, FILE* f_2, double &x_curr, double &h, double &NORM, double* &f_out, double* &y_curr, double* &y, double* &y_new, double* &y_new_2, double* &y_temp, double* &k1, double* &k2){
	//Для постоянного шага + уменьшение последнего шага в 2 раза
	if (usual_runge || usual_eiler){
		if (x_curr + h >= x && !index_h){
			h *= 0.5;
			index_h = true;
			if (h < H_MIN){
				exit (0);
			}
		}
		step(x_curr, h, f_out, y_curr, y_new, y_new_2, y_temp, k1, k2);
		x_curr += h;
		for (int i = 0; i < M; ++i){
			y_curr[i] = y_new[i];
		}
		print_solution(f_1, f_2, x_curr, h, y_curr, y);
	}
	
	//Автошаг
	if (autostep_eiler || autostep_runge){
		step(x_curr, h, f_out, y_curr, y_new, y_new_2, y_temp, k1, k2);
		norma(NORM, h, y_new, y_new_2);
		if (NORM > EPS_m && NORM < EPS_M){
			x_curr += h;
			for (int i = 0; i < M; ++i){
				y_curr[i] = y_new[i];
			}
			print_solution(f_1, f_2, x_curr, h, y_curr, y);
		} else if (NORM < EPS_m){
			x_curr += h;
			for (int i = 0; i < M; ++i){
				y_curr[i] = y_new[i];
			}
			print_solution(f_1, f_2, x_curr, h, y_curr, y);
			h *= 2;
		} else if (NORM > EPS_M){
			while (NORM > EPS_M){
				h *= 0.5;
				if (h < H_MIN){
					printf("Too small h\n");
					exit (-1);
				}
				step(x_curr, h, f_out, y_curr, y_new, y_new_2, y_temp, k1, k2);
				norma(NORM, h, y_new, y_new_2);
			}
			x_curr += h;
			for (int i = 0; i < M; ++i){
				y_curr[i] = y_new[i];
			}
			print_solution(f_1, f_2, x_curr, h, y_curr, y);
		} else {
			printf("Unknown norma value\n");
			exit (-1);
		}
	}
}

void f_function(double* &f_out, double x_curr, double* y_curr){
	//f_out[0] = -1.0 * A * y_curr[1] + x_curr * 0.0;
	//f_out[1] = A * y_curr[0];
	f_out[0] = 2.0  * x_curr + 1.0 + y_curr[0] * 0.0;
}

void y_true(double* &y, double &x_curr){
	//y[0] = cos(A * x_curr);
	//y[1] = sin(A * x_curr);
	y[0] = x_curr * x_curr + x_curr + 1.0;
}

void norma(double &NORM, double h, double* y_new, double* y_new_2){
	NORM = 0.0;
	for (int i = 0; i < M; ++i){
		NORM += pow(y_new[i] - y_new_2[i], 2);
	}
	NORM = sqrt(NORM);
}

void print_solution(FILE* f_1, FILE* f_2, double &x_curr, double &h, double* &y_curr, double* &y){
	fprintf(f_1, "x = %.3f  h = %.7f  ", x_curr, h);
	y_true(y, x_curr);
	fprintf(f_2, "%.6f %.10f ", x_curr, h);
	for (int i = 0; i < M; ++i){
		if (fabs(y[i] - y_curr[i]) > pow(10, -15)){
			degree = log(fabs(y[i] - y_curr[i])) / log(h);
			if (degree > max_degree){
				max_degree = degree;
				max_degree_x = x_curr;
				max_degree_i = i;
			}
			if (degree < min_degree){
				min_degree = degree;
				min_degree_x = x_curr;
				min_degree_i = i;
			}
		}
		fprintf(f_1, "%d:  TRUE = %.5f  SOLVE = %.5f  DEL = %.10f  ", i, y[i], y_curr[i], fabs(y[i] - y_curr[i]));
	}
	fprintf(f_1, "\n");
	for (int i = 0; i < M; ++i){
		fprintf(f_2, "%.10f %.10f ", y[i], y_curr[i]);
	}
	fprintf(f_2, "\n");
	if (fabs(x_curr - x) < pow(10, -15)){
		printf("%.6f[x = %.3f, i = %d] < degree < %.6f[x = %.3f, i = %d]\n", min_degree, min_degree_x, min_degree_i, max_degree, max_degree_x, max_degree_i);
		for (int i = 0; i < M; ++i){
			if (fabs(y[i] - y_curr[i]) > pow(10, -15)){
				degree = log(fabs(y[i] - y_curr[i])) / log(h);
				printf("x = %.3f DEL = %.10f = [h = %.10f] ^ [%.4f]\n", x_curr, fabs(y[i] - y_curr[i]), h, degree);
			}
		}
	}
}
