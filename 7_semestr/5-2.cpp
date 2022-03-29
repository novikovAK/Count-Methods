// to plot y solve print:  plot "plot.txt" using 1:2 title "solve y" with lines
// to plot y true and y solve print:  plot "plot.txt" using 1:2 title "true y" with lines, "plot.txt" using 1:3 title "solve y" with lines

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int sub_main(FILE* out);
void C(double &c, int n);
double Y(int m, int k);
double Yn(int m, int k);
double L(int m);
double A(int i, int j);
double b(int k);
double f(int k);
double Y_true(int k);

int N = 1;
double h = 1.0 / double(N);
double* length = NULL;
bool out_print = true, plot_print = true, vector_checking = true, matrix_checking = true; //file print markers 

int main(){
	FILE* out = fopen("out.txt", "w");
	if (out == NULL){
		printf("Open file out.txt error\n");
		return -1;
	}
	printf("Set N:\n");
	if (scanf("%d", &N) != 1){
		printf("Incorrect N.\n");
		return -1;
	}
	h = 1.0 / double(N);
	length = new double[N - 1];
	if (sub_main(out) != 0){
		printf("Error in the sub_main function.\n");
		return -1;
	}
	fclose(out);
	return 0;
}

int sub_main(FILE* out){
	FILE* plot = fopen("plot.txt", "w");
	if (out == NULL){
		printf("Open file plot.txt error\n");
		return -1;
	}
	double temp = 0.0;
	double c = 0.0;
	double* y_solve = new double[N];
	y_solve[N - 1] = 0.0;
	for (int m = 1; m <= N - 1; ++m){
		y_solve[m - 1] = 0.0;
		length[m - 1] = 0.0;
		for (int k = 1; k <= N - 1; ++k){
			length[m - 1] += h * Y(m, k) * Y(m, k);
		}
		length[m - 1] = sqrt(length[m - 1]);
	}
	double accuracy = 0.0, max_accuracy = 0.0, min_accuracy = 10.0, accuracy_2 = 0.0;
	int max_accuracy_number_m = 0, min_accuracy_number_m = 0;
	if (vector_checking){
		// проверяем, что длины векторов единицы, и вектора ортогональны
		for (int m = 1; m <= N - 1; ++m){
			temp = 0.0;
			for (int k = 1; k <= N - 1; ++k){
				temp += h * Yn(m, k) * Yn(m, k);
			}
			if (temp > max_accuracy){
				max_accuracy = temp;
				max_accuracy_number_m = m;
			}
			if (temp < min_accuracy){
				min_accuracy = temp;
				min_accuracy_number_m = m;
			}
		}
		if (out_print){
			fprintf(out, "report for N = %d, h = %f:\n", N, h);
			fprintf(out, "max (y[%d], y[%d]) = %1.15f\nmin (y[%d], y[%d]) = %1.15f\n", max_accuracy_number_m, max_accuracy_number_m, max_accuracy, min_accuracy_number_m, min_accuracy_number_m, min_accuracy);
		}
		max_accuracy = 0.0;
		for (int m = 1; m <= N - 1; ++m){
			for (int s = m + 1; s <= N - 1; ++s){
				temp = 0.0;
				for (int k = 1; k <= N - 1; ++k){
					temp += h * Yn(m, k) * Yn(s, k);
				}
				if (temp > max_accuracy){
					max_accuracy = temp;
					max_accuracy_number_m = m;
					min_accuracy_number_m = s;
				}
			}
		}
		if (out_print){
			fprintf(out, "max (y[%d], y[%d]) = %.3e\n", max_accuracy_number_m, min_accuracy_number_m, max_accuracy);
		}
	}
	if (matrix_checking){
		max_accuracy = 0.0;
		for (int m = 1; m <= N - 1; ++m){
			accuracy = 0.0;
			accuracy_2 = 0.0;
			for (int k = 1; k <= N - 1; ++k){
				temp = 0.0;
				for (int s = 1; s <= N - 1; ++s){
					temp += A(k, s) * Yn(m, s);
				}
				accuracy += h * pow(temp - L(m) * Yn(m, k), 2);
				accuracy_2 += h * pow(L(m) * Yn(m, k), 2);
			}
			accuracy = sqrt(accuracy / accuracy_2);
			if (accuracy > max_accuracy){
				max_accuracy = accuracy;
				max_accuracy_number_m = m;
			}
		}
		if (out_print){
			fprintf(out, "max || Ay - Ly ||h / || Ly ||h = %1.3e for %d egen function\n", max_accuracy, max_accuracy_number_m);
		}
	}
	//Решение задачи
	for (int n = 1; n <= N - 1; ++n){
		C(c, n);
		for (int k = 1; k <= N - 1; ++k){
			y_solve[k - 1] += c * Yn(n, k);
		}
	}
	//Ищем невязку ||Ау - b||
	accuracy = 0.0;
	for (int i = 1; i <= N - 1; ++i){
		temp = 0.0;
		for (int j = 1; j <= N - 1; ++j){
			temp += A(i, j) * y_solve[j - 1];
		}
		accuracy += pow(temp - b(i), 2);
	}
	if (out_print){
		fprintf(out, "||Ay - b||2 = %1.3e\n", sqrt(accuracy));
	}
	// for degree
	double p_deg = 0.0;
	accuracy = h * pow(Y_true(0) - y_solve[0] + f(0) * h * h * 0.5, 2);
	if (plot_print){
		fprintf(plot, "%.10f %.10f %.10f\n", 0.0, Y_true(0), y_solve[0] - f(0) * h * h * 0.5);
	}
	for (int k = 1; k <= N; ++k){
		accuracy += h * pow(Y_true(k) - y_solve[k - 1], 2);
		if (plot_print){
			fprintf(plot, "%.10f %.10f %.10f\n", double(k) * h, Y_true(k), y_solve[k - 1]);
		}
	}
	accuracy = sqrt(accuracy);
	if (accuracy > pow(10, -20)){
		p_deg = log(accuracy) / log(h);
		if (out_print){
			fprintf(out, "||true y - solve y||h = %1.7f = h in [ %1.7f ]\n", accuracy, p_deg);
		}
	}
	fclose(plot);
	return 0;
}

void C(double &c, int n){
	c = 0.0;
	for (int i = 1; i <= N - 1; ++i){
		c += h * b(i) * Yn(n, i);
	}
	c = c / L(n);
}

double Y(int m, int k){
	return sin((M_PI * k * (2 * m - 1)) / (2 * N - 1)) - tan((M_PI * N * (2 * m - 1)) / (2 * N - 1)) * cos((M_PI * k * (2 * m - 1)) / (2 * N - 1));
}

double Yn(int m, int k){
	return Y(m, k) / length[m - 1];
}

double L(int m){
	return (-4.0 / (h * h)) * pow(sin((M_PI * (2 * m - 1)) / (4 * N - 2)), 2);
}

double A(int i, int j){
	if (i == j){
		if (i == 1){
			return -1.0 / (h * h);
		} else {
			return -2.0 / (h * h);
		}
	} else if (i == j + 1 || j == i + 1){
		return 1.0 / (h * h);
	} else {
		return 0.0;
	}
}

double b(int k){
	if (k == 1){
		return f(1) + 0.5 * f(0);
	} else {
		return f(k);
	}
}

double f(int k){
	return ((-9.0 * M_PI * M_PI) / 4.0) * cos(1.5 * M_PI * double(k) * h);
}

double Y_true(int k){
	return cos(1.5 * M_PI * double(k) * h);
}
