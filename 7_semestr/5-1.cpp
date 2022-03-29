#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double A(int i, int j);
double Y(int m, int k);
double L(int m);

int N = 10;
double h = 0.001 / double(N);

int main(){
	FILE* out = fopen("out.txt", "w");
	if (out == NULL){
		printf("Open file out.txt error\n");
		return -1;
	}
	int max_norma_number = 0, max_scalar_number_1 = 0, max_scalar_number_2 = 0, max_length_number = 0, min_length_number = 0;
	double v = 0.0, u = 0.0, norma_numerator = 0.0, norma_denominator = 0.0, norma = 0.0, max_norma = 0.0, scalar = 0.0, max_scalar = 0.0, length = 0.0, max_length = 0.0, min_length = 1000000;
	fprintf(out, "N = %d  h = %f\n\n", N, h);
	for (int m = 0; m <= N - 2; ++m){
		for (int s = m + 1; s <= N - 2; ++s){
			scalar = 0.0;
			for (int k = 1; k <= N - 1; ++k){
				scalar += h * Y(m, k) * Y(s, k);
			}
			if (scalar > max_scalar){
				max_scalar = scalar;
				max_scalar_number_1 = m;
				max_scalar_number_2 = s;
			}
		}
		length = 0.0;
		for (int k = 1; k <= N - 1; ++k){
			length += h * Y(m, k) * Y(m, k);
		}
		length = sqrt(length);
		if (length > max_length){
			max_length = length;
			max_length_number = m;
		}
		if (length < min_length){
			min_length = length;
			min_length_number = m;
		}
		for (int k = 1; k <= N - 1; ++k){
			v = 0.0;
			for (int j = 1; j <= N - 1; ++j){
				v += A(k, j) * Y(m, j);
			}
			u = L(m) * Y(m, k);
			norma_numerator += (v - u) * (v - u);
			norma_denominator += v * v;
		}
		norma = sqrt(norma_numerator / norma_denominator);
		if (norma > max_norma){
			max_norma = norma;
			max_norma_number = m;
		}
	}
	fprintf(out, "max norma: y[ %d ]: %2.3e\n", max_norma_number, max_norma);
	fprintf(out, "max scalar: (y[ %d ], y[ %d ]) = %2.3e\n", max_scalar_number_1, max_scalar_number_2, max_scalar);
	fprintf(out, "min length: |y[ %d ]| = %2.15f\nmax length: |y[ %d ]| = %2.15f", min_length_number, min_length, max_length_number, max_length);
	fclose(out);
	return 0;
}

double Y(int m, int k){
	return sin((M_PI * k * (2 * m + 1)) / (2 * N - 1)) - tan((M_PI * N * (2 * m + 1)) / (2 * N - 1)) * cos((M_PI * k * (2 * m + 1)) / (2 * N - 1));
}

double L(int m){
	return (-4.0 / (h * h)) * pow(sin((M_PI * (2 * m + 1)) / (4 * N - 2)), 2);
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


/*

N = 1000  h = 0.000001

max norma: y[ 0 ]: 1.293e-10
max scalar: (y[ 0 ], y[ 1 ]) = 8.598e-14
min length: |y[ 998 ]| = 0.022355116513222
max length: |y[ 0 ]| = 28.449154457459617

*/
