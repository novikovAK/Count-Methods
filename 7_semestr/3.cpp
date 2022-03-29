 
// y_solve(x) and y_true(x)    plot "plot.txt" using 1:2 with lines title "solve Y", "plot_true.txt" using 1:2 with lines title "true Y"
// h(x)   plot "plot.txt" using 1:2 with lines title "h"
// y_solve(x) and y_true(x) from two eq   plot "plot.txt" using 1:3 with lines title "solve Y 1", "plot.txt" using 1:4 with lines title "solve Y 2", "plot_true.txt" using 1:2 with lines title "true Y 1", "plot_true.txt" using 1:3 with lines title "true Y 2"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void step_auto(FILE* f, FILE* p, double* &Y_new, double* &Y_new_h, double* &Y_cur, double* &Y_out, double* Y_temp, double* K_1, double* K_2, double* &F_out,  double &x_cur, double &h, int M, double &norma, double p_1, double p_2, double alpha_2, double beta_2_1);
void step(double* &Y_new, double* &Y_new_h, double* &Y_cur, double* Y_temp, double* K_1, double* K_2, double* &F_out, double &x_cur, double &h, int M, double p_1, double p_2, double alpha_2, double beta_2_1);
void F(double* &F_out, double x_cur, double* Y_cur, int M);
void Y(double* &Y_out, double x_cur, int M);
void NORMA(double &norma, double* Y_out, double* Y_new, int M);
void print(FILE* f, FILE* p, double &x_cur, double &h, double* &Y_cur, double* &Y_out, int M);

const double EPS_m = pow(10, -2), EPS_M = pow(10, -1), H_MIN = pow(10, -15);
double A = 1.0;
double p_deg = 0.0, deg_min = 10.0, deg_min_x = 0.0;
double X = 2.0 * A * M_PI;
//double X = 5.0;
int deg_min_i = 0;

int main(){
    FILE* f = fopen("out.txt", "w");
	FILE* p = fopen("plot.txt", "w");
    FILE* p_true = fopen("plot_true.txt", "w");
    int N = 10000, M = 2;
    double s = 0.0;
    double h = X / double(N), x_cur = 0.0, p_1 = 0.5, p_2 = 0.5, norma = 0.0, alpha_2 = 1.0, beta_2_1 = 1.0;
    double* Y_new = new double[M];
    double* Y_new_h = new double[M];
    double* Y_cur = new double[M];
    double* Y_0 = new double[M];
    double* Y_out = new double[M];
    double* Y_temp = new double[M];
    double* F_out = new double[M];
    double* K_1 = new double[M];
    double* K_2 = new double[M];
    Y_0[0] = 0.0;
	Y_0[1] = A;
    for (int i = 0; i < M; ++i){
        Y_cur[i] = Y_0[i];
        Y_out[i] = 0.0;
        F_out[i] = 0.0;
    }
    for (int k = 0; k <= 100; ++k){
        s = double(k) * (X / 100.0);
        fprintf(p_true, "%.6f ", s);
        Y(Y_out, s, M);
        for (int i = 0; i < M; ++i){
            fprintf(p_true, "%.6f ", Y_out[i]);
        }
        fprintf(p_true, "\n");
    }
    Y(Y_out, x_cur, M);
    fprintf(f, "x = %1.3f  h = %1.3e  ", x_cur, h);
	fprintf(p, "%2.6f %2.6f ", x_cur, h);
    for (int i = 0; i < M; ++i){
		if (fabs(Y_out[i] - Y_cur[i]) < pow(10, -15)){
			p_deg = 10;
		} else {
			p_deg = (log(fabs(Y_out[i] - Y_cur[i])) / log(h));
		}
		if (p_deg < deg_min){
			deg_min = p_deg;
			deg_min_i = i;
			deg_min_x = x_cur;
		}
        fprintf(f, "y[%d] = %1.3e  y(%1.3f) = %1.3e  del[%d] = %1.3e  p = %1.3f  ", i, Y_cur[i], x_cur, Y_out[i], i, Y_out[i] - Y_cur[i], p_deg);
		//fprintf(f, "del[%d] = %1.5e  ", i, Y_out[i] - Y_cur[i]);
		fprintf(p, "%.6f ", Y_cur[i]);
    }
    fprintf(f, "\n");
	fprintf(p, "\n");
    while (x_cur < X){
        step_auto(f, p, Y_new, Y_new_h, Y_cur, Y_out, Y_temp, K_1, K_2, F_out, x_cur, h, M, norma, p_1, p_2, alpha_2, beta_2_1);
		//step(Y_new, Y_new_h, Y_cur, Y_temp, K_1, K_2, F_out, x_cur, h, M, p_1, p_2, alpha_2, beta_2_1);
		//x_cur += h;
        //Y(Y_out, x_cur, M);
        //print(f, p, x_cur, h, Y_new, Y_out, M);
        //for (int i = 0; i < M; ++i){
        //    Y_cur[i] = Y_new[i];
		//}
    }
	printf("min p(%2.5f)[%d] = %2.5f\n", deg_min_x, deg_min_i, deg_min);
    fclose(f);
    fclose(p);
    fclose(p_true);
    return 0;
}

void step_auto(FILE* f, FILE* p, double* &Y_new, double* &Y_new_h, double* &Y_cur, double* &Y_out, double* Y_temp, double* K_1, double* K_2, double* &F_out,  double &x_cur, double &h, int M, double &norma, double p_1, double p_2, double alpha_2, double beta_2_1){
    step(Y_new, Y_new_h, Y_cur, Y_temp, K_1, K_2, F_out, x_cur, h, M, p_1, p_2, alpha_2, beta_2_1);
    NORMA(norma, Y_new, Y_new_h, M);
    if (norma < EPS_M && norma > EPS_m){
        x_cur += h;
		/*
        if (x_cur > X){
            x_cur = X;
        }
        */
        Y(Y_out, x_cur, M);
        for (int i = 0; i < M; ++i){
            Y_cur[i] = Y_new[i];
        }
        print(f, p, x_cur, h, Y_cur, Y_out, M);
        //printf("CASE GOOD: h = %2.20f   x = %2.5f\n", h, x_cur);
    } else if (norma < EPS_m){
        x_cur += h;
		h *= 2.0;
		/*
        if (x_cur > X){
            x_cur = X;
        }
        */
        Y(Y_out, x_cur, M);
        //printf("CASE INCREASE: h = %2.20f  x = %2.5f\n", h, x_cur);
        for (int i = 0; i < M; ++i){
            Y_cur[i] = Y_new[i];
        }
        print(f, p, x_cur, h, Y_cur, Y_out, M);
    } 
    else if (norma > EPS_M){
        while (norma > EPS_M){
            h *= 0.5;
            //quiprintf("CASE DECREASE: h = %2.20f  x = %2.5f\n", h, x_cur);
            if (h < H_MIN){
                printf("Error: h is too small.\n");
                exit(-1);
            }
            step(Y_new, Y_new_h, Y_cur, Y_temp, K_1, K_2, F_out, x_cur, h, M, p_1, p_2, alpha_2, beta_2_1);
            //printf("Y_new = %2.20f \n", Y_new[0]);
            NORMA(norma, Y_new, Y_new_h, M);
            //printf("norma = %2.20f \n", norma);
        }
        x_cur += h;
		/*
        if (x_cur > X){
            x_cur = X;
        }
        */
        Y(Y_out, x_cur, M);
        for (int i = 0; i < M; ++i){
            Y_cur[i] = Y_new[i];
        }
        print(f, p, x_cur, h, Y_cur, Y_out, M);
    } else {
        printf("Error: unknown norma value.\n");
        exit(-1);
    }
}

void step(double* &Y_new, double* &Y_new_h, double* &Y_cur, double* Y_temp, double* K_1, double* K_2, double* &F_out, double &x_cur, double &h, int M, double p_1, double p_2, double alpha_2, double beta_2_1){
    //for h 
    F(F_out, x_cur, Y_cur, M);
    for (int i = 0; i < M; ++i){
        K_1[i] = h * F_out[i];
    }
    for (int i = 0; i < M; ++i){
        Y_temp[i] = Y_cur[i] + beta_2_1 * K_1[i];
    }
    F(F_out, x_cur + alpha_2 * h, Y_temp, M);
    for (int i = 0; i < M; ++i){
        K_2[i] = h * F_out[i];
    }
    for (int i = 0; i < M; ++i){
        Y_new[i] = Y_cur[i] + p_1 * K_1[i] + p_2 * K_2[i];
    }
    //printf("h = %2.16f x = %2.16f\n", h, x_cur);
    //printf("K_1[0] = %2.16f\nK_2[0] = %2.16f\nK_1[1] = %2.16f\nK_2[1] = %2.16f\n\n", K_1[0], K_2[0], K_1[1], K_2[1]);
    //printf("Y_new[0] = %2.16f\nY_new[1] = %2.16f\n\n", Y_new[0], Y_new[1]);
    //for h * 0.5 part 1
    h *= 0.5;
    F(F_out, x_cur, Y_cur, M);
    for (int i = 0; i < M; ++i){
        K_1[i] = h * F_out[i];
    }
    for (int i = 0; i < M; ++i){
        Y_temp[i] = Y_cur[i] + beta_2_1 * K_1[i];
    }
    F(F_out, x_cur + alpha_2 * h, Y_temp, M);
    for (int i = 0; i < M; ++i){
        K_2[i] = h * F_out[i];
    }
    for (int i = 0; i < M; ++i){
        Y_new_h[i] = Y_cur[i] + p_1 * K_1[i] + p_2 * K_2[i];
    }
    //printf("h = %2.16f x = %2.16f\n", h, x_cur);
    //printf("h/2 1:\nK_1[0] = %2.16f\nK_2[0] = %2.16f\nK_1[1] = %2.16f\nK_2[1] = %2.16f\n\n", K_1[0], K_2[0], K_1[1], K_2[1]);
    //printf("Y_new_h[0] = %2.16f\nY_new_h[1] = %2.16f\n\n", Y_new_h[0], Y_new_h[1]);
    //for h * 0.5 part 2 
    x_cur += h;
    F(F_out, x_cur, Y_cur, M);
    for (int i = 0; i < M; ++i){
        K_1[i] = h * F_out[i];
    }
    for (int i = 0; i < M; ++i){
        Y_temp[i] = Y_cur[i] + beta_2_1 * K_1[i];
    }
    F(F_out, x_cur + alpha_2 * h, Y_temp, M);
    for (int i = 0; i < M; ++i){
        K_2[i] = h * F_out[i];
    }
    for (int i = 0; i < M; ++i){
        Y_new_h[i] = Y_new_h[i] + p_1 * K_1[i] + p_2 * K_2[i];
    }
    //printf("h = %2.16f x = %2.16f\n", h, x_cur);
    //printf("h/2 2:\nK_1[0] = %2.16f\nK_2[0] = %2.16f\nK_1[1] = %2.16f\nK_2[1] = %2.16f\n\n", K_1[0], K_2[0], K_1[1], K_2[1]);
    //printf("Y_new_h[0] = %2.16f\nY_new_h[1] = %2.16f\n\n", Y_new_h[0], Y_new_h[1]);
    x_cur -= h;
    h *= 2.0;
}

void F(double* &F_out, double x_cur, double* Y_cur, int M){
    /*
    for (int i = 0; i < M; ++i){
        F_out[i] = (-2.0) * A * (x_cur - 5.0) * exp((-1) * A * (x_cur - 5.0) * (x_cur - 5.0))  + Y_cur[i] * 0.0;
    }
    */
    //F_out[M - 1] = Y_cur[0] * Y_cur[0] + x_cur * 0.0;
    F_out[0] = A * Y_cur[1];
	F_out[M - 1] = (-1) * A * Y_cur[0] + x_cur * 0.0;
	//F_out[M - 1] = 2 * x_cur + 1 + 0.0 * Y_cur[0];
}

void Y(double* &Y_out, double x_cur, int M){
    /*
    for (int i = 0; i < M; ++i){
        Y_out[i] = exp((-1) * A * (x_cur - 5.0) * (x_cur - 5.0));
    }
    */
    //Y_out[M - 1] = 1.0 / (1.0 - x_cur);
    Y_out[0] = sin(A * x_cur);
	Y_out[M - 1] = cos(A * x_cur);
	//Y_out[M - 1] = x_cur * x_cur + x_cur + 1;
}

void NORMA(double &norma, double* Y_new, double* Y_new_h, int M){
    norma = 0.0;
    for (int i = 0; i < M; ++i){
        norma += pow(Y_new[i] - Y_new_h[i], 2);
    }
    norma = sqrt(norma);
}

void print(FILE* f, FILE* p, double &x_cur, double &h, double* &Y_cur, double* &Y_out, int M){
	if (x_cur <= X){
		fprintf(f, "x = %1.3f  h = %1.3e  ", x_cur, h);
		fprintf(p, "%2.6f %2.6f ", x_cur, h);
		for (int i = 0; i < M; ++i){
			if (fabs(Y_out[i] - Y_cur[i]) < pow(10, -15)){
				p_deg = 10;
			} else {
				p_deg = (log(fabs(Y_out[i] - Y_cur[i])) / log(h));
			}
			if (p_deg < deg_min){
				deg_min = p_deg;
				deg_min_i = i;
				deg_min_x = x_cur;
			}
			fprintf(f, "y[%d] = %1.3e  y(%1.3f) = %1.3e  del[%d] = %1.3e  p = %1.3f  ", i, Y_cur[i], x_cur, Y_out[i], i, Y_out[i] - Y_cur[i], p_deg);
			fprintf(p, "%.6f ", Y_cur[i]);
		}
		fprintf(f, "\n");
		fprintf(p, "\n");
		/*
		if (x_cur >= X){
			fprintf(f, "%f\n", h);
			fprintf(p, "%d\n", M);
			printf("\n\nthe last delta = %2.10f - %2.10f = %2.10f\n\n", Y_cur[0], Y_out[0], Y_cur[0] - Y_out[0]);
		}
		*/
	}
}
