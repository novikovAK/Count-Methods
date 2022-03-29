#include <math.h>
#include <stdlib.h>
#include <stdio.h>

double func(double t, double *y_n, double *k, double h, int i, int M){
    int j;
    double* z  = (double*)malloc(M * sizeof(double));
    for(j = 0; j < M; j++){ 
		z[j] = 0.0; 
	}
    for( j = 0; j < M; j++){
		z[j] = y_n[j] + h * k[j];
    }
    switch (i){
		case 0 :
			return cos(t);
			//return 1.0 + 3.0 * t * t ;
		case 1 :
			return sin(t);
			//return 1.0 + 4.0 * t * t * t;
		default:
			return 0;
    }
}

double Y_acc(double t, int i){
	switch ( i ){
		case 0 :
			return sin(t);
			//return 1 + t + pow(t, 3) ;
		case 1 :
			return -cos(t);
			//return 1 + t + pow(t, 4) ;
		default:
			return 0;
	}
}

double rungekutt(double* &k_n, double* y_n, double t, double h, int M){
    int i;
    double *u = (double*)malloc(M * sizeof(double));
    double* k_n1;
    double* k_n2;
    double* k_n3;
    double* k_n4;
    k_n1 = (double*)malloc(M * sizeof(double));
    k_n2 = (double*)malloc(M * sizeof(double));
    k_n3 = (double*)malloc(M * sizeof(double));
    k_n4 = (double*)malloc(M * sizeof(double));
    for(i=0; i < M; i++){ 
		u[i] = 0.0; 
	}
    for(i=0; i < M; i++){
        k_n1[i] = func(t, y_n, u, 0.0, i, M);
		//printf("k_n1[%d] = %lf\n", i, k_n1[i]);

		k_n2[i] = func(t + h * 0.5, y_n, k_n1, h * 0.5, i, M);
		//printf("k_n2[%d] = %lf\n", i, k_n2[i]);

		k_n3[i] = func(t + h * 0.5, y_n, k_n2, h * 0.5, i, M);
		//printf("k_n3[%d] = %lf\n", i, k_n3[i]);

		k_n4[i] = func(t + h      , y_n, k_n3, h      , i, M);
		//printf("k_n4[%d] = %lf\n", i, k_n4[i]);

		k_n[i] = ( k_n1[i] + 2.0 * k_n2[i] + 2.0 * k_n3[i] + k_n4[i] ) / 6.0;
		//printf("k_n[%d] = %lf\n", i, k_n[i]);
    }
    //printf( "\n");
    //printf( "\n");
    free(k_n1);
    free(k_n2);
    free(k_n3);
    free(k_n4);
    free(u   );
    return 0;
}



double pechat(FILE* out, double t, double* y_n, int M){
    fprintf(out, "%lf\t ", t);
    for (int i = 0; i < M; i++){
        fprintf(out, "%lf\t", y_n[i]);
    }
    fprintf(out, "\n");
    return 0;
}

double step(double* &k_n, double* &y_n, double t, double h, int M){
	int i;
	rungekutt(k_n, y_n, t, h, M );
	for(i = 0; i < M; i++){
		y_n[i] += h * k_n[i];
	}
    return 0;
}

int main(){
	FILE *out = fopen("outrkh.txt", "w");
	int N = 100;
    double a = 0.0;
    double b = M_PI;
    double t;
    //double h = (b-a)/N;
    double h = 0.1;
    t = a;
    double* y_n;
    double* k_n;
    int i;
    int M = 2;
    y_n  = (double*)malloc(M * sizeof(double));
    k_n  = (double*)malloc(M * sizeof(double));
    int k = 0;
    for(i = 0; i < M; i++){
        y_n[i] = Y_acc(t, i);
        k_n[i] = 0.0;
    }
    while (t < b){
        step(k_n, y_n, t, h, M);
		pechat(out, t, y_n, M);
        t += h;
        if (t > b){
            t -= h;
            h = t - b;
            t = b;
		}
    }
    printf( "\n");
    printf( "t_%d = %lf\n", k, t);
    for(i = 0; i < M; i++){
        //printf( "accuracy in t_%d for y%d = %.15lf\n", k, i, NORM(y_n, y_n1, M));
        //printf( "accuracy in t_%d for y%d = %.15lf\n", k, i, y_n[i] - y_n1[i]);
        printf( "accuracy in t_N for y%d = %1.3e\n", i, y_n[i] - Y_acc(t, h));
        printf( "\n");
    }
    printf( "number of iterations k = %d\n", k);
    free(k_n);
    free(y_n);
    fclose(out);
    return 0;
}
