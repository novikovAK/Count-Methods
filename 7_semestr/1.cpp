#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int N = 100;
double h=1.0 / double(N-1);


void function ();
double y_function(int m, int k);
double L_function (int m);
double matrix (int i, int j);
double y_function_cont(int m, double x);
void print_y_functions(int q_1, int q_2);
void print_y_function(int q);

int main (void)
{
	function ();	
    print_y_function(2);
	
}

void function ()
{
	int max_accuracy_number = 0, max_scalar_number_1 = 0, max_scalar_number_2;
	double lambda_value = 0.0, norma_1 = 0.0, norma_2 = 0.0, temp = 0.0, accuracy_value = 0.0, max_accuracy_value = 0.0, scalar_value = 0.0, max_scalar_value = 0.0;
	for (int m = 1; m <=N-1; m++)
	{
		for (int s=m+1; s<=N-1; s++) //считаем макс скал€р
		{
			scalar_value = 0.0;
			for (int k = 1; k <= N-1; k++)
			{
				scalar_value += h*y_function(m,k)*y_function(s,k); //m, s - номер функции, k- номер компоненты
			}
			if (fabs(scalar_value)>fabs(max_scalar_value)) 
			{
				max_scalar_value = scalar_value;
				max_scalar_number_1 = m;
				max_scalar_number_2 = s;
			}  
		}
		lambda_value = L_function(m);
		norma_1 = 0.0;
		norma_2 = 0.0;
		for ( int k = 1; k <=N-1; k++) //считаем макс норму
		{
			temp = 0.0;
			for (int s = 1; s <= N-1; s++)
			{
				temp += matrix(k,s)*y_function(m,s);
			}
			norma_1 += pow((temp + lambda_value*y_function(m,k)),2);
			norma_2 += pow((lambda_value*y_function(m,k)),2); 
		}
		accuracy_value = sqrt (norma_1/norma_2);
		if (accuracy_value > max_accuracy_value)
		{
			max_accuracy_value = accuracy_value;
			max_accuracy_number = m; 
		}
	}
	printf ("%d: max accuracy value = %1.5e\n", max_accuracy_number, max_accuracy_value);
	printf ("%d, %d max scalar value = %1.5e\n", max_scalar_number_1, max_scalar_number_2, max_scalar_value);
}

double y_function(int m, int k)
{
	return sin((M_PI*m*(2*k-1))/(2*(N-1)));
}
double L_function (int m)
{
	return (4.0/(h*h))*sin((M_PI*m)/(2*(N-1)))*sin((M_PI*m)/(2*(N-1))); 
} 

double matrix (int i, int j)
{
	if ((i == 1 && j == 1) || (i == N-1 && j == N-1)) 
	{
		return -3.0/(h*h);
	} else if (i == j)
	{
		return -2.0/(h*h);
	} else if (i == j+1 || j == i+1)
	{
		return 1.0/(h*h);
	} else 
	{
		return 0.0;
	}
}

double y_function_cont(int m, double x)
{
	return sin((M_PI * double(m) * (2.0 * x - 1.0)) / (2.0 * (double(N) - 1.0)));
}

void print_y_function(int q)
{
	FILE* f = fopen("out.txt", "w");
	if (f == NULL)
	{
		exit(-1);
	}
	double d = 0.01;
	for (int i = 0; i <= 100 * N; i++)
	{
		fprintf (f, "%f %f \n", d * double(i), y_function_cont(q, d * double(i)));
	}
	fclose(f); 
}


