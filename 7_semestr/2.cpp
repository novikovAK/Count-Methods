#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int N = 100;

double y_function_cont(int m, double x);
void print_y_function(int q);

int main(){
	print_y_function(1);
	return 0;
}

double y_function_cont(int m, double x)
{
	return sin((M_PI * double(m) * (2.0 * x - 1.0)) / (2.0 * (double(N) - 1.0)));
	//return sin((M_PI * double(m) * x) / (2.0 * double(N)));
}

void print_y_function(int q)
{
	FILE* f = fopen("out.txt", "w");
	if (f == NULL)
	{
		exit(-1);
	}
	double d = 0.01;
	for (int i = 0; i <= 10000; i++)
	{
		fprintf (f, "%f %f \n", d * double(i), y_function_cont(q, d * double(i)));
	}
	fclose(f); 
}
