#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int my_pow(int x, int n, int p);

int main(){
	printf("%d\n", (my_pow(294, 8, 353) * 252) % 353);
	return 0;
}

int my_pow(int x, int n, int p){ // возводим х в степень n и берём остаток при делении на p
	int a = 1;
	for (int i = 1; i <= n; ++i){
		a *= x;
		a = a % p;
	}
	return a;
}