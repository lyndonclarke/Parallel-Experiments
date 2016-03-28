#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>

/*
 Author: Lyndon Ong Yiu
 * Purpose:  Implement the DAXPY algorithm on two arrays of doubles
 *           x and y, with seperate double alpha, such that:
 *
 *				  for (i = 0; i < n; i++)
 *    			 y[i] += alpha*x[i];
 *
 * Input:    n = size of arrays  x and y
 *           alpha = double to multiply against values from x
 *			 x = array of doubles
 *			 y = array of doubles
 *			 num_threads = number of threads, specified in command line
 * Output:   The values of y after it has been run through the 
 *			 DAXPY protocol

 * Compile:  gcc -g -Wall -o dax daxpy.c
 * Run:      ./dax [num_threads]
 *           
 *
 */
 
 void daxpy(void *num);
 void readArr(double *j, int n);

 double *x, *y, alpha;
 long num_threads;
 int n;

int main(int argc, char *argv[]){
 	long i;
 	pthread_t* threads;

   	num_threads = strtol(argv[1], NULL, 10);
   	threads = malloc(num_threads*sizeof(pthread_t)); 

 	printf("To compute a DAXPY we need two arrays of doubles\nPlease enter the size of the arrays\n");
 	scanf("%d", &n);
 	printf("Please enter alpha\n");
 	scanf("%lf", &alpha);
 	x = malloc(n*sizeof(double));
   	y = malloc(n*sizeof(double));
 	printf("Now enter the values for x\n");
 	readArr(x, n);
 	printf("Now enter the values for y\n");
 	readArr(y, n);

 	for (i = 0; i < num_threads; i++){
 		pthread_create(&threads[i], NULL,(void*) daxpy,(void*) i);
 	}
 	for (i = 0; i < num_threads; i++){
 		pthread_join(threads[i], NULL);
 	}
 	printf("Thanks for your help!\nHere is the newly calculated y array\n");
 	for (i = 0; i < n; i++){
 		printf("%lf ", y[i]);
 	}
 	
 	free(threads);
 	free(x);
 	free(y);
 	return 0;
 }

 void daxpy(void *num){
 	long rank = (long) num;
 	long portion = n / num_threads;
 	long start = portion * rank;
 	long i;
 	for (i = start; i < start + portion; i++){
 		y[i] += alpha * x[i];
 	}


 }

 void readArr(double *j, int n){
 	int i;
 	printf("Enter %d doubles\n", n);
 	for (i = 0; i < n; i = i + 1){
 		scanf("%lf", &j[i]);
 	}

 }

 


