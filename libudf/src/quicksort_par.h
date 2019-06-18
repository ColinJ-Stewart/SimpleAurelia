#ifndef QUICKSORT_H
#define QUICKSORT_H

#include <udf.h>
/* #include <unsteady.h> */

/* constants */
static const int MIN_LIST_SIZE = 10;

/* prototypes */
void quickSort( Node *a[], int low, int high, char meshZone );
void insertSort( Node *a[], int n, char meshZone);
void swap (Node **a, Node **b);
int fpeek(FILE *stream);
void bubbleSort_int(int arr[], int n);
void bubbleSort_int2D(int *arr, int rows, int cols);
void bubbleSort_double(double arr[], int n);
void bubbleSort_double2D(double *arr, int rows, int cols, char direction);
void swap_int(int* a, int* b);
void swap_dub(double* a, double* b);
void swap_any(void* a, void* b, size_t s);

#endif
