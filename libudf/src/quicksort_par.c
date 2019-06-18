 #include "quicksort_par.h"
 
 /*---------- quickSort ------------------------------------------
Purpose: sorts the array of nodes a*[] in ascending order     
based on distance from the node to the tip of the exumbrella 
using the "divide and conquer" style quick sort algorithm.   
When the partitions drop below min_list_size, an insertion   
sort algorithm is used to sort the remaining elements.

Input:	Node *a[]	-	array of nodes along a surface			
		int low		-	array index of lower bound
		int high	-	array index of upper bond
----------------------------------------------------------------*/
void quickSort( Node *a[], int low, int high, char meshZone ) 
{
	#if !RP_HOST
	int i, j;
	Node *pivot;
	if ( low == high) return;
	if ( low + MIN_LIST_SIZE > high)
		insertSort(&a[low], high - low + 1, meshZone); /* use insertion sort for small (sub)lists */
	else {
		int mid = (low + high) / 2;
		if (meshZone == 'e') {
			if (N_UDMI(a[mid], 6) < N_UDMI(a[low],6)) swap(&a[mid],  &a[low]);
			if (N_UDMI(a[high],6) < N_UDMI(a[low],6)) swap(&a[high], &a[low]);
			if (N_UDMI(a[high],6) < N_UDMI(a[mid],6)) swap(&a[high], &a[mid]);
		} 
		else if (meshZone == 's') {
			if (N_UDMI(a[mid], 7) < N_UDMI(a[low],7)) swap(&a[mid],  &a[low]);
			if (N_UDMI(a[high],7) < N_UDMI(a[low],7)) swap(&a[high], &a[low]);
			if (N_UDMI(a[high],7) < N_UDMI(a[mid],7)) swap(&a[high], &a[mid]);
		}
		
		pivot = a[mid];
		swap (&a[mid], &a[high-1]); /* move pivot next to high */
		
		i = low + 1;
		j = high - 2;
		for (;;) {
			if (meshZone == 'e' ) {
				/* scan right until a[i] >= pivot */
				while (N_UDMI(a[i],6) < N_UDMI(pivot,6))
					i++;
				
				/* scan left until a[j] <= pivot */
				while (N_UDMI(a[j],6) > N_UDMI(pivot,6))
					j--;
				
				if(i < j)
					swap(&a[i], &a[j]);
				else
					break;
			}
			else if (meshZone == 's') {
				/* scan right until a[i] >= pivot */
				while (N_UDMI(a[i],7) < N_UDMI(pivot,7))
					i++;
				
				/* scan left until a[j] <= pivot */
				while (N_UDMI(a[j],7) > N_UDMI(pivot,7))
					j--;
				
				if(i < j)
					swap(&a[i], &a[j]);
				else
					break;
			}
		}
		swap (&a[i], &a[high-1]); 	/* restore pivot  */
		quickSort(a, low, i-1, meshZone); 	/* sort L sublist */
		quickSort(a, i+1, high, meshZone);	/* sort R sublist */
	}
	#endif
}

/*---------- insertSort ------------------------------------------
Purpose: sorts the array of nodes a*[] in ascending order     
based on distance from the node to the tip of the exumbrella 
using an insertion sort algorithm.

Input:	Node *a[]	-	array of nodes along a surface			
		int n		-	number of elements to be sorted
----------------------------------------------------------------*/
void insertSort( Node *a[], int n, char meshZone) 
{
	#if !RP_HOST
	int i, j; 
	Node *tmp = NULL; 

	for (j = 1; j < n; j++) { 
		tmp = a[j];                            
		i = j - 1;
		if (meshZone == 'e') {
			while (N_UDMI(tmp,6) < N_UDMI(a[i],6)) {
				a[i+1] = a[i]; 
				i--;
				if(i < 0) break;
			}
		}
		else if (meshZone == 's') {
			while (N_UDMI(tmp,7) < N_UDMI(a[i],7)) {
				a[i+1] = a[i]; 
				i--;
				if(i < 0) break;
			}
		}
		a[i+1] = tmp; 
	} 
	#endif
}

/*---------- swap ------------------------------------------------
Purpose: 

Input:	
----------------------------------------------------------------*/
void swap (Node **a, Node **b) 
{
	#if !RP_HOST
	Node *temp = *a;
	*a = *b;
	*b = temp;
	#endif
}

/*---------- fpeek -----------------------------------------------
Purpose: Returns the next character in a file stream without
removing it. Similar to C++ function "peek()".

Input:	FILE *stream	-	File stream
----------------------------------------------------------------*/
int fpeek(FILE *stream) 
{
    int c;

    c = fgetc(stream);
    ungetc(c, stream);

    return c;
}


/*---------- bubbleSort_int ---------------------------------------
Purpose: Sort array arr in ascending order

from https://www.geeksforgeeks.org/bubble-sort

Input:	FILE *stream	-	File stream
----------------------------------------------------------------*/
void bubbleSort_int(int arr[], int n) 
{ 
   int i, j; 
   for (i = 0; i < n-1; i++)       
  
       /* Last i elements are already in place   */
       for (j = 0; j < n-i-1; j++)  
           if (arr[j] > arr[j+1]) 
              swap_int(&arr[j], &arr[j+1]); 
} 

/*---------- bubbleSort_int2D ---------------------------------------
Purpose: Sort array arr in ascending order

from https://www.geeksforgeeks.org/bubble-sort

Input:	FILE *stream	-	File stream
----------------------------------------------------------------*/
void bubbleSort_int2D(int *arr, int rows, int cols) 
{ 
   int i, j, k; 
   for (i = 0; i < rows-1; i++) {	/* Last i elements are already in place   */
       for (j = 0; j < rows-i-1; j++) {
		   if (arr[j*cols+0] > arr[(j+1)*cols+0]) /* ie if (arr[j][1] > arr[j+1][1]) */
		   {
			   for (k = 0; k < cols; k++)
				swap_int(&arr[j*cols + k], &arr[(j+1)*cols + k]); /* ie swap_int(&arr[j][k], &arr[j+1][k]); */
				/* need to change swap func to generalize to >2D array */
		   }
	   }
   }
} 

/*---------- bubbleSort_double ---------------------------------------
Purpose: Sort array arr in ascending order

from https://www.geeksforgeeks.org/bubble-sort

Input:	FILE *stream	-	File stream
----------------------------------------------------------------*/
void bubbleSort_double(double arr[], int n)
{ 
   int i, j; 
   for (i = 0; i < n-1; i++)       
  
       /* Last i elements are already in place   */
       for (j = 0; j < n-i-1; j++)  
           if (arr[j] > arr[j+1]) 
              swap_dub(&arr[j], &arr[j+1]); 
} 

/**---------- bubbleSort_double2D ---------------------------------------
* Purpose: Sort 2D array of doubles arr in ascending/descending order. 
* The array is sorted by the values in the FIRST column.
*
* from https://www.geeksforgeeks.org/bubble-sort
*
* Input:	double *arr -- array of doubles 
*			int rows -- number of rows in arr
*			int cols -- number of cols in arr
*			char direction -- sorting order: 'a' for ascending,
*							  'd' for descending
*
----------------------------------------------------------------**/
void bubbleSort_double2D(double *arr, int rows, int cols, char direction)
{ 
   int i, j, k; 

   for (i = 0; i < rows-1; i++) {	/* Last i elements are already in place   */
       for (j = 0; j < rows-i-1; j++) {
		   if (( (direction == 'a') && (arr[j*cols+0] > arr[(j+1)*cols+0])) ||
		   		((direction == 'd') && (arr[j*cols+0] < arr[(j+1)*cols+0])) )
		   {
			   for (k = 0; k < cols; k++)
				swap_dub(&arr[j*cols + k], &arr[(j+1)*cols + k]); 
		   }
	   }
   }
} 

/*---------- swap_int ---------------------------------------

-------------------------------------------------------------*/
void swap_int(int* a, int* b)
{
    int tmp = *a;
    *a = *b;
    *b = tmp;
}

/*---------- swap_dub ---------------------------------------

-------------------------------------------------------------*/
void swap_dub(double* a, double* b)
{
    double tmp = *a;
    *a = *b;
    *b = tmp;
}

/*---------- swap_any ---------------------------------------

-------------------------------------------------------------*/
void swap_any(void* a, void* b, size_t s)
{
    void* tmp = malloc(s);
    memcpy(tmp, a, s);
    memcpy(a, b, s);
    memcpy(b, tmp, s);
    free(tmp);
}


