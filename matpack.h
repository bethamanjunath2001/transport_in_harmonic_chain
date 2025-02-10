/*

#############################################################################################################################################################
									Documentation
#############################################################################################################################################################

Matpack is the package written to perform basic operations on matrices such as addition, scalar multiplication, scalar multiplication, transposing, etc.
and several functions are written to perform the same. There are certain guidelines that the user must adhere to for the smooth running of the package.

The package defines a data structure called `matrix' which is supposed to represent the matrix from linear algebra. The structure of matrix and its 
contents are provided in its definition. `matrix' stores the size of the matrix in two int variables 'row' and 'column' and it stores the actual entries 
of matrix inside the object: array of pointers to 'double' arrays (double **).

Someone with sufficient knowledge of the C language can easily see that the object (double **) is just the address of a pointer and not an array. Thus, 
the datastructure should always be followed by the function `matrix_init(int num_row, int num_col, matrix* A)' that initialises it to have appropriate 
size. The user must initialise a declared variable of the type `matrix' with the said function before doing any further operations.

The `matrix_init' functassignment to ‘double **’ from incompatible ion uses calloc() function from <stdlib.h> to allocate the memory required for storing the matrix. Thus, the user is adviced to free the memory that is used to store the matrix after the matrix is no longer required or at the end of the program they write. This will avoid potential 
memory leaks when the program is run for multiple times in a row. A function `free_matrix(matrix *A)' is also written to do the same. The user can proceed 
to use the function to free the matrix defined using matrix_init function, in the following way

	Let us suppose we declare and initialise a matrix variable:
		matrix B;
		matrix_init(numberofrows, numberofcolumns, &B);			//Note that address of B i.e., (&B) is passed instead of B itself
		
		//
		.......some process involving B		
		//
		
		free_matrix(&B);							//Using this line the user can free the memory that was allocated for
										matrix B on the heap. Once this step is completed the user will not be 
										able to access the contents of the matrix B. This should only be used 
										only when the matrix B is no longer required for the program.

The user should keep track of the matrix order that was defined while initialising it. If the user tries to manipulate a matrix out of its bounds, 
it will lead to segmentation faults. 

The functions that take matrices as inputs have an inbuilt option in the form of the argument `int free' to provide a way for the user to automatically 
free the memory occupied by those input matrices. The user can choose to not free and do it later using `free_matrix' function whenever they choose to.

The documentation for the rest of the functions are given below in this header file wherever they are declared. The user can refer to those for guidance.

HAPPY CODING!! CHEERS!!!
*/



//####################################################HEADER FILE WITH FUNCTION DECLARATIONS STARTS HERE###################################################//

#ifndef MATHPACK_H
#define MATHPACK_H



//###########################################################################################//
//		    Defining data types of required for the harmonic chain		     //
//###########################################################################################//

//matrix along with the size//
typedef struct matrix{
	double** array;
	int row, col;
}matrix;


typedef double double4_t __attribute__ ((vector_size (4 * sizeof(double))));

//###########################################################################################//
//		Declaring all the functions that are used for doing matrix operations	     //
//###########################################################################################//

void matrix_init(int num_row, int num_col, matrix *A);
//function to initialise a matrix of given order and dynamically allocate the required memory in heap(only for the array elements).
/*
	args:		num_row				Number of rows
			num_col				Number of columns
			matrix *A			Pointer to a matrix that we want to allocate in the memory and initialise with zeros(calloc function)

	returns:	None
	
	Notes:		
*/



void free_matrix(matrix *A);
/*
	args:		matrix *A			Matrix whose allocated memory is to be freed
	
	returns:	None
	
	Notes:		This function is to be used to free the memory allocated to the matrices so as to avoid memory leaks

*/



void print_matrix(matrix A);								
//function to print a matrix along with order
/*
	args:		matrix A			The matrix to be printed (datatype: matrix)

	returns:	None
	
	Notes:		The function prints out the matrix of data type matrix and its order in a neat manner.
*/



matrix transpose(matrix A, int free_mat);
//function to give the transpose of a given matrix
/*
	args:		matrix A			The matrix that is to be transposed
			int free_mat			The user can decide whether to automatically free the memory of the matrix A that they give as 
							input to the function. 
							free = 1 : free the memory occupied by the matrix A
							free = any other integer : Do not free

	returns:	matrix A_T			Returns the trasposed matrix of A
	
	Notes:		None
*/


matrix prod(matrix A, matrix B, int free_mat);							
//function to multiply two matrices
/*
	args:		matrix A			First matrix (datatype: matrix)
			matrix B			Second matrix (datatype: matrix)
			int free_mat			The user can decide whether to automatically free the memory of the matrix A that they give as 
							input to the function. 
							free = 1 : free the memory occupied by the matrix A
							free = 2 : free the memory occupied by the matrix B
							free = 3 : free the memory occupied by both the matrices
							free = any other integer : Do not free
	
	returns:	matrix product			returns the product of the two matrices in the same order as they are entered in the argument
	
	Notes:		standard matrix multiplication given compatibility in the given order
*/



double dot_prod(matrix A, matrix B, int free_mat);
//function to multiply a row and a column vector to give their dot product
/*
	args:		matrix A			First matrix (datatype: matrix) row matrix
			matrix B			Second matrix (datatype: matrix) column matrix
			int free_mat			The user can decide whether to automatically free the memory of the matrix A that they give as 
							input to the function. 
							free = 1 : free the memory occupied by the matrix A
							free = 2 : free the memory occupied by the matrix B
							free = 3 : free the memory occupied by both the matrices
							free = any other integer : Do not free
	
	returns:	double product		returns the dot product of the given two vectors in the same order as they are entered in the argument
	
	Notes:		standard matrix multiplication given compatibility in the given order	
*/



matrix plus(matrix A, matrix B, int free_mat);							
//function to add two matrices
/*
	args:		matrix A			First matrix (datatype: matrix)
			matrix B			Second matrix (datatype: matrix)
			int free_mat			The user can decide whether to automatically free the memory of the matrix A that they give as 
							input to the function. 
							free = 1 : free the memory occupied by the matrix A
							free = 2 : free the memory occupied by the matrix B
							free = 3 : free the memory occupied by both the matrices
							free = any other integer : Do not free
	
	returns:	matrix sum			returns the sum of the two matrices in the same order as they are entered in the argument
	
	Notes:		Standard matrix addition given they are compatible
*/


matrix minus(matrix A, matrix B, int free_mat);							
//function to subtract two matrices
/*
	args:		matrix A			First matrix (datatype: matrix)
			matrix B			Second matrix (datatype: matrix)
			int free_mat			The user can decide whether to automatically free the memory of the matrix A that they give as 
							input to the function. 
							free = 1 : free the memory occupied by the matrix A
							free = 2 : free the memory occupied by the matrix B
							free = 3 : free the memory occupied by both the matrices
							free = any other integer : Do not free
	
	returns:	matrix diff			returns the difference (not absolute) of the two matrices in the same order
							as they are entered in the argument
	
	Notes:		Standard matrix subtraction given they are compatible
*/



matrix scalX(double a, matrix A, int free_mat);							
//function to multiply a matrix with a real number
/*
	args:		double a			doubleing point number(real number) to scale the matrix
			matrix A			matrix to be scaled (datatype: matrix)
			int free_mat			The user can decide whether to automatically free the memory of the matrix A that they give as 
							input to the function. 
							free = 1 : free the memory occupied by the matrix A
							free = any other integer : Do not free
	
	returns:	matrix prod			Returns the matrix multilied with the given scalar
	
	Notes:		scalar multiplication of a matrix with a given scalar.
*/



matrix make_row(int rowi, matrix A, int free_mat);							
//function to make a row matrix out of a row of a matrix
/*
	args:		int rowi			index of the row which we want to pull apart
			matrix A			matrix from which the row is to be pulled out to make a row vector
			int free_mat			The user can decide whether to automatically free the memory of the matrix A that they give as 
							input to the function. 
							free = 1 : free the memory occupied by the matrix A
							free = any other integer : Do not free
	
	returns:	matrix R			Returns the row with index rowi as a row matrix
	
	Notes:		None
*/


matrix make_col(int coli, matrix A, int free_mat);							
//function to make a column matrix out of a column of a matrix
/*
	args:		int coli			index of the column which we want to pull apart
			matrix A			matrix from which the column is to be pulled out to make a column vector
			int free_mat			The user can decide whether to automatically free the memory of the matrix A that they give as 
							input to the function. 
							free = 1 : free the memory occupied by the matrix A
							free = any other integer : Do not free
	
	returns:	matrix C			Returns the column with index coli as a column matrix
	
	Notes:		None
*/



matrix make_Id(int size);								
//function to make an Identity matrix of given size
/*
	args:		int size			Order of the identity matrix
	
	returns:	matrix Id			Returns the Identity matrix of the order (size X size)
	
	Notes:		None
*/












#endif
