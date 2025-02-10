#include <stdlib.h>
#include <stdio.h>
#include "matpack.h"


/*#######################################################################################*/
/*******function to allocate the required memory and initialise a matrix to zeros*********/
/*#######################################################################################*/
void matrix_init(int num_row, int num_col, matrix *A){

(*A).array = (double **)calloc(num_row, sizeof(double *));

int i;
for (i=0; i<num_row; i++){
	(*A).array[i] = (double *)calloc(num_col, sizeof(double));
}

(*A).row = num_row;
(*A).col = num_col;

}

/*#######################################################################################*/
/*******function to free the memory allocated for the matrix by matrix_init function******/
/*#######################################################################################*/
void free_matrix(matrix *A){

int i;
for (i=0; i<(*A).row; i++){
	free((*A).array[i]);
}

free((*A).array);

}



/*#######################################################################################*/
/**************function to print a matrix of the structured datatype matrix***************/
/*#######################################################################################*/
void print_matrix(matrix A){

int i,j;

if (A.row == 1){
	printf("(");
	i=0;
	for (j=0; j<A.col; j++){
		printf(" %4.8lf ", A.array[i][j]);
	}
	printf(")\n");
} else{

	for (i=0; i<A.row; i++){
		if (i==0){printf("/");} else if (i==(A.row-1)){printf("\\");} else {printf("|");}
		for (j=0; j<A.col; j++){
			printf(" %4.8lf ", A.array[i][j]);
		}
		if (i==0){printf("\\\n");} else if (i==(A.row-1)){printf("/\n");} else {printf("|\n");}
	}
}

printf("ORDER = %d X %d\n", A.row, A.col);

}

/*#######################################################################################*/
/****************function to take a matrix and give its transpose out*********************/
/*#######################################################################################*/
matrix transpose(matrix A, int free){

matrix A_T;
matrix_init(A.col, A.row, &A_T);

int i,j;
for (i=0; i<A.row; i++){
	for (j=0; j<A.col; j++){
		A_T.array[j][i] = A.array[i][j];
	}
}

if (free == 1){
	free_matrix(&A);
}

return A_T;

}




/*#######################################################################################*/
/************function to multiply two matrices of the struct datatype matrix**************/
/*#######################################################################################*/

matrix prod(matrix A, matrix B, int free){

matrix product;
matrix_init(A.row, B.col, &product);
int i, j, k;

matrix C = transpose(B, 0);

if (A.col == B.row){
	for (i=0; i<A.row; i++){
		for (j=0; j<B.col; j++){
			for (k=0; k<A.col; k++){
				product.array[i][j] += (A.array[i][k])*(C.array[j][k]);
			}
		}
	}
} else{
	printf("The multiplication between the given two matrices in the given order is incompatible.\n");
}

free_matrix(&C);

if (free == 1){
	free_matrix(&A);
}else if (free == 2){
	free_matrix(&B);
}else if (free == 3){
	free_matrix(&A);
	free_matrix(&B);
}


return product;

}


/*#################################################################################################*/
/************function to multiply two vectors (row and column) to find the dot product**************/
/*#################################################################################################*/

double dot_prod(matrix A, matrix B, int free){

double dot = 0.0;

if (A.col == B.row){
	for (int i=0; i<A.col; i++){
		dot += (A.array[0][i])*(B.array[i][0]);
	}
} else{
	printf("Dot product is only possible between a row vector and a column vector both of the same size.\n");
}

if (free == 1){
	free_matrix(&A);
}else if (free == 2){
	free_matrix(&B);
}else if (free == 3){
	free_matrix(&A);
	free_matrix(&B);
}

return dot;

}



/*#######################################################################################*/
/***************function to add two matrices of the struct datatype matrix****************/
/*#######################################################################################*/

matrix plus(matrix A, matrix B, int free){

matrix sum;
matrix_init(A.row, A.col, &sum);

int i,j;
if (A.row == B.row && A.col == B.col){
	for (i=0; i<A.row; i++){
		for (j=0; j<A.col; j++){
			sum.array[i][j] = (A.array[i][j]) + (B.array[i][j]);
		}
	}
}else {
	printf("The given two matrices are incompatible for addition.\n");
}

if (free == 1){
	free_matrix(&A);
}else if (free == 2){
	free_matrix(&B);
}else if (free == 3){
	free_matrix(&A);
	free_matrix(&B);
}

return sum;

}

/*#######################################################################################*/
/*************function to add subtract matrices of the struct datatype matrix*************/
/*#######################################################################################*/

matrix minus(matrix A, matrix B, int free){

matrix diff;
matrix_init(A.row, A.col, &diff);

int i,j;
if (A.row == B.row && A.col == B.col){
	for (i=0; i<A.row; i++){
		for (j=0; j<A.col; j++){
			diff.array[i][j] = (A.array[i][j]) - (B.array[i][j]);
		}
	}
}else {
	printf("The given two matrices are incompatible for subtraction.\n");
}

if (free == 1){
	free_matrix(&A);
}else if (free == 2){
	free_matrix(&B);
}else if (free == 3){
	free_matrix(&A);
	free_matrix(&B);
}

return diff;

}

/*#######################################################################################*/
/********************function to multiply a matrix with a given scalar********************/
/*#######################################################################################*/

matrix scalX(double a, matrix A, int free){

matrix prod;
matrix_init(A.row, A.col, &prod);

int i,j;
for (i=0; i< prod.row; i++){
	for (j=0; j< prod.col; j++){
		prod.array[i][j] = a*(A.array[i][j]);
	}
}

if (free == 1){
	free_matrix(&A);
}

return prod;

}

/*#######################################################################################*/
/*******************function to make row matrices out of row of a matrix******************/
/*#######################################################################################*/

matrix make_row(int rowi, matrix A, int free){

matrix R;
matrix_init(1, A.col, &R);

int i;
for(i=0; i<R.col; i++){
	R.array[0][i] = A.array[rowi-1][i];
}

if (free == 1){
	free_matrix(&A);
}

return R;

}

/*#######################################################################################*/
/*******************function to make column matrices out of row of a matrix***************/
/*#######################################################################################*/

matrix make_col(int coli, matrix A, int free){

matrix C;
matrix_init(A.row, 1, &C);

int i;
for(i=0; i<C.row; i++){
	C.array[i][0] = A.array[i][coli-1];
}

if (free == 1){
	free_matrix(&A);
}

return C;

}


/*#######################################################################################*/
/*******************Function to make an Identity matrix of given size*********************/
/*#######################################################################################*/

matrix make_Id(int size){

matrix Id;
matrix_init(size, size, &Id);

int i,j;
for (i=0; i<size; i++){
	for(j=0; j<size; j++){
		if (i==j){
			Id.array[i][j] = 1;
		} else{
			Id.array[i][j] = 0;
		}
	}
}

return Id;

}

































