#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <omp.h>
#include "harmonic.h"





/*#######################################################################################*/
/****************Function for printing the logo and title of this project*****************/
/*#######################################################################################*/

void logo(){ 
printf("#####################################################################\n");
printf("     ____          ______         __    __  _______  ___   _________ \n");
printf("    /_  /          \\  __ \\       / /   / / / _____/ /   | /___  ___/ \n");
printf("     / /   ____    / /  \\ \\     / /___/ / / /__    / /| |    / /     \n");
printf("    / /   /___/   / /   / /    / ____  / / ___/   / __  |   / /      \n");
printf("  _/ /_          / /___/ /    / /   / / / /____  / /  | |  / /       \n");
printf(" /____/         /_______/    /_/   /_/ /______/ /_/   |_| /_/        \n\n");
printf("#####################################################################\n");
}


//###########################################################################################//
//###############Declaring all the functions to be used in the harmonic chain################//
//###########################################################################################//
//addition
vector vector_add(vector A, vector B){
vector sum;

for (int kv=0; kv<VEC_LEN; kv++){
	sum.array[kv] = A.array[kv] + B.array[kv];
}

return sum;
}

//subtraction
vector vector_subtract(vector A, vector B){
vector diff;

for (int kv=0; kv<VEC_LEN; kv++){
	diff.array[kv] = A.array[kv] - B.array[kv];
}

return diff;
}

//scalar multiplication
vector scalar_mult(double a, vector A){
vector prod;

for (int kv=0; kv<VEC_LEN; kv++){
	prod.array[kv] = a * A.array[kv];
}

return prod;
}

//element by element multiplication
vector elem_mult(vector A, vector B){
vector prod;

for (int kv=0; kv<VEC_LEN; kv++){
	prod.array[kv] = A.array[kv]*B.array[kv];
}

return prod;
}


//#######################################################################################//
//****************Function to calculate eigen values \omega_l*****************************/
//#######################################################################################//
double omegal(int size, int ind) {

double omegal;
 
omegal = ind;
omegal *= PI;
omegal /= 2.0;
omegal /= (size+1.0);
omegal = sin(omegal);
omegal *= OMEGA_0;
omegal *= 2.0;

return omegal;

}



/*############################################################################################*/
/********Function to make the matrix in the middle of quadratic form in the Hamilonian*********/
/*############################################################################################*/
matrix make_quad_form_Omega(int size) {

matrix Omega;
matrix_init(size, size, &Omega);


int i,j;
for (i=0; i<size; i++){
	for (j=0; j<size; j++){
		if (i==j){
			Omega.array[i][j] = 2.0;
		} else if (i == j+1 || i == j-1){
			Omega.array[i][j] = -1.0;
		} else {
			Omega.array[i][j] = 0.0;
		}
	}

}

return Omega;

}



/*#######################################################################################*/
/*********************Function to make the transformation matrix U************************/
/*#######################################################################################*/

matrix make_U(int size){

matrix U;
matrix_init(size, size, &U);

int i,j;
for (i=0; i<size; i++){
	for(j=0; j<size; j++){
		U.array[i][j] = (i+1.0);
		U.array[i][j] *= (j+1.0);
		U.array[i][j] *= PI;
		U.array[i][j] /= (size+1.0);
		U.array[i][j] = sin(U.array[i][j]);
		U.array[i][j] /= sqrt((size+1.0)/2.0);
	}
}

return U;

}



/*#######################################################################################*/
/*********************Function to make the transformation matrix U************************/
/*#######################################################################################*/
matrix make_eigenvector(int size, int ind){

matrix U_l;
matrix_init(size, 1, &U_l);

int i;
for (i=0; i<size; i++){
	U_l.array[i][0] = ind;
	U_l.array[i][0] *= (i+1.0);
	U_l.array[i][0] *= PI; 
	U_l.array[i][0] /= (size+1.0);
	U_l.array[i][0] = sin(U_l.array[i][0]);
	U_l.array[i][0] /= sqrt((size+1.0)/2.0);
}

return U_l;

}




/*#######################################################################################*/
/******************Function to make the eigenvalue diagonal matrix W**********************/
/*#######################################################################################*/

matrix make_W(int size){

matrix W;
matrix_init(size, size, &W);

int i,j;
for (i=0; i<size; i++){
	for(j=0; j<size; j++){
		if (i==j){
			W.array[i][j] = omegal(size, i+1.0);
		} else{
			W.array[i][j] = 0.0;
		}
	}
}

return W;

}



/*#######################################################################################*/
/******************Function to make the inverse matrix of the matrix W********************/
/*#######################################################################################*/

matrix make_invW(int size){

matrix invW;
matrix_init(size, size, &invW);

int i,j;
for (i=0; i<size; i++){
	for(j=0; j<size; j++){
		if (i==j){
			invW.array[i][j] = 1/omegal(size, i+1.0);
		} else{
			invW.array[i][j] = 0.0;
		}
	}
}

return invW;

}



/*#######################################################################################*/
/************Function to make the diagonal matrix of cosines of omega_nt's****************/
/*#######################################################################################*/

matrix make_cosWt(int size, double time){

matrix cosWt;
matrix_init(size, size, &cosWt);

int i,j;
for (i=0; i<size; i++){
	for(j=0; j<size; j++){
		if (i==j){
			cosWt.array[i][j] = omegal(size, i+1.0);
			cosWt.array[i][j] *= time;
			cosWt.array[i][j] = cos(cosWt.array[i][j]);
		} else{
			cosWt.array[i][j] =  0.0;
		}
	}
}

return cosWt;

}




/*#######################################################################################*/
/************Function to make the diagonal matrix of the sines of omega_nt's**************/
/*#######################################################################################*/

matrix make_sinWt(int size, double time){


matrix sinWt;
matrix_init(size, size, &sinWt);

int i,j;
for (i=0; i<size; i++){
	for(j=0; j<size; j++){
		if (i==j){
			sinWt.array[i][j] = omegal(size, i+1.0);
			sinWt.array[i][j] *= time;
			sinWt.array[i][j] = sin(sinWt.array[i][j]);
		} else{
			sinWt.array[i][j] =  0.0;
		}
	}
}

return sinWt;

}




/*#######################################################################################*/
/************Function to generate random initial displacements of the particles***********/
/*#######################################################################################*/

matrix gen_rand_x0(int size, bool negatives){

matrix x0;
matrix_init(size, 1, &x0);

int i;
for (i=0; i<size; i++){
	if(negatives == true){
		x0.array[i][0] = random() % 1000000000;
		x0.array[i][0] /= (500000000.0);
		x0.array[i][0] -= 1.0;
	}else {
		x0.array[i][0] = random() % 1000000000;
		x0.array[i][0] /= 1000000000.0;
	}
}

return x0;

}


/*#######################################################################################*/
/************Function to generate random initial velocities of the particles**************/
/*#######################################################################################*/

matrix gen_rand_xdot0(int size, bool negatives){

matrix xdot0;
matrix_init(size, 1, &xdot0);

int i;
for (i=0; i<size; i++){
	if(negatives == true){
		xdot0.array[i][0] = random() % 1000000000;
		xdot0.array[i][0] /= (500000000.0);
		xdot0.array[i][0] -= 1.0;
	}else {
		xdot0.array[i][0] = random() % 1000000000;
		xdot0.array[i][0] /= 1000000000.0;
	}
}

return xdot0;

}




/*#######################################################################################*/
/*******Function to compute the final displacement of lth particle given everything*******/
/*#######################################################################################*/

double xloft(int size, int part_ind, double time, matrix x0, matrix xdot0, bool tilde) {

matrix U, U_l, invW, cosWt, sinWt;

if (tilde == true){
	U = make_Id(size);
} else{
	U = make_U(size);
}

U_l = make_row(part_ind, U, 0);
invW = make_invW(size);
cosWt = make_cosWt(size, time);
sinWt = make_sinWt(size, time);

matrix x1, x;

x1 = prod(U, x0, 0);									
x1 = prod(cosWt, x1, 3);
x1 = prod(U_l, x1, 2);				
x = prod(U, xdot0, 1);
x = prod(invW, x, 3);
x = prod(sinWt, x, 3);
x = prod(U_l, x, 3);
x = plus(x1, x, 3);

double xl = x.array[0][0];

free_matrix(&x);

return xl;

}




/*#######################################################################################*/
/*********Function to compute the final velocity of lth particle given everything*********/
/*#######################################################################################*/

double xldotoft(int size, int part_ind, double time, matrix x0, matrix xdot0, bool tilde) {

matrix U, U_l, W, cosWt, sinWt;

if (tilde == true){
	U = make_Id(size);
} else{
	U = make_U(size);
}

U_l = make_row(part_ind, U, 0);
W = make_W(size);
cosWt = make_cosWt(size, time);
sinWt = make_sinWt(size, time);

matrix xdot1, xdot;

xdot1 = prod(U, x0, 0);
xdot1 = prod(W, xdot1, 3);
xdot1 = prod(sinWt, xdot1, 3);
xdot1 = prod(U_l, xdot1, 2);
xdot = prod(U, xdot0, 1);
xdot = prod(cosWt, xdot, 3);
xdot = prod(U_l, xdot, 3);
xdot = minus(xdot, xdot1, 3);

double xldot = xdot.array[0][0];

free_matrix(&xdot);

return xldot;

}




/*#######################################################################################*/
/*******Function to compute the final displacement vector of all th particles*************/
/*#######################################################################################*/

matrix xoft(int size, double time, matrix x0, matrix xdot0, bool tilde){

matrix U, invW, cosWt, sinWt;

if (tilde == true){
	U = make_Id(size);
} else{
	U = make_U(size);
}

invW = make_invW(size);
cosWt = make_cosWt(size, time);
sinWt = make_sinWt(size, time);

matrix x1, x;

x1 = prod(U, x0, 0);									
x1 = prod(cosWt, x1, 3);
x1 = prod(U, x1, 2);				
x = prod(U,xdot0, 0);
x = prod(invW, x, 3);
x = prod(sinWt, x, 3);
x = prod(U, x, 3);
x = plus(x1, x, 3);

return x;

}




/*#######################################################################################*/
/***********Function to compute the final velocity vector of all th particles*************/
/*#######################################################################################*/

matrix xdotoft(int size, double time, matrix x0, matrix xdot0, bool tilde){

matrix U, W, cosWt, sinWt;

if (tilde == true){
	U = make_Id(size);
} else{
	U = make_U(size);
}

W = make_W(size);
cosWt = make_cosWt(size, time);
sinWt = make_sinWt(size, time);

matrix xdot1, xdot;

xdot1 = prod(U, x0, 0);
xdot1 = prod(W, xdot1, 3);
xdot1 = prod(sinWt, xdot1, 3);
xdot1 = prod(U, xdot1, 2);
xdot = prod(U, xdot0, 0);
xdot = prod(cosWt, xdot, 3);
xdot = prod(U, xdot, 3);
xdot = minus(xdot, xdot1, 3);

return xdot;

}




/*#######################################################################################*/
/****Function to calculate energy of 1 particle at a given time and initial conditions****/
/*#######################################################################################*/

double e_part(int size, int ind, double time, matrix x0, matrix xdot0){

double term = xloft(size, ind+1, time, x0, xdot0, false);
term -= xloft(size, ind, time, x0, xdot0, false);
term *= term;
double sum = term;
term = xloft(size, ind, time, x0, xdot0, false);
term -= xloft(size, ind-1, time, x0, xdot0, false);
term *= term;
sum += term;
sum *= K;
sum *= 0.25;
term = xldotoft(size, ind, time, x0, xdot0, false);
term *= term;
term *= M;
term *= 0.5;
sum += term;

return sum;

}




/*###################################################################################################*/
/****Function to calculate the current through a particular bond as a function of final conditions****/
/*###################################################################################################*/

double J_of_xnxdotoft(int size, int bond_ind, double time, matrix x0, matrix xdot0){

double term, sum;
sum = xloft(size, bond_ind+1, time, x0, xdot0, false);
sum -= xloft(size, bond_ind, time, x0, xdot0, false);
term = sum;
sum = xldotoft(size, bond_ind+1, time, x0, xdot0, false);
sum += xldotoft(size, bond_ind, time, x0, xdot0, false);
term *= sum;
term *= -0.5;
term *= K;

return term;
}



/*###################################################################################################*/
/***Function to calculate the current through a particular bond as a function of initial conditions***/
/*###################################################################################################*/

double Joft(int size, int bond_ind, double time, matrix x0, matrix xdot0, bool tilde){

matrix U, U_m, U_n;

if (tilde == true){
	U = make_Id(size);
}else {
	U = make_U(size);
}

double term, sum, pro;
double J = 0.0;

#pragma omp parallel for collapse(2) private(term, sum, pro, U_m, U_n) reduction(+:J)
for (int i=0; i<size; i++){
	for (int j=0; j<size; j++){
		U_m = make_row(i+1, U, 0);
		U_n = make_row(j+1, U, 0);
		
		term = dot_prod(U_m, xdot0, 0);
		term *= dot_prod(U_n, xdot0, 0);
		term /= omegal(size, j+1.0);
		sum = term;
		term = dot_prod(U_m, x0, 0);
		term *= dot_prod(U_n, x0, 0);
		term *= omegal(size, i+1.0);
		sum -= term;
		term = omegal(size, i+1.0);
		term += omegal(size, j+1.0);
		term *= time;
		term = sin(term);
		sum *= term;
		pro = sum;
		
		term = dot_prod(U_m, xdot0, 0);
		term *= dot_prod(U_n, xdot0, 0);
		term /= omegal(size, j+1.0);
		sum  = term;
		term = dot_prod(U_m, x0, 0);
		term *= dot_prod(U_n, x0, 0);
		term *= omegal(size, i+1.0);
		sum += term;
		term = omegal(size, i+1.0);
		term -= omegal(size, j+1.0);
		term *= time;
		term = sin(term);
		sum *= term;
		pro -= sum;
		
		term = dot_prod(U_m, xdot0, 0);
		term *= dot_prod(U_n, x0, 0);
		sum = term;
		term = dot_prod(U_m, x0, 0);
		term *= dot_prod(U_n, xdot0, 0);
		term *= omegal(size, i+1.0);
		term /= omegal(size, j+1.0);
		sum += term;
		term = omegal(size, i+1.0);
		term += omegal(size, j+1.0);
		term *= time;
		term = cos(term);
		sum *= term;
		pro += sum;
		
		term = dot_prod(U_m, xdot0, 0);
		term *= dot_prod(U_n, x0, 0);
		sum = term;
		term = dot_prod(U_m, x0, 1);
		term *= dot_prod(U_n, xdot0, 1);
		term *= omegal(size, i+1.0);
		term /= omegal(size, j+1.0);
		sum -= term;
		term = omegal(size, i+1.0);
		term -= omegal(size, j+1.0);
		term *= time;
		term = cos(term);
		sum *= term;
		pro += sum;
		
		term = 2.0*bond_ind;
		term += 1;
		term *= (i+1.0);
		term *= PI;
		term /= 2.0;
		term /= (size+1.0);
		term = sin(term);
		pro *= term;
		
		term = i+1.0;
		term *= PI;
		term /= 2.0;
		term /= (size+1.0);
		term = cos(term);
		pro *= term;
		
		term = 2.0*bond_ind;
		term += 1;
		term *= (j+1.0);
		term *= PI;
		term /= 2.0;
		term /= (size+1.0);
		term = cos(term);
		pro *= term;
		
		term = j+1.0;
		term *= PI;
		term /= 2.0;
		term /= (size+1.0);
		term = sin(term);
		pro *= term;
		
		J += pro;
	}
}

J *= K;
J *= -2.0;
J /= (size+1.0);

free_matrix(&U);

return J;
}




/*######################################################################################################*/
/**Function to calculate the IHC through a particular bond as a function of initial conditions and time**/
/*######################################################################################################*/


double Qoft(int size, int bond_ind, double time, matrix x0, matrix xdot0, bool tilde){

matrix U;
if (tilde == true){
	U = make_Id(size);
}else {
	U = make_U(size);
}
double Q = 0.0;

#pragma omp parallel for collapse(2) reduction(+:Q) num_threads(64)
for(int i=0; i<size; i++){
	for(int j=0; j<size; j++){
		matrix U_m = make_row(i+1, U, 0);
		matrix U_n = make_row(j+1, U, 0);
		
		double term, sum, pro;
		pro = 0.0;
		#pragma omp parallel sections private(term, sum) reduction(+:pro)
		{	
			#pragma omp section
			{
				term = dot_prod(U_m, xdot0, 0);
				term *= dot_prod(U_n, xdot0, 0);
				term /= omegal(size, j+1.0);
				sum = term;
				term = dot_prod(U_m, x0, 0);
				term *= dot_prod(U_n, x0, 0);
				term *= omegal(size, i+1.0);
				sum -= term;
				term = omegal(size, i+1.0);
				term += omegal(size, j+1.0);
				sum /= term;
				term *= time;
				term = cos(term);
				term -= 1.0; 
				sum *= term;
				pro += sum;
			}
			
			#pragma omp section
			{
				if(i!=j){
					term = dot_prod(U_m, xdot0, 0);
					term *= dot_prod(U_n, xdot0, 0);
					term /= omegal(size, j+1.0);
					sum  = term;
					term = dot_prod(U_m, x0, 0);
					term *= dot_prod(U_n, x0, 0);
					term *= omegal(size, i+1.0);
					sum += term;
					term = omegal(size, i+1.0);
					term -= omegal(size, j+1.0);
					sum /= term;
					term *= time;
					term = cos(term);
					term -= 1.0;
					sum *= term;
					sum *= -1.0;
					pro += sum;
				}else{
					pro += 0.0;
				}
			}
			
			#pragma omp section
			{
				term = dot_prod(U_m, xdot0, 0);
				term *= dot_prod(U_n, x0, 0);
				sum = term;
				term = dot_prod(U_m, x0, 0);
				term *= dot_prod(U_n, xdot0, 0);
				term *= omegal(size, i+1.0);
				term /= omegal(size, j+1.0);
				sum += term;
				term = omegal(size, i+1.0);
				term += omegal(size, j+1.0);
				sum /= term;
				term *= time;
				term = sin(term);
				sum *= term;
				sum *= -1.0;
				pro += sum;
			}
			
			#pragma omp section
			{
				if (i!=j){
					term = dot_prod(U_m, xdot0, 0);
					term *= dot_prod(U_n, x0, 0);
					sum = term;
					term = dot_prod(U_m, x0, 0);
					term *= dot_prod(U_n, xdot0, 0);
					term *= omegal(size, i+1.0);
					term /= omegal(size, j+1.0);
					sum -= term;
					term = omegal(size, i+1.0);
					term -= omegal(size, j+1.0);
					sum /= term;
					term *= time;
					term = sin(term);
					sum *= term;
					sum *= -1.0;
					pro += sum;
				} else{
					pro += 0.0;
				}
			}
		}
		free_matrix(&U_m);
		free_matrix(&U_n);
		
		#pragma omp parallel sections private(term) reduction(*:pro)
		{
			#pragma omp section
			{
				term = 2.0*bond_ind;
				term += 1;
				term *= (i+1.0);
				term *= PI;
				term /= 2.0;
				term /= (size+1.0);
				term = sin(term);
				pro *= term;
			}
			
			#pragma omp section
			{
				term = i+1.0;
				term *= PI;
				term /= 2.0;
				term /= (size+1.0);
				term = cos(term);
				pro *= term;
			}
			
			#pragma omp section
			{
				term = 2.0*bond_ind;
				term += 1;
				term *= (j+1.0);
				term *= PI;
				term /= 2.0;
				term /= (size+1.0);
				term = cos(term);
				pro *= term;
			}
			
			#pragma omp section
			{
				term = j+1.0;
				term *= PI;
				term /= 2.0;
				term /= (size+1.0);
				term = sin(term);
				pro *= term;
			}
		}
		
		Q += pro;
	}
}

Q *= K;
Q *= 2.0;
Q /= (size + 1.0);

free_matrix(&U);

return Q;

}



/*######################################################################################################################*/
/**Function to calculate the IHC through a particular bond as a function of initial conditions and a set of time values**/
/*######################################################################################################################*/


vector Qoftvec(int size, int bond_ind, vector time, matrix x0, matrix xdot0, bool tilde){

matrix U;
if (tilde == true){
	U = make_Id(size);
}else {
	U = make_U(size);
}

vector ZERO_VECTOR;

for (int kv=0; kv<VEC_LEN; kv++){
	ZERO_VECTOR.array[kv] = 0.0;
}

vector omp_orig = ZERO_VECTOR;

vector Q = ZERO_VECTOR;

#pragma omp declare reduction(add_vectors: vector: omp_out=vector_add(omp_out, omp_in)) initializer(omp_priv=omp_orig)

#pragma omp parallel for collapse(2) reduction(add_vectors:Q) num_threads(4)
for(int i=0; i<size; i++){
	for(int j=0; j<size; j++){

		double pro = 1.0;
		
		#pragma omp parallel sections reduction(*:pro)
		{
			#pragma omp section
			{
				double term = 2.0*bond_ind;
				term += 1;
				term *= (i+1.0);
				term *= PI;
				term /= 2.0;
				term /= (size+1.0);
				term = sin(term);
				pro *= term;
			}
			
			#pragma omp section
			{
				double term = i+1.0;
				term *= PI;
				term /= 2.0;
				term /= (size+1.0);
				term = cos(term);
				pro *= term;
			}
			
			#pragma omp section
			{
				double term = 2.0*bond_ind;
				term += 1;
				term *= (j+1.0);
				term *= PI;
				term /= 2.0;
				term /= (size+1.0);
				term = cos(term);
				pro *= term;
			}
			
			#pragma omp section
			{
				double term = j+1.0;
				term *= PI;
				term /= 2.0;
				term /= (size+1.0);
				term = sin(term);
				pro *= term;
			}
		}
		
		matrix U_m = make_row(i+1, U, 0);
		matrix U_n = make_row(j+1, U, 0);
		
		vector pro_vector = ZERO_VECTOR;
		
			
		#pragma omp parallel sections reduction(add_vectors:pro_vector) 
		{	
			
			#pragma omp section
			{
				double term = dot_prod(U_m, xdot0, 0);
				term *= dot_prod(U_n, xdot0, 0);
				term /= omegal(size, j+1.0);
				double sum = term;
				term = dot_prod(U_m, x0, 0);
				term *= dot_prod(U_n, x0, 0);
				term *= omegal(size, i+1.0);
				sum -= term;
				term = omegal(size, i+1.0);
				term += omegal(size, j+1.0);
				sum /= term;
				vector term_vector;
				for (int kv=0; kv<VEC_LEN; kv++){ 
					term_vector.array[kv] = term * (time.array[kv]);
					term_vector.array[kv] = cos(term_vector.array[kv]);
					term_vector.array[kv] -= 1.0;
					term_vector.array[kv] *= sum;
				}
				pro_vector = vector_add(pro_vector, term_vector);
			}
			
			#pragma omp section
			{
				if(i!=j){
					double term = dot_prod(U_m, xdot0, 0);
					term *= dot_prod(U_n, xdot0, 0);
					term /= omegal(size, j+1.0);
					double sum  = term;
					term = dot_prod(U_m, x0, 0);
					term *= dot_prod(U_n, x0, 0);
					term *= omegal(size, i+1.0);
					sum += term;
					term = omegal(size, i+1.0);
					term -= omegal(size, j+1.0);
					sum /= term;
					vector term_vector;
					for (int kv=0; kv<VEC_LEN; kv++){
						term_vector.array[kv] = term * (time.array[kv]);
						term_vector.array[kv] = cos(term_vector.array[kv]);
						term_vector.array[kv] -= 1.0;
						term_vector.array[kv] *= sum;
						term_vector.array[kv] *= -1.0;
					}
	 				pro_vector = vector_add(pro_vector, term_vector);
				}else{
					pro_vector = vector_add(pro_vector, ZERO_VECTOR);
				}
			}
			
			#pragma omp section
			{
				double term = dot_prod(U_m, xdot0, 0);
				term *= dot_prod(U_n, x0, 0);
				double sum = term;
				term = dot_prod(U_m, x0, 0);
				term *= dot_prod(U_n, xdot0, 0);
				term *= omegal(size, i+1.0);
				term /= omegal(size, j+1.0);
				sum += term;
				term = omegal(size, i+1.0);
				term += omegal(size, j+1.0);
				sum /= term;
				vector term_vector;
				for (int kv=0; kv<VEC_LEN; kv++) {
					term_vector.array[kv] = term * (time.array[kv]);
					term_vector.array[kv] = sin(term_vector.array[kv]);
					term_vector.array[kv] *= sum;
					term_vector.array[kv] *= -1.0;
				}
				pro_vector = vector_add(pro_vector, term_vector);
			}
			
			#pragma omp section
			{
				if (i!=j){
					double term = dot_prod(U_m, xdot0, 0);
					term *= dot_prod(U_n, x0, 0);
					double sum = term;
					term = dot_prod(U_m, x0, 0);
					term *= dot_prod(U_n, xdot0, 0);
					term *= omegal(size, i+1.0);
					term /= omegal(size, j+1.0);
					sum -= term;
					term = omegal(size, i+1.0);
					term -= omegal(size, j+1.0);
					sum /= term;
					vector term_vector;
					for (int kv=0; kv<VEC_LEN; kv++){
						term_vector.array[kv] = term * (time.array[kv]);
						term_vector.array[kv] = sin(term_vector.array[kv]);
						term_vector.array[kv] *= sum;
						term_vector.array[kv] *= -1.0;
					}
					pro_vector = vector_add(pro_vector, term_vector);
				} else{
					pro_vector = vector_add(pro_vector, ZERO_VECTOR);
				}
			}
		}
		
		free_matrix(&U_m);
		free_matrix(&U_n);
		
		pro_vector = scalar_mult(pro, pro_vector);
		
		Q = vector_add(Q , pro_vector);
	}
}

for (int kv=0; kv<VEC_LEN; kv++){
	Q.array[kv] *= K;
	Q.array[kv] *= 2.0;
	Q.array[kv] /= (size + 1.0);
}

free_matrix(&U);

return Q;

}

/*##################################################################################################*/
/****Function to sample displacement and velocity vectors of given size using Box-Muller sampling****/
/*##################################################################################################*/
double4_t* boxmuller_sampler(int size, int num_samples, double temp){


int nb = 4;
int na = (num_samples + nb -1)/ nb;

double4_t* gauss_samples = (double4_t *)calloc(na*size*2, sizeof(double4_t));
int inf_counter = 0;



for (int j=0; j<na; j++){
	for (int i=0; i<size; i++){
		double4_t u1, u2, xtilde, xtildedot;
		for (int kb=0; kb<nb; kb++){
			u1[kb] = ((double)rand())/ RAND_MAX;
			u2[kb] = ((double)rand())/ RAND_MAX;
		}
		
		for (int kb=0; kb<nb; kb++){
			u1[kb] = log(u1[kb]);
		}
		u1 *= (-2.0);
		u1 *= K_B;
		u1 *= temp;
		u1 /= M;
		for (int kb=0; kb<nb; kb++){
			u1[kb] = sqrt(u1[kb]);
		}
		
		u2 *= 2.0;
		u2 *= PI;

		for (int kb=0; kb<nb; kb++){
			xtilde[kb] = sin(u2[kb]);
		}
		xtilde *= u1;
		xtilde /= omegal(size, i+1.0);
		
		for (int kb=0; kb<nb; kb++){
			xtildedot[kb] = cos(u2[kb]);
		}
		xtildedot *= u1;
		
		for (int kb=0; kb<nb; kb++){
			inf_counter += (xtilde[kb] == INFINITY) ? 1 : 0 ;
		}
		for (int kb=0; kb<nb; kb++){
			inf_counter += (xtildedot[kb] == INFINITY) ? 1 : 0 ;
		}
		/*
		*/
		
		gauss_samples[size*(2*j + 0) + i] = xtilde;
		gauss_samples[size*(2*j + 1) + i] = xtildedot;
	}
	printf("%d initial conditions sampled\n", (j+1)*nb);
}

printf("\nnumber of infinities = %d\n", inf_counter);

return gauss_samples;


}



/*######################################################################################################*/
/****Function to calculate the total Hamiltonian of the system given the system configuration vectors****/
/*######################################################################################################*/
double sys_hamiltonian(int size, double time, matrix x0, matrix xdot0, bool tilde){

double ham;
matrix x, xdot, term, Omega;

x = xoft(size, time, x0, xdot0, tilde);
xdot = xdotoft(size, time, x0, xdot0, tilde);

term = transpose(xdot, 0);
term = prod(term, xdot, 3);
term = scalX(M/2.0, term, 1);
ham = term.array[0][0];
free_matrix(&term);

if (tilde == true){
	Omega = make_W(size);
	Omega = prod(Omega, Omega, 1);
	Omega = scalX(M/2.0, Omega, 1);
}else {
	Omega = make_quad_form_Omega(size);
	Omega = scalX(K/2.0, Omega, 1);
}
Omega = prod(Omega, x, 1);
term = transpose(x, 1);
term = prod(term, Omega, 3);
ham += term.array[0][0];
free_matrix(&term);

return ham;
}





























