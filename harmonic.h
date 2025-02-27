/*

#############################################################################################################################################################
									Documentation
#############################################################################################################################################################

This package is written to do the computationally intensive tasks dealing with one dimensional Harmonic chain of particles hooked to walls on both of 
its ends. We primarily aim at calculating the quantities related to heat transport in this type of a chain. Lot of mathematical expressions that need to be 
calculated are in terms of matrices of very large sizes. To perform the matrix computations the package matpack is used which can perform very basic matrix 
operations such as multiplication, addition, etc. The package defines matrices as arrays of pointers and has special functions to allocate memory for them and free the memory when their usage is done. The user should take a look at the header file of matpack.h also before jumping to use the package harmonic.h.

Numerous functions in this package take the inputs of type `matrix' from matpack and return the objects of the same type too. The user must properly declare 
and allocate memory for the matrix they want to create themselves for use, using matrix_init function from matpack.

***However for the matrices that user defines solely to read the outputs of functions in this package, the user need not allocate memory using matrix_init 
function (matpack). They can just assign the matrix that they want to read into the matrix variable. 
For the functions in this package where there are intermediate matrix variables needed in the calculations, the memory is properly allocated and freed always.
When the functions in this package output a matrix variable the user should define a matrix variable before and read the output of the function by assigning it to the defined matrix using '='.
They should remember to free the matrix when the job is done.
The user is adviced to free the memory they allocated for any matrix once they are no longer required. The ways to do these operations are described in the header file matpack.h.


//////////////////////////////////////////////////////////////////////////////////////////////////////
Let us look at an example:
////Example starts

matrix A, B;			//initialise two matrix variables A and B.
matrix_init(m, n, &A);		//Allocating a (mXn) sized memory to hold mXn matrix in it and initialising it to all zero values.

//
Do whatever is required with A. Read values out of it and into it. Use it in matrix operations along with other matrices.
//

B = generic_function(list of all the input variables separated by commas)	//When we are reading the output matrix we need not initialise it using matrix_init(). Thus did not initialise B.

//
Do whatever is to be done with matrix B.
//

free_matrix(A);			
free_matrix(B);			//freeing the memory allocated for the matrix A and matrix B

////Example ends
/////////////////////////////////////////////////////////////////////////////////////////////////////
*** 

The functions in this package are declared and documented with all the necessary information in the header file given below. The user can refer to the same 
while trying to use the package.

The units used for each base quantity are as follows. These are chosen so that they will represent the atomic scale of the chains.

Mass M			10^-26 kg						~order of the mass of a carbon atom
Length L		Angstrom unit						~order of the interatomic distances in graphite
Time t			100 femtosecond						~order of the periods of typical molecular vibrations
Temperature T		Kelvin							~SI unit

All the quantities and constants involved are taken to be in these or suitable combinations of these units.

CHEERS!!! FELLOW COMPUTATIONAL PHYSICIST!!!
*/


//####################################################HEADER FILE WITH FUNCTION DECLARATIONS STARTS HERE###################################################//

#ifndef HARMONIC_H
#define HARMONIC_H

#include "matpack.h"
#define PI 3.141592653589793238462643383279502884197					//Numerical value of the constant pi
#define N 234
#define K 1.0										//Spring constant k
#define M 1.0										//Mass of each particle m
#define OMEGA_0 1.0									//Natural frequency of each spring omega_0= k/m
#define TIME 100.0 									//Time elapsed
#define K_B 1.380649e-3									//Boltzmann's constant
#define TEMP 375.0									//Temperature of the system initially
#define VEC_LEN 400

//###########################################################################################//
//##################Defining data types of required for the harmonic chain###################//
//###########################################################################################//
//This datatype is used in our calculations so that multiple operations can be done on multiple scalar variable at once. The operations are initialised and explained below.
typedef struct vector{
	double array[VEC_LEN];
}vector;


vector vector_add(vector A, vector B);

vector vector_subtract(vector A, vector B);

vector scalar_mult(double a, vector A);

vector elem_mult(vector A, vector B);



//###########################################################################################//
//###############Declaring all the functions to be used in the harmonic chain################//
//###########################################################################################//

void logo();										
//function to print the logo
/*
	args: 		None
	
	returns:	None
*/


double omegal(int size, int ind);							
//function to calculate eigenvalues (omega_l)
/*
	args:		int size			:size of the system i.e., number of particles.
			int ind				:index of the eigenvalue which we want to find
			
	returns:	double omegal		:The real variable to store the eigenvalue with the index ind
	
	Notes:		\omega_l = OMEGA_0*sin((ind+1)*PI/(2*(size+1))			:expression used for calculating eigenvalues
*/



matrix make_quad_form_Omega(int size);
//function to make the matrix Omega which is the tridiagonal matrix inside the quadratic form of the Hamiltonian of harmonic chain. Omega_{m n} = 2*delta_{m n} - delta_{m (n-1)} - delta{m (n+1)}.
/*
	args: 		int size 	:size of the system i.e.,number of particles.
	
	returns: 	matrix Omega 	:matrix Omega
*/	
	


matrix make_U(int size);								
//function to make the diagonal transformation matrix U of Omega
/*
	args: 		int size 	:size of the system i.e.,number of particles.
	
	returns: 	matrix U 	:transformation matrix U comprising of all the eigenvectors(modes) of the system.
*/



matrix make_eigenvector(int size, int ind);
//function to make the lth normalised eigenvector of Omega for a given system size.
/*
	args: 		int size 	:size of the system i.e.,number of particles.
			int ind 	:index of the eigenvector
	
	returns: 	matrix U_l 	:lth eigenvector in the form of a column matrix of given size
*/



matrix make_W(int size);								
//function to make the diagonal matrix of eigenvalues W for a given system size.
/*
	args: 		int size 	:size of the system i.e.,number of particles.
	
	returns: 	matrix W 	:diagonal matrix comprising of all the eigenfrequencies of the system.
*/


matrix make_invW(int size);								
//function to make the inverse matrix of W
/*
	args: 		int size 	:size of the system i.e.,number of particles.
	
	returns: 	matrix invW 	:inverse matrix of W i.e., same diagonal matrix but with eigenvalues replaces by their inverses.
*/


matrix make_cosWt(int size, double time);						
//function to make the diagonal matrix of cosines of omega_n*t
/*
	args: 		int size 		:size of the system i.e.,number of particles.
			double time	:time instant at which we want to calculate this time dependent matrix function.
	
	returns: 	matrix cosWt 		:diagonal matrix with each element given by cos(\omega_t) which is a function of time and \omega's are eigenvalues.
*/


matrix make_sinWt(int size, double time);						
//function to make the diagonal matrix of sines of omega_n*t
/*
	args: 		int size 		:size of the system i.e.,number of particles.
			double time	:time instant at which we want to calculate this time dependent matrix function.
	
	returns: 	matrix sinWt 		:diagonal matrix with each element given by sin(\omega_t) which is a function of time and \omega's are eigenvalues.
*/	


matrix gen_rand_x0(int size, bool negatives);								
//function to randomly generate initial displacement vector of the system
/*
	args: 		int size 		:size of the system i.e.,number of particles.
			bool negatives		:This decides if the random numbers generated are in the interval
							(0, 1)  for negatives = false
							(-1, 1) for negatives = true
	
	returns: 	matrix x0 		:generates a column vector with random initial displacements based on the 
						seeding given at the start of the program
*/


matrix gen_rand_xdot0(int size, bool negatives);							
//function to randomly generate initial velocity vector of the system
/*
	args: 		int size 		:size of the system i.e.,number of particles.
			bool negatives		:This decides if the random numbers generated are in the interval
							(0, 1)  for negatives = false
							(-1, 1) for negatives = true
	
	returns: 	matrix xdot0 		:generates a column vector with random initial velocities based on
						the seeding given at the start of the program
*/


double xloft(int size, int part_ind, double time, matrix x0, matrix xdot0, bool tilde);		
//function to calculate dispacement of one particle at a later time given some initial conditions
/*
	args:		int size			:size of the system i.e., number of particles.
			int part_ind			:index of the particle for which we want to find out final displacement
			double time		:final time at which displacements and velocities are to be calculated
			matrix x0			:column vector consisting initial displacements of the particles
			matrix xdot0			:column vector consisting initial velocities of the particles
			bool tilde			:Boolean variable to decide the function is defined on the displacement vectors transformed
							 into the eigenbasis of Omega (tilde) or not. If yes both the inputs and outputs will be in
							 the same basis
			
	returns:	double xl			:two vectors x(t) and v(t) will be written as a single datastructure and can be accessed separately
							 after the double matrix is read into a variable
	
	Notes:		//x_l(t) = (2/N+1)(U_l*cosWt*U*x0 + U_l*sinWt*invW*U*xdot0)			:expression used to calculate x(t)
			We proceed with the calculation exactly as in the function xnxdott but with 1 key distinction that in the final left hand matrix
			product instead of taking the whole of U matrix we only take the row having the same index number as our particle.
*/


double xldotoft(int size, int part_ind, double time, matrix x0, matrix xdot0, bool tilde);		
//function to calculate velocity of one particle at at time given some initial conditions
/*
	args:		int size			:size of the system i.e., number of particles.
			int part_ind			:index of the particle for which we want to find out final displacement
			double time		:final time at which displacements and velocities are to be calculated
			matrix x0			:column vector consisting initial displacements of the particles
			matrix xdot0			:column vector consisting initial velocities of the particles
			bool tilde			:Boolean variable to decide the function is defined on the displacement vectors transformed
							 into the eigenbasis of Omega (tilde) or not. If yes both the inputs and outputs will be in
							 the same basis
			
	returns:	double xldot		:two vectors x(t) and v(t) will be written as a single datastructure and can be accessed separately
							 after the double matrix is read into a variable
	
	Notes:		//dx/dt (t) = -(2/N+1)(U*sinWt*W*U*x0 - U*cosWt*U*xdot0)		:expression used to calculate xdot(t)
			We proceed with the calculation exactly as in the function xnxdott but with 1 key distinction that in the final left hand matrix
			product instead of taking the whole of U matrix we only take the row having the same index number as our particle.
*/



matrix xoft(int size, double time, matrix x0, matrix xdot0, bool tilde);
//function to calculate final displacement vector at a later time given some initial condition vectors
/*
	args:		int size				:size of the system, i.e., number of particles
			double time			:final time variable at which we want to calculate
			matrix x0				:initial displacement vector
			matrix xdot0				:initial velocity vector
			bool tilde				:Boolean variable to decide the function is defined on the displacement vectors transformed
								 into the eigenbasis of Omega (tilde) or not. If yes both the inputs and outputs will be in
								 the same basis
			
	returns:	matrix xoft				:final displacement vector of all particles
	
	Notes:		x(t) = (U*cosWt*U*x0 + U*sinWt*invW*U*xdot0)	:matrix expression used to calculate final displacement vector
*/



matrix xdotoft(int size, double time, matrix x0, matrix xdot0, bool tilde);
//function to calculate final velocity vector at a later time given some initial condition vectors
/*
	args:		int size				:size of the system, i.e., number of particles
			double time			:final time variable at which we want to calculate
			matrix x0				:initial displacement vector
			matrix xdot0				:initial velocity vector
			bool tilde				:Boolean variable to decide the function is defined on the displacement vectors transformed
								 into the eigenbasis of Omega (tilde) or not. If yes both the inputs and outputs will be in
								 the same basis
								 
	returns:	matrix xdotoft				:final velocity vector of all particles
	
	Notes:		dx/dt = -(U*sinWt*W*U*x0 - U*cosWt*U*xdot0)	:matrix expression used to calculate final displacement vector
*/


double e_part(int size, int ind, double time, matrix x0, matrix xdot0);
//function to calculate energy of a single particle in the chain at a time t given the initial condition vectors
/*
	args:		int size				:size of the system
			int ind					:index of the particle for which the energy is to be calculated
			double time				:Time at which to calculate the energy
			matrix x0				:initial displacement vector
			matrix xdot0				:initial velocity vector
			
	returns:	double e_l				:The energy of the particle that will be returned by the function
	
	Notes:		e_l = 1/2 m x_ldot^2 + k/4 ((x_l-x_lminus1)^2 + (x_lplus1 - x_l)^2)		:expression used for calculating energy
*/


double heat_spec_eigen(int size, int bond_ind, int eigen, double time);
//function to calculate the (m=n) part of the IHC(actually, total IHC as the other part vanishes) for the special case when x0 = eigenvector of the system and xdot0 = zero vector
/*
	args:		int size				:size of the system
			int bond_ind				:index of the particle on the left side of the bond through which we are calculating IHC
			int eigen				:Eigenvector that is taken as initial displacement vector
			double time				:Time until which we have to integrate the current
			
	returns:	double h				:Integrated heat current for the specific eigenvector as initial condition and all the
								 initial velocities being zero
								 
	Notes:		None
*/


double J_of_xnxdotoft(int size, int bond_ind, double time, matrix x0, matrix xdot0);
//function to calculate the current through a particular bond at a given instant of time by substituting displacement and velocities of the bonding particles at time t directly as inputs
/*
	args:		int size				:size of the system
			int bond_ind				:index of the particle on the left side of the bond through which we are calculating IHC
			double time			:Time at which we have to find the current
			matrix x0				:initial displacement vector
			matrix xdot0				:initial velocity vector

	returns:	double J				:heat current through a particular bond at a particular instant of time given the initial
								 conditions x0 and xdot0
								 
	Notes:		None
*/


double Joft(int size, int bond_ind, double time, matrix x0, matrix xdot0, bool tilde);
//function to calculate the current through a particular bond at a given instant of time using the expression directly in terms of time
/*
	args:		int size				:size of the system
			int bond_ind				:index of the particle on the left side of the bond through which we are calculating IHC
			double time			:Time at which we have to find the current
			matrix x0				:initial displacement vector
			matrix xdot0				:initial velocity vector

	returns:	double J				:heat current through a particular bond at a particular instant of time using the expression
								 for J in terms of time
								 
	Notes:		None
*/


double Qoft(int size, int bond_ind, double time, matrix x0, matrix xdot0, bool tilde);
//function to calculate the integrated heat current Q(t) through a particular bond over a period of `time' using the expression for time integral of J(t)
/*
	args:		int size				:size of the system
			int bond_ind				:index of the particle on the left side of the bond through which we are calculating IHC
			double time			:Time until which we have to integrate the current
			matrix x0				:initial displacement vector
			matrix xdot0				:initial velocity vector

	returns:	double J				:Integrated heat current through a particular bond over period of time using the expression
								 for Q(t).
								 
	Notes:		Successfully been verified as the correct expression by comparing it with the energy of the particle.
*/


vector Qoftvec(int size, int bond_ind, vector time, matrix x0, matrix xdot0, bool tilde);
//function to calculate the integrated heat current Q(t) through a particular bond over a period of `time' using the expression for time integral of J(t) for a vector comprising several time values at once.
/*
	args:		int size				:size of the system
			int bond_ind				:index of the particle on the left side of the bond through which we are calculating IHC
			vector time				:vector of time values until which we have to integrate the current
			matrix x0				:initial displacement vector
			matrix xdot0				:initial velocity vector

	returns:	double J				:Integrated heat current through a particular bond over period of time using the expression
								 for Q(t).
								 
	Notes:		Successfully been verified as the correct expression by comparing it with the energy of the particle.
*/


double4_t* boxmuller_sampler(int size, int num_samples, double temp);
//function to sample initial position and momentum vectors from Boltzmann distribution --> gaussian distribution for each variable separately. As discussed in theory
/*
	args:		int size				:size of the system or the size of the position and momentum vectors
			int num_samples				:Number of samples of initial conditions
			double temp			:temperature of the system initially
			
	returns: 	matrix x0				:displacement vector of given size
			matrix xdot0				:velocity vector of given size
			
	Notes:		This uses the Box - Muller sampling for sampling a pair of gaussian samples each time and them stitching them together to
			get random displacement and velocity vectors.
*/

double sys_hamiltonian(int size, double time, matrix x0, matrix xdot0, bool tilde);
//function to calculate the total Hamitonian of the system given certain initial conditions
/*
	args:		int size				:size of the system
			double time				:Time elapsed after the given initial conditions
			matrix x0				:Initial displacement vector
			matrix xdot0				:Initial velocity vector
			tilde					:Boolean variable to decide whether the input vectors are transformed in to tilde 
								 versions i.e., the eigen basis of Omega matrix or not
								 
	returns:	double ham				:Total Hamiltonian or energy of the system
	
	Notes: 		The expression in terms of matrices will be used

*/






#endif
