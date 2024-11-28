#include <iostream>
#include <map>
#include <queue>
#include <set>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <utility>
#include "eigen-3.3.9/Eigen/Sparse"
#include "eigen-3.3.9/Eigen/dense"
#include "eigen-3.3.9/Eigen/Eigenvalues"
#include "eigen-3.3.9/Eigen/Core"

using namespace std;
using namespace Eigen;

const double pi = M_PI;

double pot_1D(double x){
	double r = sqrt(x * x );

	return 0.5 * r * r; // Potential for 1D Harmonic Oscillator
}

int main(int argc, char* argv[])
    {
	int state_number = 20; //Number of excited states 

	double hbar = 1; //Plank constant = 1
	double m = 1; // Mass = 1

	int n = 500; // Grid base
	int N = 2 * n + 1; // Total number of grid points

	double dx = 0.05; // variable parameter Delta x as mentioned in theory
	double dk = (2 * pi)/ (N*dx); // definition of Wavenumber in terms of Grids
    // To allocate the memory for Kinetic Energy 
	double *T_i = (double*)malloc(sizeof(double) * (N + 2));
	for(int di = 0; di <= N; di++)
	{
		T_i [di] = 0;
		for(int j = 0; j <= n; j++)
		{
			double T_k = (hbar * j * dk) * (hbar * j * dk) / (2.0 * m);
			T_i[di] += 2.0 * cos(2.0 * pi * j * di / N) / N * T_k;
		}
	}
    // To construct the Hamiltonian Matrix of 1D Harmonic Oscillator of N*N dimension in Grid Basis
	MatrixXd H(N,N); 
	printf("Construction of Hamiltonian of 1D Harmonic Oscillator\n");
		printf("**************************************************************************************************************************************\n");
	for(int i = -n; i <= n; i++) // x_i = dx * i
	{ 
		for(int j = -n; j <= n; j++) // x_j = dx * j
		{ 
			int di = i < j? j - i: i - j; // for checking the greater than condition and assiging value of 'di' accordingly

			double &matrix_elements = H((i+n) , (j+n)); // Hamiltonian Constructed till now by KE 
			matrix_elements = 0;
			matrix_elements += T_i[di]; // Kinetic Energy Calculation

			if(i == j) // Using the algorithm of Kronecker Delta as mentioned in theory
				H((i+n),(j+n))+= pot_1D(dx * i); // Total Hamiltonian constructed by addition of Kinetic Energy(done previously) + Potential Energy
		}
	}

	printf("Calculating the Eigenvalues\n");
	//Calculating the eigen values 
	SelfAdjointEigenSolver<MatrixXd> es(H);
	VectorXd e_val = es.eigenvalues();
	//Storing all Eigenvalues according to the eigensates
    FILE* output = fopen("eigen_values.txt", "w");
	printf("Energy Eigenvalues are :\n");
	for(int i = 0; i <= state_number; i++)
	{
		printf("%f\n", e_val(i));
		fprintf(output, "%f\n", e_val(i));
	}
    //Storing Eigenstates which have amplitude at the grid points
	FILE* out = fopen("Eigenstates.txt", "w");
	for(int i = 0; i <=state_number; i++)
	{
		VectorXd e_vec = es.eigenvectors().col(i);
		for(int l = -n; l <= n; l++)
		   {
				
				fprintf(out,"%f\n", e_vec((l+n)));
			}
			
	}
	//Storing the values of Potential Energy of Harmonic Oscillator at the grid points same for determining the eigenstates
			FILE* pot= fopen("potential_energy.txt", "w");
	    	for(int i = -n; i <= n; i++)
			{
		
				fprintf(pot, "%f\n", pot_1D(dx * i) );
			}
     	}
