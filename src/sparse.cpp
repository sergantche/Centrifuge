#define _HAS_ITERATOR_DEBUGGING 0
#define _SECURE_SCL 0

#include "sparse.h"

Sparse::Sparse(int s) 
{
	size = s;
	matrix.resize(size);
	for (int i=0; i<size; i++) matrix[i].clear();
	b.resize(size, 0);
};

Sparse::~Sparse()
{
};

//Solve system using SOR method, X - initial guess, omega - relaxation factor
//Requirements all diagonal elements must be non-zero
std::vector<double> Sparse::Solve_by_SOR(std::vector<double>& X, double omega, double eps)
{	
	int k = 0;
	//SOR algorithm
	for (;;) {		
		k++;
		if (k>10000) throw;
	//Calculate next iteration
		for (int i = 0; i<size; i++) {
			double t = 0;
			for (Line::iterator it = matrix[i].begin(); it!=matrix[i].end(); it++) {
				int j = it->first;
				if (j!=i) t+=(it->second)*X[j]; //t += matrix[i][j]*result[j]
			};
			//Calculate result[i]
			if (matrix[i][i] == 0) throw ;
			X[i] = (1-omega)*X[i] + (omega/matrix[i][i])*(b[i] - t);
		};
	//Check if convergence is reached
		double maxR = 0;
		double R;
		for (int i =0; i<size; i++) {
			R = -b[i];
			for (Line::iterator it = matrix[i].begin(); it!=matrix[i].end(); it++) {
				R+=(it->second)*X[it->first];
			};
			R = (R>0)?R:-R;
			maxR = (R>maxR)?R:maxR;
		};
		if (maxR < eps) break;
	};
	return X;
};

void Sparse::Output()
{
	std::map<int, double>::iterator it;
	for (int i=0; i<size; i++) 
	{
		for (int j=0; j<size; j++)
		{
			printf("%.3lf ", matrix[i][j]);
		};
		printf("\n");
	};
	return;
};

double Sparse::get(int i, int j)
{
	std::map<int, double>::const_iterator it;
	if (it == matrix[i].end())
	{
		return 0;
	} else {
		return matrix[i][j];
	};
};

void Sparse::set(int i, int j, double val)
{
	matrix[i][j] = val;
	return;
};

double Sparse::get_b(int i)
{
	return b[i];
};

void Sparse::set_b(int i, double val)
{
	b[i] = val;
	return;
};
