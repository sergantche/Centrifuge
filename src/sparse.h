#define _HAS_ITERATOR_DEBUGGING 0
#define _SECURE_SCL 0

#include <vector>
#include <map>
#include <iostream>

//Represents sparse matrix
class Sparse
{
public:
	//Define hashed matrix line type (position, value)
	typedef std::map<int, double> Line;
	//Size of matrix
	int size;
	//Array of lines , each line represented by hash
	std::vector<Line> matrix;
	//Array of right side values
	std::vector<double> b;
public:
	//Constructor (matrix size)
	Sparse(int);
	//Destructor
	~Sparse();
	//Solving method
	std::vector<double> Solve_by_SOR(std::vector<double>&, double, double);
	//Determinant
	//TO DO double Determinant();
	//Output
	void Output();
	//Accessors for matrix
	double get(int,int);
	void set(int,int, double);
	//Accessors for right side
	double get_b(int);
	void set_b(int, double);
};