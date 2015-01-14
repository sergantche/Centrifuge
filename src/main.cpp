#define _HAS_ITERATOR_DEBUGGING 0
#define _SECURE_SCL 0

#include "stdio.h"
#include "grid.h"
#include "sparse.h"
#include <iostream>
#include <sstream>
#include <math.h>

int main(int argc, char *argv[])
{		
	

	int k = 1;
	//Define calculation parameters
	Conditions conditions(60, 40);
	Grid grid(conditions);

	//Debug
	//int dn = 3;
	//Sparse A(dn);
	//A.set(0,0, 1);
	//A.set(1,1, 1);
	//A.set(2,2, 1);
	//grid.b.resize(dn, 0);
	//grid.b[1] = 15;
	//grid.x.resize(dn, 0);

	//grid.values.clear();
	//grid.columns.clear();
	//grid.rowIndex.clear();
	//int currentRowIndex = 1; // 1 based indexing
	//grid.rowIndex.push_back(currentRowIndex);
	//for (int i = 0; i< dn; i++) {		
	//	for (std::map<int, double>::iterator it = A.matrix[i].begin(); it!=A.matrix[i].end(); it++) {
	//			int j = it->first;
	//			//Mat(i,j) = A.get(i,j);
	//			grid.values.push_back(A.get(i,j));
	//			grid.columns.push_back(j+1);
	//			currentRowIndex++;
	//		};
	//	grid.rowIndex.push_back(currentRowIndex);
	//};
	//grid.SolveSystem();

	double omega = 0.01;
	grid.dt = 10.0;
	grid.eps = 1E-4;
	grid.q_r = 0.97;
	grid.q_z = 0.97;	
	grid.Init();	

	//Output initial data
	grid.FinalCalc();
	grid.Output("init.dat");
	grid.Task_Out();

	//Define timestep
	

	//Parameters
	const int MaxIter = 1000;
	const int SolutionSnapshotIter = 1;
	const double toleranceEps = 0.001;

	//Simulation main cycle
	std::ofstream historyFile("convergence.dat");
	double norm = grid.Get_Norm();
	double prevNorm = norm;
	for (int iter = 0; iter < MaxIter; iter++) {
		//Define timestep
		grid.dt = 1;

		//Run main step
		grid.CalcMatrix();
		/*int n = grid.rowIndex.size() - 1;
		std::cout<<"x = zeros("<<n<<", 1);\n";
		std::cout<<"A = sparse("<<n<<");\n";
		for (int i = 0; i < n; i++) {
			for (int j = grid.rowIndex[i]; j < grid.rowIndex[i+1]; j++) {
				std::cout<<"A("<<i+1<<","<<grid.columns[j-1]<<") = "<<grid.values[j-1]<<";\n";
			};
		};
		grid.CalcRightSide();
		std::cout<<"b = zeros("<<n<<", 1);\n";
		for (int i = 0; i < n; i++) {
			std::cout<<"b("<<i+1<<") = "<<grid.b[i]<<";\n";
		};*/

		int res = grid.Step();

		//Output snapshots
		if (iter % SolutionSnapshotIter == 0) {
			std::stringstream filename((std::string()));
			filename << "I" << iter << ".dat";
			grid.FinalCalc();
			grid.Output(filename.str().c_str());
		};

		//Compute norm
		prevNorm = norm;
		norm = grid.Get_Norm();
		double relativeNormChange = std::abs(norm - prevNorm)/norm;

		//Save convergence history
		historyFile<<iter<<" "<<relativeNormChange<<std::endl;
		
		//Stop criteria
		if (std::abs(norm - prevNorm) <= toleranceEps * norm) {
			std::cout<<"Stop criteria reached"<<std::endl;
			std::cout<<"Stopping iterations at relative norm change "<<relativeNormChange<<std::endl;
			break;
		};		
	};	

	//Close opened files
	historyFile.close();
	return 0;
};
