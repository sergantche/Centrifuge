#pragma once

#include <vector>
#include <iostream>
#include <fstream>
#include "mkl.h"
#include "mkl_dss.h" 

typedef double Type;

double max(double , double );

double Abs(double);

//Cell object class
class Cell
{
public:
	//Pressure
	double P;
	//Density
	double p;
	//Velocity
	double u;	//Radial
	double v;	//Tangential
	double w;	//Axial
	//Temperature
	double T;
	//Not undemenshional Pressure
	double P_True;
	//Not undemenshional Density
	double p_True;
	//Velocity
	double u_True;	//Not undemenshional Radial
	double v_True;	//Not undemenshional Tangential
	double w_True;	//Not undemenshional Axial
	//Not undemenshional Temperature
	double T_True;
	//Coordinats
	double r_left;
	double r_right;
	double z_top;
	double z_bottom;
	double r;
	double z;
};

//Class that contains : initial conditions, axial and radial
//mesh cells distribution, parametrs of cylinder, etc.
class Conditions {
public:
	//Parametrs of cylinder and number of cells in each direction
	double Radius;
	int R_num;
	double Height;
	int Z_num;
	//Distribution of cells
	double Axial_Cell_Density(double);
	double Radial_Cell_Density(double);
	//Initial conditions
	Cell Init_Conditions(double, double);
	//Rotor wall temperature distribution
	double Temp_Distribution(double);
	//Boundary conditions values for fictitious cells
	double V_boundary(int, int);
	//Constructor
	Conditions();
	Conditions(int, int);
};

//Class represents grid object
class Grid 
{
public:
	//Size of grid matrix and coefficient of geometric progression
	double q_r;
	double q_z;
	int r_num;
	int z_num;
	//Task describing class instance
	Conditions cond;
	//Array of cells
	std::vector<std::vector<Cell> > cells;
	std::vector<std::vector<Cell> > cells_prev;
	std::vector<double> r_left;
	std::vector<double> r_right;
	std::vector<double> z_top;
	std::vector<double> z_bottom;
	std::vector<double> r;
	std::vector<double> z;
	//Time step
	double dt;
	//Size of mesh
	double min_h;
	//Matrix
	//Matrix Mat;

	//Matrix in CSR format
	std::vector<double> values;
	std::vector<long long> columns;
	std::vector<long long> rowIndex;
	void CalcMatrix();
	double CalcConditionNumber();
	//Vector of variables
	std::vector<double> x;
	//Right side
	std::vector<double> b;
	void CalcRightSide();
	//Preconditioner
	//ILUT<Matrix>* precond;
	void Step_Back();
	int Step();
	int Run2();
public:
	//Total work time
	double total_time;
	//Constructor
	Grid(Conditions& C);
	//Destructor
	~Grid();
	//Init grid cells, initial conditions etc.
	void Init();	

	//Solve system
	int SolveSystem();

	//Calculations
	double eps;		
	//Calculate norm
	double Get_Norm();
	//For not demensionless variables
	void FinalCalc();
	//Output results somehow
	void Out(int);
	void Output(const char*);
	void Task_Out();
	//Save grid to file
	void Save(const char *);
	//Load grid from file
	void Load(const char *);
	double Mass();
	//Calculate position of (i,j) element in vector (for system) and vice versa
	int nP(int, int);
	int nT(int, int);
	int nw(int, int);
	int nu(int, int);
	int nv(int, int);
	
	//Mesh dependent coefficients
	std::vector<double> B1;
	std::vector<double> B2;
	std::vector<double> B3;
	std::vector<double> B4;
	std::vector<double> B5;
	std::vector<double> B6;
	std::vector<double> B7;
	std::vector<double> B8;
	std::vector<double> B9;
	std::vector<double> B10;
	std::vector<double> B11;
	std::vector<double> B12;
	std::vector<double> B13;
	std::vector<double> B14;
	std::vector<double> B20;
	std::vector<double> B21;
	std::vector<double> B22;
	std::vector<double> B23;
	std::vector<double> B24;
	std::vector<double> B25;
	std::vector<double> B26;
	std::vector<double> B27;
	std::vector<double> B32;
	std::vector<double> B34;
	//
	double P1(int, int);
	double P2(int, int);
	double P3(int, int);
	double P4(int, int);
	double P5(int, int);
	double P6(int, int);
	double P0(int, int);

	double V1(int, int);
	double V2(int, int);
	double V3(int, int);
	double V4(int, int);
	double V5(int, int);
	double V6(int, int);
	double V7(int, int);
	double V0(int, int);

	double T0(int, int);
	double T1(int, int);
	double T2(int, int);
	double T3(int, int);
	double T4(int, int);
	double T5(int, int);
	double T6(int, int);
	double T7(int, int);

	double U0(int, int);
	double U1(int, int);
	double U2(int, int);
	double U3(int, int);
	double U4(int, int);
	double U5(int, int);
	double U6(int, int);
	double U7(int, int);
	double U8(int, int);
	double U9(int, int);
	double U10(int, int);
	double U11(int, int);

	double W0(int, int);
	double W1(int, int);
	double W2(int, int);
	double W3(int, int);
	double W4(int, int);
	double W5(int, int);
	double W6(int, int);
	double W7(int, int);
	double W8(int, int);
	double W9(int, int);
	double W10(int, int);
	double W11(int, int);
};

