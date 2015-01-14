#define _HAS_ITERATOR_DEBUGGING 0
#define _SECURE_SCL 0

#include "grid.h"
#include "sparse.h"
#include "stdio.h"
#include "params.h"
#include "mkl_lapack.h"
#include <math.h>
#include <iostream>

//Cell operation
double max(double x, double y) {
	if (x>y) {
		return x;
	} else {
		return y;
	};
};

double Abs(double x)
{
	return (x>0)?x:-x;
};

//Conditions
Conditions::Conditions()
{
};

Conditions::Conditions(int r_num, int z_num) 
{
	Height = La2/a;
	Radius = 1;
	R_num = r_num;
	Z_num = z_num;
};

Cell Conditions::Init_Conditions(double x, double y)
{
	Cell A;
	A.P = 0;
	A.p = 0;
	A.T = 0;//????
	A.u = 0;
	A.v = 0;
	A.w = 0;
	return A;
};

double Conditions::Temp_Distribution(double z)
{
	//return 0;
	//Linear distribution
	return (5.0/T_0)*(1-(2*a*z)/La2);
};

double Conditions::V_boundary(int i, int j) {
	return 0;
};


//Grid operational))
Grid::Grid(Conditions& C)
{
	cond = C;
	r_num = C.R_num;
	z_num = C.Z_num;
	//Cells
	cells.resize(r_num+2);
	for (int i = 0; i<=r_num+1; i++) cells[i].resize(z_num+2);
	cells_prev.resize(r_num+2);
	for (int i = 0; i<=r_num+1; i++) cells_prev[i].resize(z_num+2);
	r.resize(r_num+2);
	r_left.resize(r_num+2);
	r_right.resize(r_num+2);
	z.resize(z_num+2);
	z_top.resize(z_num+2);
	z_bottom.resize(z_num+2);
	int N = (r_num+2)*(z_num+2) * 5;	
	b.resize(N);
};

Grid::~Grid()
{
};

//Global variabals for computational convenience


//Caclulate ny
double ny(double r) {
	return exp(-A2*(1-(r*r)));//in not dimensionless variables ny=exp(-A2*(1-(r*r/a*a)))
};
double ny_True(double r) {
	return exp(-A2*(1-(r*r/a*a)));
};

//Other coefficients
double U_2 = 2*A2;
double B = 1 + (lambda/my);


void Grid::Step_Back()
{
	for (int i = 0; i<=r_num+1; i++) {
		for (int j = 0; j<=z_num+1; j++) {
			cells[i][j] = cells_prev[i][j];
		};
	};
	return;
};

//Initial conditions + mesh distribution
void Grid::Init()
{
	total_time = 0;
	min_h = (cond.Height>cond.Radius)?(cond.Height):(cond.Radius);
	//radial index i
	//axial index j
	/*
	//Building mesh according to distribution requaried
	for (int i=0; i<=r_num+1; i++) {
		for (int j=0; j<=z_num+1; j++) {
			//Coordinats
			r_left[i] = (cond.Radius/r_num)*(i-1);					
			r_right[i] = r_left[i]+(cond.Radius/r_num);				
			z_bottom[j] = (cond.Height/z_num)*(j-1);		
			z_top[j] = z_bottom[j] + (cond.Height/z_num);	
			r[i] = (r_left[i] + r_right[i])/2;
			z[j] = (z_top[j] + z_bottom[j])/2;
			//Define minimal cell size
			if ((r_right[i] - r_left[i])<min_h) min_h = r_right[i] - r_left[i];
			if ((z_top[j] - z_bottom[j])<min_h) min_h = z_top[j] - z_bottom[j];
		};

	};
	*/
	//For non uniform grid
	//R distribution
	double h = (cond.Radius*(1-R_left)) * (q_r-1)/(pow(q_r, r_num)-1);
	double r_l = R_left;
	r_right[0] = r_l;
	r_left[0] = r_right[0]-h/q_r;
	r[0] = (r_left[0] + r_right[0])/2;
	for (int i = 1; i<=r_num; i++) 
	{		
		r_left[i] = r_l;
		r_right[i] = r_l + h;		
		r[i] = (r_left[i] + r_right[i])/2;
		r_l = r_right[i];
		h *= q_r;
	};
	r_right[r_num] = cond.Radius;
	r_left[r_num+1] = r_right[r_num];
	r_right[r_num+1] = r_left[r_num+1] + h;
	r[r_num+1] = (r_left[r_num+1] + r_right[r_num+1])/2;
	//Z distrinbution (uniform)
	/*
	for (int j=0; j<=z_num+1; j++) {	
			z_bottom[j] = (cond.Height/z_num)*(j-1);		
			z_top[j] = z_bottom[j] + (cond.Height/z_num);				
			z[j] = (z_top[j] + z_bottom[j])/2;
	};
	*/
	double h_z = cond.Height * (q_z-1)/(2*(pow(q_z, z_num/2)-1));
	double z_down = cond.Height/2;
	for (int i = z_num/2+1; i<=z_num+1; i++)
	{
		z_bottom[i] = z_down;
		z_top[i] = z_bottom[i] + h_z;
		z[i] = (z_bottom[i] + z_top[i])/2;
		z_down = z_top[i];
		h_z *= q_z;
	};
	h_z = cond.Height * (q_z-1)/(2*(pow(q_z, z_num/2)-1));
	double z_up = cond.Height/2;
	for (int i = z_num/2; i>=0; i--)
	{
		z_top[i] = z_up;
		z_bottom[i] = z_top[i] - h_z;
		z[i] = (z_bottom[i] + z_top[i])/2;
		z_up = z_bottom[i];
		h_z *= q_z;
	};

	//for (int i = 1; i<=r_num; i++) printf("%.05lf\n", r_right[i] - r_left[i]);
	//exit(0);

	//Initialization of internal cells
	for (int i=0; i<=r_num+1; i++) {
		for (int j=0; j<=z_num+1; j++) {
			//Initialization of cell(i,j)
			//Initial conditions			
			cells[i][j] = cond.Init_Conditions(1,1);///???
			cells[i][j].r_left = r_left[i];
			cells[i][j].r_right = r_right[i];
			cells[i][j].z_bottom = z_bottom[j];
			cells[i][j].z_top = z_top[j];
			cells[i][j].r = r[i];
			cells[i][j].z = z[j];
			cells_prev[i][j] = cells[i][j];
		};
	};

	//Set size for vector of mesh dependent coefficients
	B1.resize(r_num+2);
	B2.resize(r_num+2);
	B20.resize(r_num+2);
	B22.resize(r_num+2);
	B25.resize(r_num+2);
	B32.resize(r_num+2);
	B27.resize(r_num+2);
	B34.resize(r_num+2);
	B5.resize(r_num+2);
	B6.resize(r_num+2);
	B9.resize(r_num+2);
	B10.resize(r_num+2);
	B11.resize(r_num+2);
	B12.resize(r_num+2);
	B3.resize(z_num+2);
	B4.resize(z_num+2);
	B21.resize(z_num+2);
	B7.resize(z_num+2);
	B8.resize(z_num+2);
	B23.resize(z_num+2);
	B24.resize(z_num+2);
	B26.resize(z_num+2);
	B13.resize(z_num+2);
	B14.resize(z_num+2);
	//Calculate mesh dependent coefficients
	for (int i=1; i<=r_num; i++) {
		B1[i] = ((r[i]+r[i+1]))/(r[i]*(r[i+1]-r[i-1])*(r[i+1]-r[i]));
		B2[i] = ((r[i-1]+r[i]))/(r[i]*(r[i+1]-r[i-1])*(r[i]-r[i-1]));
		B20[i] = B1[i]+B2[i]+1/(r[i]*r[i]);
		B22[i] = B1[i] + B2[i];
		B25[i] = r[i]/(Re*ny(r[i]));
		B32[i]=B*2*A2+B5[i]+B6[i]+1/(r_left[i]*r_left[i]);
		B27[i]=1/(2*A2*(r[i]-r[i-1]));
		B34[i]=B*2*A2*r_left[i]/(Re*ny(r_left[i])*(r_left[i+1]-r_left[i-1]));
		B5[i]=r[i]/(r_left[i]*(r[i]-r[i-1])*(r_left[i+1]-r_left[i]));
		B6[i]=r[i-1]/(r_left[i]*(r[i]-r[i-1])*(r_left[i]-r_left[i-1]));
		B9[i]=r_left[i+1]*ny(r_left[i+1])/(ny(r[i])*r[i]*(r_left[i+1]-r_left[i]));
		B10[i]=r_left[i]*ny(r_left[i])/(ny(r[i])*r[i]*(r_left[i+1]-r_left[i]));
		B11[i]=r_left[i+1]*ny(r_left[i+1])*(1/(2*A2*(r[i+1]-r[i])))/(ny(r[i])*r[i]*(r_left[i+1]-r_left[i]));
		B12[i]=r_left[i]*ny(r_left[i])*B27[i]/(ny(r[i])*r[i]*(r_left[i+1]-r_left[i]));
	};
	for (int j=1; j<=z_num; j++) {
		B3[j] = 2/((z[j+1]-z[j-1])*(z[j+1]-z[j]));
		B4[j] = 2/((z[j+1]-z[j-1])*(z[j]-z[j-1]));
		B21[j] = B3[j]+B4[j];
		B7[j] = 1/((z[j] - z[j-1])*(z_top[j] - z_bottom[j]));
		B8[j] = 1/((z[j] - z[j-1])*(z_bottom[j] - z_bottom[j-1]));
		B23[j] = 1/(U_2*(z[j] - z[j-1]));
		B24[j] = A2*B/(z[j]-z[j-1]);
		B26[j] = B7[j]+B8[j];
		B13[j]=(1/(U_2*(z[j+1] - z[j])))/(z_bottom[j+1]-z_bottom[j]);
		B14[j]=(1/(U_2*(z[j+1] - z[j])))/(z_bottom[j+1]-z_bottom[j]);
	};

	//Allocate memory
	int n = (z_num + 2) * (r_num + 2) * 5;
	b.resize(n);
	x.resize(n);
};


//Somehow
void Grid::Out(int n)
{
	char name[] = "srez_centr_000.txt";
	sprintf(name, "srez_centr_%ld.txt", n);
	FILE *fout = fopen(name, "w");    
	for (int i = 0; i<=r_num+1; i++)
	{
		fprintf(fout, "%.5lf %.5lf %.5lf %.5lf %.5lf %.5lf \n", a*cells[i][z_num/2].r, cells[i][z_num/2].T_True, cells[i][z_num/2].u_True, cells[i][z_num/2].w_True, cells[i][z_num/2].P_True, cells[i][z_num/2].p_True);
	};
	fclose(fout);

	double n_up;
	n_up = (int)floor(log((3*pow(q_z, z_num/2)+1)/4)/log(q_z));
	char name_up[] = "srez_up_000.txt";
	sprintf(name_up, "srez_up_%ld.txt", n);
	fout = fopen(name_up, "w");    
	for (int i = 0; i<=r_num+1; i++)
	{
		fprintf(fout, "%.5lf %.5lf %.5lf %.5lf %.5lf %.5lf \n", a*cells[i][z_num/2+n_up].r, cells[i][z_num/2+n_up].T_True, cells[i][z_num/2+n_up].u_True, cells[i][z_num/2+n_up].w_True, cells[i][z_num/2+n_up].P_True, cells[i][z_num/2+n_up].p_True);
	};
	fclose(fout);

	char name_down[] = "srez_down_000.txt";
	sprintf(name_down, "srez_down_%ld.txt", n);
	fout = fopen(name_down, "w");    
	for (int i = 0; i<=r_num+1; i++)
	{
		fprintf(fout, "%.5lf %.5lf %.5lf %.5lf %.5lf %.5lf \n", a*cells[i][z_num/2-n_up+1].r, cells[i][z_num/2-n_up+1].T_True, cells[i][z_num/2-n_up+1].u_True, cells[i][z_num/2-n_up+1].w_True, cells[i][z_num/2-n_up+1].P_True, cells[i][z_num/2-n_up+1].p_True);
	};
	fclose(fout);

};
void Grid::Task_Out()
{
	char name[] = "Task.txt";
	FILE *fout = fopen(name, "w");  
	fprintf(fout, "H = %.1lf sm R = %.1lf sm\n", La2, a);
	fprintf(fout, "R_num = %ld Z_num = %ld \n", r_num, z_num);
	fprintf(fout,"Cv = %lf\n", Cv);
	fprintf(fout,"T_0 = %lf K\n", T_0);
	fprintf(fout,"Om*a = %lf sm/s\n", om_a);
	fprintf(fout,"my = %lf\n", my);
	fprintf(fout,"P_wall = %lf din/sm^2\n", P_wall);
	fprintf(fout,"MV = %ld\n", MV);
	fprintf(fout,"Pr = %lf\n", Pr);
	fprintf(fout,"A2 = %lf\n", A2);
	fprintf(fout,"p_wall = %lf\n", p_wall);
	fprintf(fout,"Re = %lf\n", Re);
	fprintf(fout,"R_left/a = %.2lf\n", R_left);
	double n_up;
	n_up = (int)floor(log((3*pow(q_z, z_num/2)+1)/4)/log(q_z));
	fprintf(fout, "Output_Info\nSrez variables:r, T, u, w, P, ro\n");
	fprintf(fout, "H_up = %.5lf sm\n", a*cells[1][z_num/2+n_up].z);
	fprintf(fout, "H_down = %.5lf sm\n", a*cells[1][z_num/2-n_up+1].z);
	fprintf(fout, "H_centr = %.5lf sm\n", a*cells[1][z_num/2].z);
	fclose(fout);	
};

void Grid::Output(const char *fname)
{
	FILE *fout = fopen(fname, "w");
	fprintf(fout,"TITLE=%s\n","\"Data  TEST\"");
    fprintf(fout,"VARIABLES=\"r\" \"z\" \"u\" \"w\" \"v\" \"P\" \"T\" \"pl\" \n");
	fprintf(fout,"ZONE T=\"  \"\nI=%d J=%d F=POINT\n",z_num+2,r_num+2);
    fprintf(fout,"DT=(SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE)\n");

	for (int i=0;i<=r_num+1;i++)
    {
		for (int j=0;j<=z_num+1;j++)
		{ 
			fprintf(fout,"%13.5lg %13.5lg %13.5lg %13.5lg %13.5lg %13.5lg %13.5lg %13.5lg\n", a*cells[i][j].r, a*cells[i][j].z, cells[i][j].u_True, cells[i][j].w_True, cells[i][j].v_True, cells[i][j].P_True, cells[i][j].T_True, cells[i][j].p_True);
		};
    }; 

	//printf("%lf\n",cells[1][1].r);
	/*
	for (int i = 1; i<=r_num; i++) {
		for (int j = 1; j<=z_num; j++) {
			std::cout<<"("<<cells[i][j].r<<","<<cells[i][j].z<<") ";
		};
		std::cout<<"\n";
	};*/

	fclose(fout);
	return;
};

//Calculate nP and nT from ( i , j )
int Grid::nP(int i, int j) 
{
	return (i*(z_num+2) + j)*5;
};

int Grid::nT(int i, int j)
{
	return (i*(z_num+2) + j)*5 + 1;
};

int Grid::nu(int i, int j)
{
	return (i*(z_num+2) + j)*5 + 2;
};

int Grid::nw(int i, int j)
{
	return (i*(z_num+2) + j)*5 + 3;
};

int Grid::nv(int i, int j)
{
	return (i*(z_num+2) + j)*5 + 4;
};

//Calculate coefficients for systems
double Grid::P1(int i, int j)
{
	return 1.0;
};
double Grid::P2(int i, int j)
{
	return -1.0;
};
double Grid::P3(int i, int j)
{
	return dt*cells[i][j].r_right*ny(cells[i][j].r_right)/
		(cells[i][j].r*ny(cells[i][j].r)*(cells[i][j].r_right-cells[i][j].r_left));
};
double Grid::P4(int i, int j)
{
	return -dt*cells[i][j].r_left*ny(cells[i][j].r_left)/
		(cells[i][j].r*ny(cells[i][j].r)*(cells[i][j].r_right-cells[i][j].r_left));
};
double Grid::P5(int i, int j)
{
	return dt/(cells[i][j].z_top-cells[i][j].z_bottom);
};
double Grid::P6(int i, int j)
{
	return -dt/(cells[i][j].z_top-cells[i][j].z_bottom);
};
double Grid::P0(int i, int j)
{
	return cells[i][j].P-cells[i][j].T;
};

double Grid::V0(int i, int j)
{
	return cells[i][j].v;
};

double Grid::V1(int i, int j)
{
	return 1+dt*(B20[i]+B21[j])/(Re*ny(cells[i][j].r));
};
double Grid::V2(int i, int j)
{
	return dt;
};
double Grid::V3(int i, int j)
{
	return dt;
};
double Grid::V4(int i, int j)
{
	return -dt*B1[i]/(Re*ny(cells[i][j].r));
};
double Grid::V5(int i, int j)
{
	return -dt*B2[i]/(Re*ny(cells[i][j].r));
};
double Grid::V6(int i, int j)
{
	return -dt*B3[j]/(Re*ny(cells[i][j].r));
};
double Grid::V7(int i, int j)
{
	return -dt*B4[j]/(Re*ny(cells[i][j].r));
};

double Grid::T0(int i, int j)
{
	return cells[i][j].T;
};
double Grid::T1(int i, int j)
{
	return 1+Gamma*dt*(B21[j]+B22[i])/(Re*Pr*ny(cells[i][j].r));
};
double Grid::T2(int i, int j)
{
	return -dt*(Gamma-1)*A2*cells[i][j].r;
};
double Grid::T3(int i, int j)
{
	return -dt*(Gamma-1)*A2*cells[i][j].r;
};
double Grid::T4(int i, int j)
{
	return -dt*Gamma*B1[i]/(Re*Pr*ny(cells[i][j].r));
};
double Grid::T5(int i, int j)
{
	return -dt*Gamma*B2[i]/(Re*Pr*ny(cells[i][j].r));
};
double Grid::T6(int i, int j)
{
	return -dt*Gamma*B3[j]/(Re*Pr*ny(cells[i][j].r));
};
double Grid::T7(int i, int j)
{
	return -dt*Gamma*B4[j]/(Re*Pr*ny(cells[i][j].r));
};

double Grid::U0(int i, int j)
{
	return cells[i][j].u;
};
double Grid::U1(int i, int j)
{
	return 1+dt*(B32[i]+B21[j])/(Re*ny(cells[i][j].r_left));
};
double Grid::U2(int i, int j)
{
	return -dt*B27[i];
};
double Grid::U3(int i, int j)
{
	return dt*B27[i];
};
double Grid::U4(int i, int j)
{
	return -dt;
};
double Grid::U5(int i, int j)
{
	return -dt;
};
double Grid::U6(int i, int j)
{
	return cells[i][j].r_left*dt;
};
double Grid::U7(int i, int j)
{
	return cells[i][j].r_left*dt;
};
double Grid::U8(int i, int j)
{
	return dt*(B34[i]-B5[i]/(Re*ny(cells[i][j].r_left)));
};
double Grid::U9(int i, int j)
{
	return dt*(-B34[i]-B6[i]/(Re*ny(cells[i][j].r_left)));
};
double Grid::U10(int i, int j)
{
	return -dt*B3[j]/(Re*ny(cells[i][j].r_left));
};
double Grid::U11(int i, int j)
{
	return -dt*B4[j]/(Re*ny(cells[i][j].r_left));
};
double Grid::W0(int i, int j)
{
	return cells[i][j].w;
};
double Grid::W1(int i, int j)
{
	return 1+dt*(B22[i]+B26[j])/(Re*ny(cells[i][j].r));
};
double Grid::W2(int i, int j)
{
	return -dt*B23[j];
};
double Grid::W3(int i, int j)
{
	return dt*B23[j];
};
double Grid::W4(int i, int j)
{
	return dt*B24[j]*B25[i];
};
double Grid::W5(int i, int j)
{
	return dt*B24[j]*B25[i];
};
double Grid::W6(int i, int j)
{
	return -dt*B24[j]*B25[i];
};
double Grid::W7(int i, int j)
{
	return -dt*B24[j]*B25[i];
};

double Grid::W8(int i, int j)
{
	return -dt*B1[i]/(Re*ny(cells[i][j].r));
};
double Grid::W9(int i, int j)
{
	return -dt*B2[i]/(Re*ny(cells[i][j].r));
};
double Grid::W10(int i, int j)
{
	return -dt*B7[j]/(Re*ny(cells[i][j].r));
};
double Grid::W11(int i, int j)
{
	return -dt*B8[j]/(Re*ny(cells[i][j].r));
};

void Grid::CalcMatrix() {
	//Amount of cells in grid
	int N = (z_num+2) * (r_num+2);	
	Sparse A(5*N);		
	//Construct system matrix
	//Equations for inner cells
	for (int i = 1; i<=r_num; i++) {
		for (int j = 1; j<=z_num; j++) {
			//Coefficients for equation for P
			int line = nP(i,j);
			A.set(line, nP(i,j), P1(i,j));
			A.set(line, nT(i,j), P2(i,j));
			A.set(line, nu(i+1,j), P3(i,j));
			A.set(line, nu(i,j), P4(i,j));
			A.set(line, nw(i,j+1), P5(i,j));
			A.set(line, nw(i,j), P6(i,j));									
			//Coefficients for equation for V
			line = nv(i,j);
			A.set(line, nv(i,j), V1(i,j));
			A.set(line, nu(i,j), V2(i,j));
			A.set(line, nu(i+1,j), V3(i,j));
			A.set(line, nv(i+1,j), V4(i,j));
			A.set(line, nv(i-1,j), V5(i,j));
			A.set(line, nv(i,j+1), V6(i,j));
			A.set(line, nv(i,j-1), V7(i,j));			
			//Coefficients for equation for T
			line = nT(i,j);
			A.set(line, nT(i,j), T1(i,j));
			A.set(line, nu(i+1,j), T2(i,j));
			A.set(line, nu(i,j), T3(i,j));
			A.set(line, nT(i+1,j), T4(i,j));
			A.set(line, nT(i-1,j), T5(i,j));
			A.set(line, nT(i,j+1), T6(i,j));
			A.set(line, nT(i,j-1), T7(i,j));			
			//Coefficients for equation for u
			if (i != 1) {
				line = nu(i,j);
				A.set(line, nu(i,j), U1(i,j));
				A.set(line, nP(i-1,j), U2(i,j));
				A.set(line, nP(i,j), U3(i,j));
				A.set(line, nv(i,j), U4(i,j));
				A.set(line, nv(i-1,j), U5(i,j));
				A.set(line, nT(i,j), U6(i,j));
				A.set(line, nT(i-1,j), U7(i,j));
				A.set(line, nu(i+1,j), U8(i,j));
				A.set(line, nu(i-1,j), U9(i,j));
				A.set(line, nu(i,j+1), U10(i,j));
				A.set(line, nu(i,j-1), U11(i,j));				
			} else {
				A.set(nu(1,j), nu(1,j), 1);
			};
			double tmp = U7(1,1);
			//Coefficients for equation for w
			if (j != 1) {
				line = nw(i,j);
				A.set(line, nw(i,j), W1(i,j));
				A.set(line, nP(i,j-1), W2(i,j));
				A.set(line, nP(i,j), W3(i,j));
				A.set(line, nu(i+1,j), W4(i,j));
				A.set(line, nu(i,j), W5(i,j));
				A.set(line, nu(i+1,j-1), W6(i,j));
				A.set(line, nu(i,j-1), W7(i,j));
				A.set(line, nw(i+1,j), W8(i,j));
				A.set(line, nw(i-1,j), W9(i,j));
				A.set(line, nw(i,j+1), W10(i,j));
				A.set(line, nw(i,j-1), W11(i,j));				
			} else {
				A.set(nw(i,1), nw(i,1), 1);							
			};	
		};
	};
	//Equations for corner cells
	//Just let them be equal zero
	A.set(nP(0,0), nP(0,0), 1);
	A.set(nT(0,0), nT(0,0), 1);
	A.set(nu(0,0), nu(0,0), 1);
	A.set(nw(0,0), nw(0,0), 1);
	A.set(nv(0,0), nv(0,0), 1);

	A.set(nP(0,z_num+1), nP(0,z_num+1), 1);
	A.set(nT(0,z_num+1), nT(0,z_num+1), 1);
	A.set(nu(0,z_num+1), nu(0,z_num+1), 1);
	A.set(nw(0,z_num+1), nw(0,z_num+1), 1);
	A.set(nv(0,z_num+1), nv(0,z_num+1), 1);	

	A.set(nP(r_num+1,0), nP(r_num+1,0), 1);
	A.set(nT(r_num+1,0), nT(r_num+1,0), 1);
	A.set(nu(r_num+1,0), nu(r_num+1,0), 1);
	A.set(nw(r_num+1,0), nw(r_num+1,0), 1);
	A.set(nv(r_num+1,0), nv(r_num+1,0), 1);
	
	A.set(nP(r_num+1,z_num+1), nP(r_num+1,z_num+1), 1);
	A.set(nT(r_num+1,z_num+1), nT(r_num+1,z_num+1), 1);
	A.set(nu(r_num+1,z_num+1), nu(r_num+1,z_num+1), 1);
	A.set(nw(r_num+1,z_num+1), nw(r_num+1,z_num+1), 1);
	A.set(nv(r_num+1,z_num+1), nv(r_num+1,z_num+1), 1);
	
	//Equations for side cells. Applying boundary conditions.
	int i , j;
	//For left edge
	i = 0;
	for (j = 1; j<=z_num; j++) {
		//P(0, j) = P(1, j)
		A.set(nP(0,j), nP(0,j), 1);
		A.set(nP(0,j), nP(1,j), -1);
		//T(0, j) = T(1, j);
		A.set(nT(0,j), nT(0,j), 1);
		A.set(nT(0,j), nT(1,j), -1);
		//
		A.set(nw(0,j), nw(0,j), 1);
		A.set(nw(0,j), nw(1,j), -1);
		//
		A.set(nv(0,j), nv(0,j), 1);
		A.set(nv(0,j), nv(1,j), -1);
		// u = 0
		A.set(nu(i,j), nu(i,j), 1);	
	};
	//For rotor wall
	i = r_num+1;
	for (j = 1; j<=z_num; j++) {
		//P(i, j) = P(i-1, j)
		A.set(nP(i,j), nP(i,j), 1);
		A.set(nP(i,j), nP(i-1,j), -1);
		//T(0, j) = f(z)
		A.set(nT(i,j), nT(i,j), 1);
		// w = u = v = 0
		A.set(nu(i,j), nu(i,j), 1);
		A.set(nw(i,j), nw(i,j), 1);
		A.set(nv(i,j), nv(i,j), 1);
	};
	//For bottom endcap
	j = 0;
	for (i = 1; i<=r_num; i++) {
		//P(i, 0) = P(i, 1)
		A.set(nP(i,0), nP(i,0), 1);
		A.set(nP(i,0), nP(i,1), -1);
		//T(i, 0) = f(0);
		A.set(nT(i,0), nT(i,0), 1);		
		// w = u = v = 0
		A.set(nu(i,j), nu(i,j), 1);		
		// w = 0
		A.set(nw(i,j), nw(i,j), 1);		
		//
		A.set(nv(i,j), nv(i,j), 1);	
	};
	//For top endcap
	j = z_num+1;
	for (i = 1; i<=r_num; i++) {
		//P(i, j) = P(i, j-1)
		A.set(nP(i,j), nP(i,j), 1);
		A.set(nP(i,j), nP(i,j-1), -1);
		//T(i, j) = f(Height);
		A.set(nT(i,j), nT(i,j), 1);		
		// w = u = v = 0
		A.set(nu(i,j), nu(i,j), 1);		
		A.set(nw(i,j), nw(i,j), 1);		
		A.set(nv(i,j), nv(i,j), 1);		
	};		

	//Convert to sparse row format (CSR)
	values.clear();
	columns.clear();
	rowIndex.clear();
	int currentRowIndex = 1; // 1 based indexing
	rowIndex.push_back(currentRowIndex);
	for (int i = 0; i<5*N; i++) {		
		for (std::map<int, double>::iterator it = A.matrix[i].begin(); it!=A.matrix[i].end(); it++) {
				int j = it->first;
				//Mat(i,j) = A.get(i,j);
				values.push_back(A.get(i,j));
				columns.push_back(j+1);
				currentRowIndex++;
			};
		rowIndex.push_back(currentRowIndex);
	};

	//Compute incomplete LU preconditioner
	//this->precond = new ILUT<Matrix>(Mat, 10, 1e-10);

	return;
};

void Grid::CalcRightSide() {
	//Amount of cells in grid	
	for (int i = 0; i < rowIndex.size() - 1; i++) b[i] = 0;
	//Construct system matrix
	//Equations for inner cells
	for (int i = 1; i<=r_num; i++) {
		for (int j = 1; j<=z_num; j++) {
			//Coefficients for equation for P
			int line = nP(i,j);			
			b[line] = P0(i,j);
			//Coefficients for equation for V
			line = nv(i,j);
			b[line] = V0(i,j);
			//Coefficients for equation for T
			line = nT(i,j);			
			b[line] = T0(i,j);
			//Coefficients for equation for u
			if (i != 1) {
				line = nu(i,j);
				b[line] = U0(i,j);
			} else {				
				b[nu(1,j)] = 0;
			};			
			//Coefficients for equation for w
			if (j != 1) {
				line = nw(i,j);
				b[line] = W0(i,j);
			} else {				
				b[nw(i,1)] = 0;			
			};	
		};
	};
	
	//Equations for side cells. Applying boundary conditions.
	int i , j;
	//For left edge
	i = 0;
	for (j = 1; j<=z_num; j++) {		
		b[nu(i,j)] = 0;
	};
	//For rotor wall
	i = r_num+1;
	for (j = 1; j<=z_num; j++) {
		double val = b[nT(i,j)] = cond.Temp_Distribution(cells[i-1][j].z);
		// w = u = v = 0
		b[nu(i,j)] = 0;
		b[nw(i,j)] = 0;
		b[nv(i,j)] = 0;
	};
	//For bottom endcap
	j = 0;
	for (i = 1; i<=r_num; i++) {
		b[nT(i,0)] = cond.Temp_Distribution(0);
		// w = u = v = 0		
		b[nu(i,j)] = 0;
		// w = 0
		b[nw(i,j)] = 0;
		//		
		b[nv(i,j)] = 0;
	};
	//For top endcap
	j = z_num+1;
	for (i = 1; i<=r_num; i++) {		
		b[nT(i,j)] = cond.Temp_Distribution(cond.Height);
		// w = u = v = 0		
		b[nu(i,j)] = 0;		
		b[nw(i,j)] = 0;		
		b[nv(i,j)] = 0;
	};
	// you fill right-hand side and initial guess	

	return;
};

int Grid::Step()
{
	int N = (z_num+2) * (r_num+2);	
	//Solve the system	
	for (int i = 0; i<=r_num+1; i++) {
		for (int j = 0; j<=z_num+1; j++) {			
			x[nP(i,j)] = cells[i][j].P;
			x[nT(i,j)] = cells[i][j].T;
			x[nu(i,j)] = cells[i][j].u;
			x[nw(i,j)] = cells[i][j].w;			
			x[nv(i,j)] = cells[i][j].v;
			cells_prev[i][j] = cells[i][j];
		};
	};
	CalcRightSide();
	int res = SolveSystem();

	//Write acquired values back to cells
	for (int i = 0; i<=r_num+1; i++)
		for (int j = 0; j<=z_num+1; j++)
		{
			cells[i][j].P = x[nP(i,j)];
			cells[i][j].T = x[nT(i,j)];
			cells[i][j].u = x[nu(i,j)];
			cells[i][j].w = x[nw(i,j)];
			cells[i][j].v = x[nv(i,j)];
		};
	return res;
};

int Grid::Run2() {
	double real_dt = a*dt/om_a;
	total_time += real_dt;
	printf("time step: %lf s; total_time: %.5lf s\n", real_dt, total_time);	
	//Calc advanced time values
	return Step();
};

void Grid::FinalCalc()
{
	for(int i=0;i<=r_num+1;i++)
		for(int j=0;j<=z_num+1;j++)
		{
			cells[i][j].u_True = om_a*cells[i][j].u;
			cells[i][j].w_True = om_a*cells[i][j].w;
			cells[i][j].v_True = om_a*cells[i][j].r+om_a*cells[i][j].v;
			cells[i][j].T_True = T_0*(1+cells[i][j].T);
			cells[i][j].P_True = P_wall*ny(cells[i][j].r)*(1+cells[i][j].P);
			cells[i][j].p_True=cells[i][j].P_True*MV/(R*cells[i][j].T_True);
		};
	
};



double Grid::Get_Norm()
{	
	double norm = 0;
	for (int i = 1; i<=r_num; i++)
		for (int j = 1; j<=z_num; j++)
		{
			norm += cells[i][j].u*cells[i][j].u;
			norm += cells[i][j].v*cells[i][j].v;
			norm += cells[i][j].w*cells[i][j].w;
			norm += cells[i][j].P*cells[i][j].P;
			norm += cells[i][j].T*cells[i][j].T;
		};
	return norm;
};

void Grid::Save(const char *name)
{
	std::ofstream fout(name);
	fout<<r_num<<" "<<z_num<<" "<<q_r<<" "<<q_z<<" "<<total_time<<"\n";
	for (int i = 0; i<=r_num+1; i++)
	{
		for (int j = 0; j<=z_num+1; j++)
		{
			fout<<cells[i][j].u<<" "<<cells[i][j].w<<" "<<cells[i][j].v<<" "<<cells[i][j].P<<" ";
			fout<<cells[i][j].T<<" "<<cells[i][j].p<<" ";
			fout<<"\n";
		};
	};
	return;
};

void Grid::Load(const char *name)
{
	std::ifstream fin(name);	
	double _total_time;
	fin>>r_num>>z_num>>q_r>>q_z>>_total_time;
	//Cells
	cells.resize(r_num+2);
	for (int i = 0; i<=r_num+1; i++) cells[i].resize(z_num+2);
	r.resize(r_num+2);
	r_left.resize(r_num+2);
	r_right.resize(r_num+2);
	z.resize(z_num+2);
	z_top.resize(z_num+2);
	z_bottom.resize(z_num+2);
	Init();
	for (int i = 0; i<=r_num+1; i++)
	{
		for (int j = 0; j<=z_num+1; j++)
		{
			fin>>cells[i][j].u>>cells[i][j].w>>cells[i][j].v>>cells[i][j].P;
			fin>>cells[i][j].T>>cells[i][j].p;
		};
	};	
	total_time = _total_time;
	return;
};
double Grid::Mass()
{
	double mass = 0;
	for (int i = 1; i<=r_num; i++)
	{
		for (int j = 1; j<=z_num; j++)		
		mass+=(cells[i][j].r_right-cells[i][j].r_left)*(cells[i][j].z_top-cells[i][j].z_bottom)*cells[i][j].p_True;
		
	};
	return mass;

};


//Solve the system	
int Grid::SolveSystem() {	
	// ==== Pardiso parameters ====	
	//1.  Internal solver memory pointer pt
	_MKL_DSS_HANDLE_t pt[64];
	for (int i = 0; i<64; i++) pt[i] = 0;

	//2. Maximal number of factors in memory
	long long maxfct = 1;

	//3. The number of matrix (from 1 to maxfct) to solve
	long long mnum = 1;

	//4. Matrix type
	long long mtype = 11; // Real and nonsymmetric matrix

	//5. Controls the execution of the solver
	long long phase;

	//6. Number of right hand sides.
	long long nrhs = 1;

	/* Pardiso control parameters. */
	long long iparm[64];	
	for (int i = 0; i < 64; i++) {
		iparm[i] = 0;
	}
	
	iparm[0] = 1; /* No solver default */
	iparm[1] = 2; /* Fill-in reordering from METIS */
	/* Numbers of processors, value of OMP_NUM_THREADS */
	iparm[2] = 1;
	iparm[3] = 0; /* No iterative-direct algorithm */
	iparm[4] = 0; /* No user fill-in reducing permutation */
	iparm[5] = 0; /* Write solution into x */
	iparm[6] = 0; /* Not in use */
	iparm[7] = 2; /* Max numbers of iterative refinement steps */
	iparm[8] = 0; /* Not in use */
	iparm[9] = 13; /* Perturb the pivot elements with 1E-13 */
	iparm[10] = 1; /* Use nonsymmetric permutation and scaling MPS */
	iparm[11] = 0; /* Not in use */
	iparm[12] = 1; /* Maximum weighted matching algorithm is switched-on (default for non-symmetric) */
	iparm[13] = 0; /* Output: Number of perturbed pivots */
	iparm[14] = 0; /* Not in use */
	iparm[15] = 0; /* Not in use */
	iparm[16] = 0; /* Not in use */
	iparm[17] = -1; /* Output: Number of nonzeros in the factor LU */
	iparm[18] = -1; /* Output: Mflops for LU factorization */
	iparm[19] = 0; /* Output: Numbers of CG Iterations */	

	long long msglvl = 1; /* Print statistical information in file */
	long long error = 0; /* Initialize error flag */
	
	long long n = rowIndex.size() - 1; // Number of equations

	long long* perm = new long long[n];

	/* -------------------------------------------------------------------- */
	/* .. Reordering and Symbolic Factorization. This step also allocates */
	/* all memory that is necessary for the factorization. */
	/* -------------------------------------------------------------------- */
	phase = 11;
	pardiso_64(pt, &maxfct, &mnum, &mtype, &phase,
		&n, values.data(), rowIndex.data(), columns.data(), perm, &nrhs, iparm, &msglvl, b.data(), x.data(), &error);
	if (error != 0) {
		printf("\nERROR during symbolic factorization: %d", error);
		exit(1);
	}
	printf("\nReordering completed ... ");
	printf("\nNumber of nonzeros in factors = %d", iparm[17]);
	printf("\nNumber of factorization MFLOPS = %d", iparm[18]);

	/* -------------------------------------------------------------------- */
	/* .. Numerical factorization. */
	/* -------------------------------------------------------------------- */
	phase = 22;
	pardiso_64(pt, &maxfct, &mnum, &mtype, &phase,
		&n, values.data(), rowIndex.data(), columns.data(), perm, &nrhs, iparm, &msglvl, b.data(), x.data(), &error);
	if (error != 0) {
		printf("\nERROR during numerical factorization: %d", error);
		exit(2);
	}
	printf("\nFactorization completed ... ");

	/* -------------------------------------------------------------------- */
	/* .. Back substitution and iterative refinement. */
	/* -------------------------------------------------------------------- */
	phase = 33;
	iparm[7] = 2; /* Max numbers of iterative refinement steps. */		
	pardiso_64(pt, &maxfct, &mnum, &mtype, &phase,
		&n, values.data(), rowIndex.data(), columns.data(), perm, &nrhs, iparm, &msglvl, b.data(), x.data(), &error);
	if (error != 0) {
		printf("\nERROR during solution: %d", error);
		exit(3);
	}
	printf("\nSolve completed ... ");
	printf("\nThe solution of the system is: ");
	for (int i = 0; i < n; i++) {
		//printf("\n x [%d] = % f", i, x[i] );
	}
	printf ("\n");
	/* -------------------------------------------------------------------- */
	/* .. Termination and release of memory. */
	/* -------------------------------------------------------------------- */
	phase = -1; /* Release internal memory. */
	pardiso_64(pt, &maxfct, &mnum, &mtype, &phase,
		&n, values.data(), rowIndex.data(), columns.data(), perm, &nrhs, iparm, &msglvl, b.data(), x.data(), &error);	

	return error;
};

double Grid::CalcConditionNumber() {
	double estimate = 0;
	//dgecon("1", 
	return estimate;
};