#define _HAS_ITERATOR_DEBUGGING 0
#define _SECURE_SCL 0

#include "stdio.h"
#include "grid.h"
#include "sparse.h"
#include <math.h>

Conditions conditions(60,30);
Grid grid(conditions);

double define_timestep() {
	//define timestep
	int m = 0;
	for (;;) {
		grid.Step_Back();
		grid.CalcMatrix();
		if (grid.Step() != 0) {
			grid.dt /= 2.0;
		} else {
			if (m<5) {
				grid.dt *= 1.5;
				m++;
			} else {
				break;
			};
		};
	};
	double real_dt = grid.dt;
	printf("Finally step is: %lf\n", real_dt);
	FILE *fout = fopen("timesteps.txt", "a");
	fprintf(fout,"total time: %lf , dt = %lf\n", grid.total_time, real_dt);
	fclose(fout);
	return real_dt;
};

int main(int argc, char *argv[])
{	
	int k = 1;
	//Define calculation parameters
	double omega = 0.01;
	grid.dt = 10.0;
	grid.eps = 1E-4;
	grid.q_r = 0.97;
	grid.q_z = 0.97;
	grid.Init();
	grid.Task_Out();
	define_timestep();	
	//grid.Load("65.sav");
	k = (int)floor(grid.total_time*10);
	double old_norm = grid.Get_Norm();		
	int change_step = 0;
	char name[] = "Norm.txt";
	for (int i=0; grid.total_time<1000; i++) {		
		if (grid.Run2() != 0) {
			define_timestep();
			change_step = 0;
		} else change_step++;
		if (change_step == 1000) {
			define_timestep();
			change_step=0;
		};

		int time = (int)floor(grid.total_time*10);
		if (time == k) {
			old_norm = grid.Get_Norm();
			FILE *fout = fopen(name, "a");
			fprintf(fout, "%lf %lf\n", grid.total_time, old_norm);
			fclose(fout);
			grid.FinalCalc();
			char fname[] = "000.txt";
			sprintf(fname, "%ld.txt", k);			
			grid.Output(fname);
			sprintf(fname, "%ld.sav", k);
			grid.Save(fname);
			grid.Out(k);
			define_timestep();
			k++;
		};
		if (i == 15) omega = 0.01;
		if (i == 10000) {
			grid.FinalCalc();
			grid.Output("test.txt");
		};
		printf("%ld\n",i);
	};

	return 0;
};
