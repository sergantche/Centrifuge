#ifndef PARAMS
#define _HAS_ITERATOR_DEBUGGING 0
#define _SECURE_SCL 0

//This is a not SGS system!!!

//SGS system!!!

double a = 10;//cm
double La2 = 40;//cm
double Cv = 3.65E6;///
double T_0 = 320;//K
double om_a = 2.6E4;//cm/s
double my = 1.47E-4;
double P_wall =1.995E5; //din/cm^2 (150torr)
double MV=349;//Molecular Weight
double R=8.31E7;//Universal Gas Constant
double Gamma = 1.065;
double Pr = 0.7;
double I_0 = Cv*T_0;
double A2 = (om_a*om_a)/(2*(Gamma - 1)*I_0);//23.613;
double p_wall=P_wall/((Gamma - 1)*I_0);//gr/cm^3
double Re = om_a*a*p_wall/my;
double lambda = -2*my/3.0;//0.1;
double R_left=0.5;



/*
double a = 10;//cm
double La2 = 40;//cm
double A2 = 1.855;//23.613;
double Re = 10.66E4;
double Gamma = 1.0935;
double Pr = 0.7;
double T_0 = 320;//K
double om_a = 6E4;//6E4;//cm/s
double p_wall=1.097E-3;//2.617E-3;//gr/cm^3
double P_wall =1.995E5; //din/cm^2 (150torr)
double my = 1.47E-4;
double lambda = -2*my/3.0;//0.1;
double MV=146;//349;//Molecular Weight
double R=8.31E7;//Universal Gas Constant
*/

#define PARAMS
#endif


