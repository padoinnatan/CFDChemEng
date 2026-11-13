#include "udf.h"

#define mw_no    30.0061   //kg/kmol
#define mw_h2o   18.01534  //kg/kmol
#define mw_no2   46.00550  //kg/kmol
#define mw_air   28.9644   //kg/kmol
#define H        0.003     //m
#define Rg       8314.462  //Pa*m³/K/kmol
#define E        10        //W/m²

DEFINE_SR_RATE(my_rate, f, t, r, mw, yi, rr)
{

	Thread *t0 = t->t0;
	cell_t c0  = F_C0(f,t);
	
	double w_no  = yi[0];
	double w_h2o = yi[1];
	double w_no2 = yi[2];
	double w_air = yi[3];
	
	double Nsum = w_no/mw_no + w_h2o/mw_h2o + w_no2/mw_no2 + w_air/mw_air;
	
	double y_no  = (w_no/mw_no)/Nsum;
	double y_h2o = (w_h2o/mw_h2o)/Nsum;
	double y_no2 = (w_no2/mw_no2)/Nsum;
	double y_air = (w_air/mw_air)/Nsum;

	double T_w = F_T(f,t);           //K
	double P_w = F_P(f,t) + 101325;  //Pa
	
	double p_no  = P_w*y_no;   //Pa
	double p_h2o = P_w*y_h2o;  //Pa
	double p_no2 = P_w*y_no2;  //Pa
	double p_air = P_w*y_air;  //Pa
	
	double C_no   = p_no/(Rg*T_w);   //kmol/m³
	double C_h2o  = p_h2o/(Rg*T_w);  //kmol/m³
	double C_no2  = p_no2/(Rg*T_w);  //kmol/m³
	double C_air  = p_air/(Rg*T_w);  //kmol/m³
	
	double fNO  = 0.9722;
	double fNO2 = 1.0125;
	
	double k_no  = 4.18/fNO;     //m/s
	double k_no2 = 6.73e1/fNO2;  //m/s

	double K_no  = 8.48e8/fNO;   //m³/kmol
	double K_no2 = 3.02e8/fNO2;  //m³/kmol
	double K_h2o = 5.07e4;       //m³/kmol

	double alfa = 2.37e-3;       //m²/W
	
	double DEN = 1 + (K_no*C_no) + (K_no2*C_no2) + (K_h2o*C_h2o);
	
	if (STREQ(r->name, "reaction-1"))
		*rr = ((k_no*C_no)*(-1+sqrt(1+(alfa*E)))/DEN);                 //kmol/m²/s
	else if (STREQ(r->name, "reaction-2"))
		*rr = (((k_no2*C_no2)-(k_no*C_no))*(-1+sqrt(1+(alfa*E)))/DEN); //kmol/m²/s
	
}