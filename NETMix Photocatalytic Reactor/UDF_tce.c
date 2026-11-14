#include "udf.h"

#define k_rate 3.8e-11 // kmol m^2/s * W
#define K_value 7.0e3 // m^3/kmol


#define mw_tce 131.39 // kg/kmol
#define mw_h2o 18.01528 // kg/kmol
#define mw_air 28.9644 //kg/kmol

#define Rgi 8314.462 // Pa*m3/K/kmol

DEFINE_SR_RATE(UDF_tce, f, t, r, mw, yi, rr)
{
#if !RP_HOST 

    Thread *t0 = t->t0;
    cell_t c0 = F_C0(f, t);

	// TCE mass fraction
    double w_tce = yi[0]; 
	// Water mass fraction
    double w_h2o = yi[1]; 
	// Air mass fraction
    double w_air = yi[2]; 

    // Irradiance
    double I = F_UDSI(f, t, 0);

    // Normalized sum of mass fractions
    double Wsum = w_tce + w_h2o + w_air;

    // Normalized mass fractions of each component
    double x_tce = w_tce / Wsum;
    double x_h2o = w_h2o / Wsum;
    double x_air = w_air / Wsum;

    // Molar fraction calculation
    double Nsum = x_tce / mw_tce + x_air / mw_air + x_h2o / mw_h2o;
    double y_tce = (x_tce / mw_tce) / Nsum;
    double y_h2o = (x_h2o / mw_h2o) / Nsum;
    double y_air = (x_air / mw_air) / Nsum;

    // Reaction rate calculation
    double T_w = F_T(f, t); //K
    double p_tce = 101325. * y_tce; //Pa
    double p_h2o = 101325. * y_h2o; //Pa

    double C_tce = p_tce / (Rgi*T_w); //kmol/m3
    double C_h2o = p_h2o / (Rgi*T_w); //kmol/m3

    *rr = (I * k_rate * K_value * C_tce) / (1.0 + K_value * C_tce);

#endif // !RP_NODE
}
