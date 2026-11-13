#include "udf.h"
#include "math.h"

DEFINE_MASS_TRANSFER(co2_mass_transfer, c, mixture_thread, from_phase_index, from_species_index, to_phase_index, to_species_index)
{
    /* Constants */
    double db = 0.00286;     // Bubble diameter [m]
    double D_L = 1.4663e-9;  // CO2 diffusivity in water [m^2/s]
    double mu_L = 1.1676e-3; // Dynamic viscosity of liquid [Pa·s]
    double rho_L = 1000.0;   // Density of liquid [kg/m^3]
    double He = 1.0432;      // Henry's constant (dimensionless)

    /* Threads */
    Thread *gas_thread = THREAD_SUB_THREAD(mixture_thread, from_phase_index);
    Thread *liquid_thread = THREAD_SUB_THREAD(mixture_thread, to_phase_index);

    /* Local properties */
    double aG = C_VOF(c, gas_thread);                      // Gas holdup
    double rhoG = C_R(c, gas_thread);                      // Gas density [kg/m^3]
    double X_CO2_G = C_YI(c, gas_thread, 0);               // CO2 mass fraction in gas
    double X_CO2_L = C_YI(c, liquid_thread, 0);            // CO2 mass fraction in liquid
    double uG = C_U(c, gas_thread);                        // Gas velocity [m/s]
    double uL = C_U(c, liquid_thread);                     // Liquid velocity [m/s]

    /* Relative velocity and dimensionless numbers */
    double u_rel = fabs(uG - uL);                          // Relative velocity
    double nu_L = mu_L / rho_L;                            // Kinematic viscosity [m²/s]
    double Sc = nu_L / D_L;                                // Schmidt number
    double Re = u_rel * db / nu_L;                         // Reynolds number
    double Pe = Re * Sc;                                   // Peclet number

    /* Sherwood number - Eq. 5.11 from Hissanaga (correct form) */
    double Sh = 1.2;
    if (Re > 0.0)
    {
        double correction = 1.0 - 2.828 / sqrt(Re);
        if (correction > 0.0)
            Sh += sqrt(0.64 + 1.273 * correction * Pe);
    }

    /* Mass transfer coefficient kL [m/s] */
    double kL = Sh * D_L / db;

    /* Specific interfacial area a [1/m] */
    double a = 6.0 * aG / db;

    /* Mass transfer source term [kg/m^3/s] */
    double mass_transfer = kL * a * (He * rhoG * X_CO2_G - rho_L * X_CO2_L);

    return mass_transfer;
}
