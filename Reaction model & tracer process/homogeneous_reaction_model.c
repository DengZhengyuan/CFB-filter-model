// -------------------------------------------------------
//     UDS_0 -> %mass of ozone in gas phase
//     X UDS_1 -> solids holdup in mixture phase
//     X UDS_2 -> slip velocity in mixture phase
// -------------------------------------------------------
//     UDM_0 -> slip velocity, m/s
//     UDM_1 -> reaction rate, kg/m^3/s
//     UDM_2 -> ozone concentration, kg/m^3/s
// --------------
//     X UDM_3 -> coefficient of Hr
//     X UDM_4 -> coefficient of Hd
// -------------------------------------------------------

#define diam_p 70.e-6        // particle diameter [m]
#define oz_conc 1.6653e-4    // initial ozone concentration, 100 ppm -> %mass
#define k_v 49.2             // reaction constant based on catalysts volume, [/s]
#define D_ozone 1.48535e-5   // mass diffusivity of ozone in air [m^2/s]
#define density_gas 1.225    // density of air, [kg/m^3]
#define density_solids 1780. // density of solids, [kg/m^3]

#include "udf.h"
#include "stdio.h"
#include "math.h"

// -------------------------------------------------------------
//                    reaction term of ozone
// -------------------------------------------------------------
DEFINE_SOURCE(rxn_ozone, cell, thread, dS, eqn)
{
    real rxn_rate, Y_ozone, coef_Hr, k;
    real epsilon_g, epsilon_s, rho_g;

    Y_ozone = C_UDSI(cell, thread, 0); // mole of ozone in the cell
    epsilon_g = C_VOF(cell, thread);   // volume fraction of gas
    epsilon_s = 1. - epsilon_g;        // volume fraction of solids
    rho_g = C_R(cell, thread);         // gas density

    //---- calculate the reaction rate ----
    coef_Hr = 1.;
    k = coef_Hr * k_v;
    rxn_rate = -k * rho_g * epsilon_s * Y_ozone;
    dS[eqn] = -k * rho_g * epsilon_s;
    C_UDMI(cell, thread, 1) = rxn_rate;
    C_UDMI(cell, thread, 2) = Y_ozone / oz_conc;
    // C_UDMI(cell, thread, 3) = coef_Hr;

    return rxn_rate;
}

// -------------------------------------------------------------
//                     diffusivity of ozone
// -------------------------------------------------------------
DEFINE_DIFFUSIVITY(diff_ozone_laminar, cell, thread, i)
{
    Thread *thread_mix, *thread_solids, *thread_gas;

    real diff, D_AB;
    real rho_g, u_slip;
    real epsilon_g, epsilon_s;
    real Y_oz_g;

    //---------------------------------
    //--------- slip velocity ---------
    // find the solids thread
    thread_mix = THREAD_SUPER_THREAD(thread);
    thread_solids = THREAD_SUB_THREAD(thread_mix, 1);
    thread_gas = THREAD_SUB_THREAD(thread_mix, 0);

#if RP_2D
    u_slip = pow(
        pow(C_U(cell, thread) - C_U(cell, thread_solids), 2.) +
            pow(C_V(cell, thread) - C_V(cell, thread_solids), 2.),
        .5);
#endif

#if RP_3D
    u_slip = pow(
        pow(C_U(cell, thread) - C_U(cell, thread_solids), 2.) +
            pow(C_V(cell, thread) - C_V(cell, thread_solids), 2.) +
            pow(C_W(cell, thread) - C_W(cell, thread_solids), 2.),
        .5);
#endif
    // UDM_0 -> slip velocity
    C_UDMI(cell, thread, 0) = u_slip;

    //--------- slip velocity ---------
    //---------------------------------

    // volume fraction of solids and gas
    epsilon_g = C_VOF(cell, thread);
    epsilon_s = 1. - epsilon_g;
    Y_oz_g = C_UDSI(cell, thread_gas, 0); // %mass of ozone in the gas phase

    //---------------------------------
    // C_UDSI(cell, thread_mix, 1) = epsilon_s;
    // C_UDSI(cell, thread_mix, 2) = u_slip;
    //---------------------------------

    // define some parameters and properties
    rho_g = C_R(cell, thread);
    D_AB = D_ozone * (1. - pow(epsilon_s, 0.5)) / epsilon_g;

    if (Y_oz_g < (0.01 * oz_conc))
    {
        diff = 0.;
    }
    else
    {
        diff = rho_g * D_AB;
    }

    return diff;
}
