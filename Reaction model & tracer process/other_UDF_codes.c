// -------------------------------------------------------
//     X UDM_4 -> coefficient of Hd
// -------------------------------------------------------

#define diam_p 70.e-6 // particle diameter [m]

#include "udf.h"
#include "stdio.h"
#include "math.h"

// -------------------------------------------------------------
//                      Sherwood number
// -------------------------------------------------------------

double Sh_num_Scala2013(double U_slip, double e_s)
{
    double mu_g = 1.79e-5;
    double rho_g = 1.225;
    double tau = 1.2, b = .7, c = .5, d = 1. / 3., X_ab = 1.;
    double e_g, Re_p, Sc, Sh_local, Sh, u_g, kd;
    if (e_s < 1.e-6 || U_slip < 1.e-6)
    {
        kd = 0.;
    }
    else
    {
        e_g = 1. - e_s;
        u_g = U_slip * e_g;
        Re_p = (rho_g * u_g * diam_p) / mu_g;
        Sc = mu_g / (rho_g * D_ozone);
        Sh_local = (2. * e_g / tau) / (1. - pow(e_s * X_ab, 1. / 3.)) + b * pow(Re_p / e_g, c) * pow(Sc, d);
        Sh = (pow(Re_p, 2.) * pow(Sc, 2.)) / (12. * e_s * X_ab * e_g / tau) * (pow(1. + (24. * e_s * X_ab * e_g / tau) / (pow(Re_p, 2.) * pow(Sc, 2.)) * Sh_local, .5) - 1.);
        kd = Sh * D_ozone / diam_p;
    }
    return kd;
}

// -------------------------------------------------------------
//                      Gidaspow drag + Hd
// -------------------------------------------------------------
DEFINE_EXCHANGE_PROPERTY(Gidaspow_Hd_drag, cell, thread_mix, s_phase, f_phase)
{
    Thread *thread_g, *thread_s;

    real rho_g, rho_s, mu_g, epsilon_g, epsilon_s, cell_vol, u_slip;
    real C_D, num_Re, k_g_s_HG, coef_Hd, k_g_s;

    thread_g = THREAD_SUB_THREAD(thread_mix, s_phase); // gas phase
    thread_s = THREAD_SUB_THREAD(thread_mix, f_phase); // solid phase

    rho_g = C_R(cell, thread_g);           // density of gas
    rho_s = C_R(cell, thread_s);           // density of solid
    mu_g = C_MU_L(cell, thread_g);         // viscosity of gas
    epsilon_g = C_VOF(cell, thread_g);     // voidage
    epsilon_s = 1. - epsilon_g + 1.e-8;    // solids holdup
    cell_vol = C_VOLUME(cell, thread_mix); // cell volume
    u_slip = C_UDMI(cell, thread_g, 0);    // slip velocity

    num_Re = rho_g * u_slip * diam_p / mu_g; // Reynolds number
    // drag coefficient
    C_D = 24. * (1. + 0.15 * pow(epsilon_g * num_Re, 0.687)) / (epsilon_g * num_Re + 1.e-8); 

    if (epsilon_g > 0.8)
    { // Wen-Yu
        k_g_s_HG = 0.75 * epsilon_g * epsilon_s * rho_g * u_slip * C_D * pow(epsilon_g, -2.65) / diam_p;
    }
    else
    { // Ergun
        k_g_s_HG = 150. * (mu_g * pow(epsilon_s, 2.) / (epsilon_g * pow(diam_p, 2.))) + 1.75 * u_slip * (rho_g * epsilon_s) / diam_p;
    }

    coef_Hd = 1.;
    k_g_s = coef_Hd * k_g_s_HG;

    // C_UDMI(cell, thread_mix, 4) = coef_Hd;

    return k_g_s;
}

// ----------------------------------------------------------
//         set velocity inlet of solids at solids inlet
// ----------------------------------------------------------
DEFINE_PROFILE(Us_at_solids_inlet, thread, i)
{
#if !RP_HOST
    Domain *domain;
    face_t face;
    cell_t cell;
    Thread *thread_mix, *thread_gas, *thread_solids;

    real flux = 0.0;
    real vfr_s_outlet = 0.0;

    domain = Get_Domain(1);
    thread_mix = Lookup_Thread(domain, ID_outlet);
    thread_gas = THREAD_SUB_THREAD(thread_mix, 0);
    thread_solids = THREAD_SUB_THREAD(thread_mix, 1);
#endif

#if !RP_HOST
    begin_f_loop(face, thread_solids)
    {
        flux += F_FLUX(face, thread_solids);
    }
    end_f_loop(face, thread_solids)

    flux = PRF_GRSUM1(flux);
    vfr_s_outlet = flux / density_solids;

    begin_f_loop(face, thread)
    {
        F_PROFILE(face, thread, i) = vfr_s_outlet / (A_inlet_solids * e_s_0);
    }
    end_f_loop(face, thread)
#endif
}